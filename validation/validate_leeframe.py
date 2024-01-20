from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import stats 
import os
import sys
# sys.path.append('.')
from arc_length.force_control_solver import force_control # import force control formulation of arc-length solver
from pathlib import Path
import argparse

# Dealing with ufl legacy
try:
    from ufl import diag, Jacobian, shape
except:
    from ufl_legacy import diag, Jacobian, shape

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 1
parameters['reorder_dofs_serial'] = False

ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

mesh = Mesh()
with XDMFFile('../examples/force_control/beam/beam_2D/mesh/lee_frame.xdmf') as infile: # os.getcwd()+'/examples/force_control/beam/beam_2D/mesh/lee_frame.xdmf'
    infile.read(mesh)

Ue = VectorElement("CG", mesh.ufl_cell(), 1, dim=2) # displacement
Te = FiniteElement("CG", mesh.ufl_cell(), 1) # rotation
V = FunctionSpace(mesh, MixedElement([Ue, Te]))   

v_ = TestFunction(V)
u_, theta_ = split(v_)
dv = TrialFunction(V)
v = Function(V, name="Generalized displacement")
u, theta = split(v)

VR = TensorFunctionSpace(mesh, "DG", 0, shape=(2, 2))

V0 = FunctionSpace(mesh, "DG", 0)


Vu = V.sub(0).collapse()
total_displ = Function(Vu, name="Total displacement")

Vr = V.sub(1).collapse()
total_rot = Function(Vr, name="Total Rotation")

Jac = Jacobian(mesh)
gdim = mesh.geometry().dim()
Jac = as_vector([Jac[i, 0] for i in range(gdim)])
g01 = Jac/sqrt(dot(Jac, Jac))
g02 = as_vector([-g01[1],g01[0]])

r01 = outer(g01,as_vector([1,0]))
r02 = outer(g02, as_vector([0,1]))

R0 = r01+r02

#-----------------------Define Functions for beams-----------------------------------#
def tgrad(u): # directional derivative w.r.t. beam centerline
    return dot(grad(u), g01)

def rotation_matrix(theta): # 2D rotation matrix -- there is no need to do rotation parametrization for 2D beams
    return as_tensor([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]])

Rot = rotation_matrix(theta)

# Beam Geometry
Lx = 120
Ly = 120
Lf = 24

#----------------------Dirichlet Boundary Conditions----------------------------------#
def Pinned_Bottom(x, on_boundary):
    return near(x[0], 0, 1e-6) and near(x[1], 0, 1e-6)
    
def Pinned_Top(x, on_boundary):
    return near(x[0], Lx, 1e-6) and near(x[1],Ly, 1e-6)

def load(x, on_boundary):
    return near(x[0],Lf, 1e-6)

facets = MeshFunction("size_t", mesh, 0)
facets.set_all(0)
AutoSubDomain(load).mark(facets,1)

Pinned_bot = DirichletBC(V.sub(0), Constant((0.0, 0.0)), Pinned_Bottom, method='pointwise') # displacement + rotation
Pinned_top = DirichletBC(V.sub(0), Constant((0.0, 0.0)), Pinned_Top, method='pointwise') # displacement + rotation
bcs = [Pinned_bot,Pinned_top]

#------------------------------2D beam Kinematics-----------------------------------------#
defo = dot(R0.T,dot(Rot.T, g01 + tgrad(u)) - g01)
curv =  tgrad(theta)

#------------------------------Beam Properties---------------------------------------#
# Geometric Properties
S = 6.0 # cross-sectional area
I = 2.0 # Area moment
E = 720.0 # Elastic modulus
nu = 0.3 # Poisson's Ratio

G = E/(2*(1+nu)) # Shear Modulus
kappa = 1.0 # Shear correction

# Beam Stiffness
ES = E*S
GS = G*S
GS_2 = G*S*kappa
GS_3 = G*S*kappa
EI = E*I

#------------------------Weak Form-----------------------------------------#
# Constitutive Equation
C_N = diag(as_vector([ES, GS_2]))

# Applied Load:
F_max = Constant((0.0,-1.0))
M_max = Constant(0.0)
load = Expression("t", t=0, degree = 0)

dS = Measure("dS", domain=mesh, subdomain_data=facets) # Note its dS not ds -- interior facet
dx = Measure("dx", domain=mesh)

elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx

F_int = derivative(elastic_energy, v, v_)
F_ext = avg(load * (-M_max*theta_ + dot(F_max, u_))) * dS(1) # avg function is used to integrte over disconticuous interior facets
residual = F_int - F_ext
tangent_form = derivative(residual, v, dv)

#-----------Get coordinated of force application for plotting---------------------------#
# Get dof coordinates:
x_dofs = V.sub(0).sub(0).dofmap().dofs()
y_dofs = V.sub(0).sub(1).dofmap().dofs()
theta_dofs = V.sub(1).dofmap().dofs()
dofs = V.tabulate_dof_coordinates()
dof_coords = dofs.reshape((-1, 2))

eps = 1e-5
dof = []
# Identify y-dof at the center for force application
for ii in y_dofs:
    if abs(dof_coords[ii,1]-Ly) <= eps:
        if abs(dof_coords[ii,0]-Lf) <=eps:
            force_dof = ii

#----------------------Solving the Nonlinear Problem--------------------------------------#
# Solver Parameters
psi = 1.0
abs_tol = 1.0e-6
lmbda0 = 0.5
max_iter = 10
solver = 'mumps' # optional

# Set up arc-length solver
solver = force_control(psi=psi, abs_tol= abs_tol, lmbda0=lmbda0, max_iter=max_iter, u=v,
                       F_int=F_int, F_ext=F_ext, bcs=bcs, J=tangent_form, load_factor=load, solver=solver)

# Lists to store solution vectors
disp = [v.vector()[:]]
lmbda = [0.0]
force_disp = [0.0] # displacement at force application

#---------------------------Compare Solution with literature solution--------------------------#
# Get solution from literature (https://www.sciencedirect.com/science/article/pii/S014102962034356X)

paper_eq = np.loadtxt('../examples/force_control/beam/beam_2D/lit_soln/path_lee.dat')
paper_disp = -paper_eq[:,2]
paper_load = paper_eq[:,0]

#Setup command line arguments
parser = argparse.ArgumentParser(description='Optional paraview saving')
parser.add_argument('-p', '--paraview', 
                    default = False, required = False, action = 'store_true',
                    help = 'True/False to save paraview files') 
args = parser.parse_args()

# Setup directory and save file if needed
if args.paraview:
    # Create directory if it doesn't exist
    Path("paraview").mkdir(parents=True, exist_ok=True)
    # Initialize file
    out_file = XDMFFile('paraview/validate_leesframe.xdmf')
    out_file.parameters['functions_share_mesh'] = True
    out_file.parameters['flush_output'] = True

ii=0
# solver iteration
while (-v.vector()[force_dof] <= paper_disp[-1]) or (solver.converged == False):
    solver.solve()
    if solver.converged and (-v.vector()[force_dof] <= paper_disp[-1]): # We only want to save the solution step if the solver covergesand within range of paper solution
        total_displ.vector()[:] = v.sub(0,True).vector()[:]
        total_rot.vector()[:] = v.sub(1,True).vector()[:]
        disp.append(v.vector()[:])
        force_disp.append(-disp[-1][force_dof])
        lmbda.append(load.t)
        ii+=1

        if args.paraview:
            out_file.write(total_displ, ii)
            out_file.write(total_rot, ii)
if args.paraview:
    out_file.close()
# Plot and compare solutions
plt.figure(figsize=(7,7))
plt.scatter(-paper_eq[:,2], paper_eq[:,0], label = 'Paper implementation', facecolors = 'None', edgecolors = 'r', marker = 's')

plt.scatter(force_disp, lmbda, marker = '.', label = 'Our Implementation')
plt.title('Equilibrium Path')
plt.ylabel('Load')
plt.xlabel('Displacement')
plt.legend()

#---------------Compare critical buckling load----------------------------------------------------------#
our_critical = np.max(lmbda)
lit_critical = np.max(paper_eq[:,0])

diff_critical = np.abs((lit_critical-our_critical)/lit_critical)
print(f'\n\nPercent Difference in Critical Buckling Load: {diff_critical*100:.4f} %')
if diff_critical*100 < 0.1:
    val = [True]
else:
    val = [False]

#---------------Create parameterized splines to compare the curves--------------------------------------#
t = np.linspace(0,1,200)

lit_t = np.linspace(0,1,len(-paper_eq[:,2]))
lit_y = np.vstack((-paper_eq[:,2],paper_eq[:,0]))

lit_spline = interp1d(lit_t,lit_y)

our_t = np.linspace(0,1,len(force_disp))
our_y = np.vstack((np.array(force_disp),lmbda))

our_spline = interp1d(our_t,our_y)



lit= lit_spline(t)
our = our_spline(t)


corr = stats.pearsonr(our[1,:],lit[1,:])

if 1-corr[0] < 0.1:
    val.append(True)
else:
    val.append(False)

print(f'Comparing Pearson correlation between force-displacement curves: {corr[0]:.4f}')

if np.all(val):
    print('Lee\'s frame validation complete!')
else:
    print('Lee\'s frame did not pass validation test')
    print('Passed test:', val)

# Create directory if it doesn't exist
Path("plots").mkdir(parents=True, exist_ok=True)

plt.savefig('plots/validate_leeframe.png')