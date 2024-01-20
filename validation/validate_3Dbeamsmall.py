from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import sys
from arc_length.rotation_parametrization import ExponentialMap # import rotation parameterization for 3D beams
from arc_length.force_control_solver import force_control # import force control formulation of arc-length solver
from pathlib import Path
import argparse

# Dealing with ufl legacy
try:
    from ufl import diag, Jacobian, shape
except:
    from ufl_legacy import diag, Jacobian, shape

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 3
parameters['reorder_dofs_serial'] = False

mesh = Mesh()
with XDMFFile('../examples/force_control/beam/beam_3D/mesh/mesh.xdmf') as infile:
    infile.read(mesh)

Ue = VectorElement("CG", mesh.ufl_cell(), 2, dim=3) # displacement
Te = VectorElement("CG", mesh.ufl_cell(), 1, dim=3) # rotation
V = FunctionSpace(mesh, MixedElement([Ue, Te]))   

v_ = TestFunction(V)
u_, theta_ = split(v_)
dv = TrialFunction(V)
v = Function(V, name="Generalized displacement")
u, theta = split(v)

VR = TensorFunctionSpace(mesh, "DG", 0, shape=(3, 3))
R_old = Function(VR, name="Previous rotation matrix")
R_old.interpolate(Constant(((1, 0, 0), (0, 1, 0), (0, 0, 1))))

V0 = VectorFunctionSpace(mesh, "DG", 0, dim=3)
curv_old = Function(V0, name="Previous curvature strain")


Vu = V.sub(0).collapse()
total_displ = Function(Vu, name="Total displacement")

Vr = V.sub(1).collapse()
total_rot = Function(Vr, name="Total Rotation")

rot_param = ExponentialMap()
R = rot_param.rotation_matrix(theta)
H = rot_param.curvature_matrix(theta)

Jac = Jacobian(mesh)
gdim = mesh.geometry().dim()
Jac = as_vector([Jac[i, 0] for i in range(gdim)])
g01 = Jac/sqrt(dot(Jac, Jac))
    
g02 = cross(as_vector([0,0,1]),g01)
g02 /= sqrt(dot(g02,g02))
    
g03 = cross(g01,g02)
g03 /= sqrt(dot(g03,g03))

r_01 = outer(g01,as_vector([1,0,0]))
r_02 = outer(g02,as_vector([0,1,0]))
r_03 = outer(g03,as_vector([0,0,1]))

R0=r_01+r_02+r_03 # Initial rotation tensor from global coordinate frame to beam orthonormal triad

#-----------------------------------Define Directional Derivative----------------------------------#
def tgrad(u):
    return dot(grad(u), g01)

L = 10
def left_end(x, on_boundary):
    return near(x[0], 0)
def right_end(x, on_boundary):
    return near(x[0], L)

facets = MeshFunction("size_t", mesh, 0)
AutoSubDomain(left_end).mark(facets, 1)
AutoSubDomain(right_end).mark(facets, 2)

#-------------------------Clamp Left----------------------------------------------------------------#
bcs = [DirichletBC(V, Constant((0,)*6),left_end, method = 'pointwise')]

#-------------------Use Incremental Method----------------------------------------------------------#
method = "total"

if method == "total":
    defo = dot(R0,dot(R.T, g01 + tgrad(u)) - g01)
    curv = dot(H.T, tgrad(theta))
elif method == "incremental":
    R_new = R * R_old
    defo = dot(R0.T,dot(R_new.T, g01 + tgrad(total_displ+u)) - g01)
    curv = curv_old + dot(R_old.T * H.T, tgrad(theta))

#---------------Beam Properties------------------------------#
# Geometric Properties
radius = Constant(0.2)
S = pi * radius ** 2
I = pi * radius ** 4 / 4

# Beam Stiffness
ES = Constant(1e4)
GS = Constant(1e4)
GS_2 = GS
GS_3 = GS
EI = Constant(1e2)
EI_2 = EI
EI_3 = EI
GJ = Constant(1e2)

#----------------------Weak Form----------------------------#
# Constitutive Equation
C_N = diag(as_vector([ES, GS_2, GS_3]))
C_M = diag(as_vector([GJ, EI_2, EI_3]))

# Applied Loads:
M_max = Constant(1.0) #Constant(200*pi)
F_max = Constant(1.0)
load = Expression("t", t=0, degree=0)

ds = Measure("ds", domain=mesh, subdomain_data=facets)
dx = Measure("dx", domain=mesh)

elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + dot(curv, dot(C_M, curv))) * dx

F_int = derivative(elastic_energy, v, v_)
F_ext = load * (-M_max* dot(H, theta_)[1] + F_max * u_[1]) * ds(2)
residual = F_int - F_ext
tangent_form = derivative(residual, v, dv)

# Solver Parameters
psi = 1.0 
abs_tol = 1.0e-6 
lmbda0 = 1e-3/2
max_iter = 10
solver = 'mumps' # Optional -- Type of linear solver

# Set up arc-length solver
solver = force_control(psi=psi, abs_tol=abs_tol, lmbda0=lmbda0, max_iter=max_iter, u=v, 
                       F_int=F_int, F_ext=F_ext, bcs=bcs, J=tangent_form, load_factor=load, solver=solver)

# Start solver 
disp = [total_displ.vector()[:]] # displacement solutions
rot = [v.sub(1,True).vector()[:]] # Rotation solutions
lmbda = [0] # save load factor

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
    out_file = XDMFFile('paraview/validate_3Dbeamsmall.xdmf')
    out_file.parameters['functions_share_mesh'] = True
    out_file.parameters['flush_output'] = True

ii=0
while (solver.load_factor.t < 1e-3 and solver.converged) or (solver.converged == False):
    solver.solve()
    if solver.converged: 
        ii+=1      
        if method == 'incremental':
            total_displ.vector()[:] += v.sub(0,True).vector()[:]
            total_rot.vector()[:] += v.sub(1,True).vector()[:]
            R_old.assign(project(R * R_old, VR))
            curv_old.assign(project(curv, V0))

            disp.append(total_displ.vector()[:])
            rot.append(v.sub(1,True).vector()[:])
            lmbda.append(solver.load_factor.t)

            solver.u_n = v.vector().copy()-solver.u_n_1
            solver.u_n_1.zero()

            if args.paraview:
                out_file.write(total_displ, ii)
                out_file.write(total_rot, ii)
        else:
            total_displ.vector()[:] = v.sub(0,True).vector()[:]
            total_rot.vector()[:] = v.sub(1,True).vector()[:]
            disp.append(v.sub(0,True).vector()[:])
            rot.append(v.sub(1,True).vector()[:])

            if args.paraview:
                out_file.write(total_displ, ii)
                out_file.write(total_rot, ii)

if args.paraview:
    out_file.close()
#----------Comparing with analytical curvature (e_2 direction)------------------#
# get curvature from FEA solution
curve_e2 = project(curv, V0)(0,0,0)[1]

# Compute analytical solution from moment-curvature relationship: \kappa = M/EI
applied_moment = -float(M_max)*load.t
curve_ana = applied_moment/float(EI)

percent_difference = np.abs((curve_e2-curve_ana)/curve_ana)*100

print('\n\nComparing analytical and FEA beam curvature (e2 direction)')
print(f'FEA solution: {curve_e2:.4e}')
print(f'Analytical Solution: {curve_ana:.4e}')
print(f'Percent difference: {percent_difference: .4e}')

if percent_difference < 1.0:
    val = [True]
else:
    val = [False]

#-------------Comparing analytical shear (e_2)-----------------------------------------#
v_reac = Function(V)
bcRy = DirichletBC(V.sub(0).sub(1), Constant(1.), left_end)
bcRy.apply(v_reac.vector())
shear_fea = assemble(action(residual, v_reac))

shear_ana = -load.t*float(F_max)

shear_percent_difference = np.abs(((shear_ana-shear_fea)/shear_ana)*100)

print('\n\nComparing analytical and FEA reaction shear (e2 direction)')
print(f'FEA solution: {shear_fea:.4e}')
print(f'Analytical Solution: {shear_ana:.4e}')
print(f'Percent difference: {shear_percent_difference: .4e}')

if shear_percent_difference < 1.0:
    val.append(True)
else:
    val.append(False)

#-------------Comparing analytical shear (e_3)-----------------------------------------#
v_reac = Function(V)
bcRz = DirichletBC(V.sub(1).sub(2), Constant(1.), left_end)
bcRz.apply(v_reac.vector())
moment_fea = assemble(action(residual, v_reac))

moment_ana=-load.t*float(F_max)*L

moment_percent_difference = np.abs(((moment_ana-moment_fea)/shear_ana)*100)

print('\n\nComparing analytical and FEA reaction moment (e3 direction)')
print(f'FEA solution: {moment_fea:.4e}')
print(f'Analytical Solution: {moment_ana:.4e}')
print(f'Percent difference: {moment_percent_difference: .4e}')

if moment_percent_difference < 1.0:
    val.append(True)
else:
    val.append(False)

if np.all(val):
    print('3D beam small deformation validation complete!')
else:
    print('3D beam small deformation did not pass validation test')
    print('Passed test:', val)
