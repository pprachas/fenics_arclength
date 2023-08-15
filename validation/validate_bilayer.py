from dolfin import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from scipy import optimize
from arc_length.displacement_control_solver import displacement_control # import displacement control formulation of arc-length solver
import sys
from pathlib import Path
import argparse

init_time = time.time() # use to time code
# FEniCS solver parameters
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

# Stiffness ratio
SR = 100

# Elasticity parameters (Plane strain assumptions)
E, nu = 10.0, 0.45

# Lame Constants
mu_f, lmbda_f = Constant(SR*E/(2*(1 + nu))), Constant(SR*E*nu/((1 + nu)*(1 - 2*nu)))
mu_s, lmbda_s = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))
# Define mesh based on analytical wavelength
Hf = 1.0 # film height
ana_wavelength = float(2*np.pi*Hf*((mu_f)/(3*mu_s))**(1/3))
Hs = ana_wavelength*4.0
L = ana_wavelength*3 # length of the domain
SR = 100 # stiffness ratio
d = Hf/2 # mesh size

nx = int(L/d)
ny = int((Hs+Hf)/d)
mesh = RectangleMesh(Point(0,0), Point(L,Hs+Hf), nx, ny, "right")
V = VectorFunctionSpace(mesh, "Lagrange", 2)
plt.figure(figsize=(7,7))
plot(mesh)

# Define Variational Form
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u = Function(V)                 # Solution

# Mark boundary subdomians
domain_markers = MeshFunction("size_t", mesh, mesh.topology().dim())
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1)

def Film(x, on_boundary):
    return x[1] >= Hs
    
def Substrate(x, on_boundary):
    return x[1] < Hs
    

AutoSubDomain(Film).mark(domain_markers,0)
AutoSubDomain(Substrate).mark(domain_markers,1)

# Dirichlet boundary conditions
def Left(x, on_boundary):
    return on_boundary and near(x[0], 0, 1e-6)

def Right(x, on_boundary):
    return on_boundary and near(x[0], L, 1e-6 )

def Bottom(x, on_boundary):
        return on_boundary and near(x[1], 0, 1e-6)

#-----------------Applied Displacement-----------------------#
apply_disp = Expression("t",t = 0.0, degree = 0) 
#---------Nonhomogenous Dirichlet Boundary Conditions---------#
bc1 = DirichletBC(V.sub(0), Constant(0), Left)
bc2 = DirichletBC(V.sub(0), apply_disp, Right)
bc3 = DirichletBC(V.sub(1), Constant(0), Bottom)
bcs = [bc1, bc2, bc3]

# Define function for point load (the perturbation)
class PointLoad(UserExpression):
    def __init__(self, x0, f, tol,**kwargs):
        super().__init__(**kwargs)
        self.x0 = x0
        self.f = f
        self.tol = tol
    def eval(self, values, x):
        if near (x[0], self.x0[0],self.tol) and near(x[1], self.x0[1],self.tol):
            values[0] = self.f[0]
            values[1] = self.f[1]
        else:
            values[0] = 0
            values[1] = 0
    def value_shape(self):
        return (2,)

# Kinematics
I = variable(Identity(mesh.topology().dim()))  # Identity tensor
F = variable(I + grad(u))                        # Deformation gradient
C = variable(F.T*F)                              # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Define Variational Form
dx = Measure('dx', domain=mesh, subdomain_data=domain_markers)

# Stored strain energy density (compressible neo-Hookean model)
psi_f = (mu_f/2)*(Ic - 2) - mu_f*ln(J) + (lmbda_f/2)*(ln(J))**2
psi_s = (mu_s/2)*(Ic - 2) - mu_s*ln(J) + (lmbda_s/2)*(ln(J))**2

# Applied Traction and Body Force
pert = 8e-2
T = Constant((0,0)) # Traction
B = Constant((0,0)) # Body Force
P = PointLoad(x0 = [L/2,Hf+Hs], f = [0,pert], tol = 1e-4, degree = 0, element=V.ufl_element()) # Point load

F_int = derivative(psi_f*dx(0) + psi_s*dx(1), u, v)
F_ext = derivative(dot(B, u)*dx + dot(T,u)*ds + dot(-P,u)*ds , u, v)
residual = F_int-F_ext
J = derivative(residual, u, du)

# Solver Parameters
psi = 1.0
abs_tol = 1.0e-6
lmbda0 = -0.07
max_iter = 30

# Set up arc-length solver
solver = displacement_control(psi=psi, abs_tol=abs_tol, lmbda0=lmbda0, max_iter=max_iter, u=u,
                       F_int=F_int, F_ext=F_ext, bcs=bcs, J=J, displacement_factor=apply_disp)

disp = [u.vector().copy()]
lmbda = [0]

# Function space to compute reaction force at each iteration
v_reac = Function(V)
bcRx = DirichletBC(V.sub(0), Constant(1.0), Left) # take reaction force from the left boundary
f_reac = [0.0]

strain_crit = float((1/4)*((3*mu_s/mu_f))**(2/3))# critial strain

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
    out_file = XDMFFile('paraview/validate_bilayer.xdmf')
    out_file.parameters['functions_share_mesh'] = True
    out_file.parameters['flush_output'] = True

ii = 0
while np.abs(apply_disp.t) < strain_crit*L*1.1 and solver.converged:
    solver.solve()
    if solver.converged:
        # Store whole displacement field
        disp.append(u.vector().copy())
        # Store displacement factor
        lmbda.append(apply_disp.t)
        # Compute and store reaction force
        bcRx.apply(v_reac.vector())
        f_reac.append(assemble(action(residual,v_reac)))

        ii+=1

        if args.paraview:
            out_file.write(u, ii)

if args.paraview:
    out_file.close()
# Post Processing (validate with analytical solution)
# Here we plot and compare the final deformed shape, equilibrium path, and the wavelength. We obtain the analytical solutions from: https://royalsocietypublishing.org/doi/epdf/10.1098/rsta.2016.0163 and https://groups.seas.harvard.edu/hutchinson/papers/WrinklingPhenonmena-JAMF.pdf

test_y = np.diff(np.array(f_reac))

test_x = np.diff(-np.array(lmbda)/L)

fea_soln = np.nanargmax(np.abs((np.diff(np.diff(((np.array(f_reac[1:])/(Hs+Hf))/(-np.array(lmbda[1:]))/L)))))) # find the inflection point from FEA solution and use that as critical strain

# Create directory if it doesn't exist
Path("plots").mkdir(parents=True, exist_ok=True)

# Plot comparison with analytical solutions
plt.figure(figsize=(7,5))
plt.plot(-np.array(lmbda)/L, np.array(f_reac)/(Hs+Hf), c='k', marker = 'o', label = 'Equilibrium Path')
plt.plot(-np.array(lmbda[fea_soln+2])/L, np.array(f_reac[fea_soln+2])/(Hs+Hf), marker = 'd', c =(0.4,0,0) , ls = 'None', markersize = 10, label='FEA critical strain')
plt.axvline(strain_crit, ls = ':', c = (0.7,0,0), label = 'Analytical Critical Strain')
plt.xlabel('Applied Strain')
plt.ylabel('Stress per unit depth')
plt.title('Equilibrium path')
plt.legend()

plt.savefig('plots/validate_bilayer_stressstrain.png')

percent_diff_crit_strain = ((strain_crit-(-np.array(lmbda[fea_soln+2])/L))/strain_crit) * 100
print('Percent Error between analytical critical strain and FEA critical strain:',percent_diff_crit_strain,'%')

if percent_diff_crit_strain < 2.0:
    val = [True]
else:
    val= [False]
# Here we extract and compare the post-bifurcation wrinkling wavelength from FEA and analytical solution.

post_bif = np.argwhere(-np.array(lmbda)/L > strain_crit).reshape(-1) # get index after bifurcation
# get solution after onset of bifurcation:
disp_bif = disp[post_bif[0]]
u_bif = Function(V) 
u_bif.vector()[:] = disp_bif[:]


# Get dof coordinates:
x_dofs = V.sub(0).dofmap().dofs()
y_dofs = V.sub(1).dofmap().dofs()
theta_dofs = V.sub(1).dofmap().dofs()
dofs = V.tabulate_dof_coordinates()
dof_coords = dofs.reshape((-1, 2))

x_nodal_coord = dof_coords[x_dofs][:,0]
y_nodal_coord = dof_coords[y_dofs][:,1]

# Get nodal values of the top layer 
top_layer = np.argwhere(np.abs(y_nodal_coord-(Hs+Hf)) < 1e-6)

# Plot displacement field
disp_x = x_nodal_coord[top_layer].reshape(-1)#(x_nodal_coord + u_bif.vector()[x_dofs])[top_layer].reshape(-1)
disp_y = (y_nodal_coord + u_bif.vector()[y_dofs])[top_layer].reshape(-1)


def fit_func(x, amp, wave, phi, offset):
    return amp * np.cos((2*np.pi)/wave * x + phi) + offset

params, params_covariance = optimize.curve_fit(fit_func, disp_x, disp_y,
                                               p0=[1.0,ana_wavelength,0.0,Hf+Hs])

plt.figure(figsize = (7,7))

plt.title('FEA wavelength vs Analytical Wavelength')
plt.scatter(disp_x,disp_y, marker = '.', c = 'r', s = 50, label='FEA data')
plt.scatter(disp_x,fit_func(disp_x, params[0], params[1], params[2], params[3]), c = 'None', marker = 's',edgecolor = 'k', s = 50, linewidth = 2, label = 'Fitted FEA data')
plt.scatter(disp_x,fit_func(disp_x, params[0], ana_wavelength, params[2], params[3]), c = 'None', marker = 'd', edgecolor = (0.5,0,0), s = 50, linewidth = 2, label = 'Analytical Wavelength')
plt.xlabel('x displacement')
plt.ylabel('y displacement')
plt.legend(loc = (1.01,0.5))
plt.tight_layout()
plt.savefig('plots/validate_bilayer_wavelength.png')

percent_diff_wavelength=((params[1]-ana_wavelength)/ana_wavelength)*100
print('Percent Error between analytical wavelength and fitted FEM wavelength:',percent_diff_wavelength, '%')

if percent_diff_wavelength < 2.0:
    val.append(True)
else:
    val.append(False)

val = np.array([percent_diff_crit_strain, percent_diff_wavelength])
if np.all(val):
    print('Bilayer wrinkling validation complete!')
else:
    print('Bilayer wrinkling did not pass validation test')
    print('Passed test:', val)

print('Elapse time:', (time.time()-init_time)/60, 'minutes')