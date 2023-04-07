# %% [markdown]
# # Bilayer Wrinkling
# 
# The problem provided in this example is a bilayer wrinkling problem. The domain is a thin stiff film resting on top of a substrate. A small point-force perturbation is applied to the substrate to instigate the instability. The domain is compressed well after the onset of instability. 
# 
# <img src="imgs/bilayer.png" width="700">
# 

# %%
%matplotlib inline
from dolfin import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import os
os.chdir('../..') # change directory to main module directory
from arc_length.displacement_control_solver import displacement_control # import displacement control formulation of arc-length solver


parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

# %% [markdown]
# ## Define domain and function spaces
# 
# We also split the domain into substrate and film. We To make sure all the assumptions of a bilayer system are valid, we construct the doamin width and depth based on the analytical wavelength.

# %%
# Stiffness ratio
SR = 100

# Elasticity parameters (Plane strain assumptions)
E, nu = 10.0, 0.45
# Lame Constants
mu_f, lmbda_f = Constant(SR*E/(2*(1 + nu))), Constant(SR*E*nu/((1 + nu)*(1 - 2*nu)))
mu_s, lmbda_s = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))
# define mesh based on analytical wavelength
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

# %%
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

# Check Domain
plt.figure(figsize=(7,7))
plt.title('Domain')
plot(domain_markers)


# %% [markdown]
# ## Define Dirichlet Boundary Conditions
# 
# **Note that for the case of displacement control, the FEniCS expression for the applied displacement mujst have be positive to prevent convergence issues.**
# 
# For example:
# 
# ```
# apply_disp = Expression("t", t = 0.0, degree = 0)
# ```
# is valid and will not have convergence issues while
# ```
# apply_disp = Expression("-t", t = 0.0, degree = 0)
# ```
# can cause convergence issues.
# 
# The direction of applied loading will be determined by the initial load step.
# 
# 

# %%
# Dirichlet boundary conditions
def Left(x, on_boundary):
    return on_boundary and near(x[0], 0, 1e-6)

def Right(x, on_boundary):
    return on_boundary and near(x[0], L, 1e-6 )

def Bottom(x, on_boundary):
        return on_boundary and near(x[1], 0, 1e-6)

#-----------------Applied Displacement-----------------------#
apply_disp = Expression("t",t = 0.0, degree = 0) 
#--------------Nonhomogenous Dirichlet Boundary Conditions------------------------------#
bc1 = DirichletBC(V.sub(0), Constant(0), Left)
bc2 = DirichletBC(V.sub(0), apply_disp, Right)
bc3 = DirichletBC(V.sub(1), Constant(0), Bottom)
bcs = [bc1, bc2, bc3]

# %% [markdown]
# ## Define Function for Point Load
# Note that this is an approximation of a point load since the FEniCS UserExpression will be projected into a discontinuous space. To approach a point load the area near the point load will need to have fine mesh. In our case it is not neccessary since the point load is just a perturbation.

# %%
# Function for Point 
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

# %% [markdown]
# ## Kinematics and Weak form

# %%
# Kinematics
I = variable(Identity(mesh.topology().dim()))  # Identity tensor
F = variable(I + grad(u))                        # Deformation gradient
C = variable(F.T*F)                              # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# %%
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

# %% [markdown]
# ## Solver
# To use our solver we first have to define the type of solver (i.e. displacement control or force control) and solver parameters before using the solver. Note that the correct type of solver has to first be imported (see first cell).
# ### Solver parameters
# Here the parameters for both types of solvers:
# 
# >* `psi` : the scalar arc-legth parameter. When `psi` = 1, the method becomes the shperical arc-length method and when `psi` = 0 the method becomes the cylindrical arc-length method
# >* `tol` : tolerance for the linear solver
# >* `lmbda0` : the initial load parameter
# >* `max_iter` : maximum number of iterations for the linear solver
# >* `solver` : (optional): type of linear solver for the FEniCS linear solve function -- default FEniCS linear solver is used if no argument is used.
# 
# Aside from these solver parameters, the arguments need to solve the FEA problem must also be passed into the solver:
# >* `u` : the solution function
# >* `F_int` : First variation of strain energy (internal nodal forces)
# >* `F_ext` : Externally applied load (external applied force)
# >* `J` : The Jacobian of the residual with respect to the deformation (tangential stiffness matrix)
# >* `displacement_factor` : The incremental load factor 
# 
# The solver can be called by:
# `solver = force_control(psi,tol,lmbda0,max_iter,u,F_int,F_ext,bcs,J,load_factor,solver)`
# 
# ### Using the solver
# 1. Initialize the solver by calling solver.initialize()
# 2. Iteratively call solver.solve() until desired stopping condition

# %%
# Solver Parameters
psi = 1.0
tol = 1.0e-6
lmbda0 = -0.07
max_iter = 30

# Set up arc-length solver
solver = displacement_control(psi=psi, tol=tol, lmbda0=lmbda0, max_iter=max_iter, u=u,
                       F_int=F_int, F_ext=F_ext, bcs=bcs, J=J, displacement_factor=apply_disp)

# %%
disp = [u.vector().copy()]
lmbda = [0]
# Function space to compute reaction force at each iteration
v_reac = Function(V)
bcRx = DirichletBC(V.sub(0), Constant(1.0), Left) # take reaction force from the left boundary
f_reac = [0.0]

strain_crit = float((1/3)*((3*mu_s/mu_f))**(2/3))# critial strain
while np.abs(apply_disp.t) < strain_crit*L*1.2 and solver.converged:
    solver.solve()
    if solver.converged:
        # Store whole displacement field
        disp.append(u.vector().copy())
        # Store displacement factor
        lmbda.append(apply_disp.t)
        # Compute and store reaction force
        bcRx.apply(v_reac.vector())
        f_reac.append(assemble(action(residual,v_reac)))

# %% [markdown]
# ## Optional: Stop at at specific criterion
# To strongly impose the stop condition, we add Newton solver at the end with a specific prescribed displacement and replace the final solution with the solution from the Newton Solver. In this case, we enforced the stop criterion with the applied displacement is slightly over the analytical critical strain.

# %%
apply_disp.t = -strain_crit*L*1.2 # stop simulation after bifurcation
problem = NonlinearVariationalProblem(residual, u, bcs, J)
newton_solver = NonlinearVariationalSolver(problem)
newton_solver.solve()

disp.pop()
lmbda.pop()
f_reac.pop()
disp.append(u.vector())
lmbda.append(apply_disp.t)
f_reac.append(assemble(action(residual,v_reac)))

# %% [markdown]
# ## Post Processing
# Here we plot the final deformed shape, equilibrium path, and the wavelength. We obtain the analytical solutions from: https://royalsocietypublishing.org/doi/epdf/10.1098/rsta.2016.0163 and https://groups.seas.harvard.edu/hutchinson/papers/WrinklingPhenonmena-JAMF.pdf

# %%
plt.figure(figsize=(7,7))
plot(u, mode = 'displacement', cmap='Reds', edgecolor = 'k', linewidth = 0.2)
plt.title('Final Deformation')

# %%
strain_crit = float((1/4)*((3*mu_s/mu_f))**(2/3))

test_y = np.diff(np.array(f_reac))

test_x = np.diff(-np.array(lmbda)/L)

fea_soln = np.nanargmax(np.abs((np.diff(np.diff(((np.array(f_reac[1:])/(Hs+Hf))/(-np.array(lmbda[1:]))/L)))))) # find the inflection point from FEA solution and use that as critical strain

plt.figure(figsize=(7,5))
plt.plot(-np.array(lmbda)/L, np.array(f_reac)/(Hs+Hf), c='k', marker = 'o', label = 'Equilibrium Path')
plt.plot(-np.array(lmbda[fea_soln+2])/L, np.array(f_reac[fea_soln+2])/(Hs+Hf), marker = 'd', c =(0.4,0,0) , ls = 'None', markersize = 10, label='FEA critical strain')
plt.axvline(strain_crit, ls = ':', c = (0.7,0,0), label = 'Analytical Critical Strain')
plt.xlabel('Applied Strain')
plt.ylabel('Stress per unit depth')
plt.title('Equilibrium path')
plt.legend()

percent_diff_crit_strain = ((strain_crit-(-np.array(lmbda[fea_soln+2])/L))/strain_crit) * 100
print('Percent Error between analytical critical strain and FEA critical strain:',percent_diff_crit_strain,'%')

# %% [markdown]
# Here we extract and compare the post-bifurcation wrinkling wavelength from FEA and analytical solution.

# %%
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
# Get nodal values 

top_layer = np.argwhere(np.abs(y_nodal_coord-(Hs+Hf)) < 1e-6)

# Plot displacement field
disp_x = x_nodal_coord[top_layer].reshape(-1)#(x_nodal_coord + u_bif.vector()[x_dofs])[top_layer].reshape(-1)
disp_y = (y_nodal_coord + u_bif.vector()[y_dofs])[top_layer].reshape(-1)

# Fit to sine functions
from scipy import optimize

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

print('Percent Error between analytical wavelength and fitted FEM wavelength:',((params[1]-ana_wavelength)/ana_wavelength)*100, '%')


# %% [markdown]
# ## Optional: Creating an animation from solution snapshots

# %% [markdown]
# 

# %%
from matplotlib import animation, rc

plt.rcParams["animation.html"] = "jshtml"

u_plot = Function(V)

fig = plt.figure(figsize=(12,5))
ax = fig.gca()

ax.set_xlim([0,L])
ax.set_ylim([(Hf+Hs)/3,Hs+Hf+5])

def drawframe(n):
    fig.clf()
    ax = fig.gca()
    ax.set_xlim([0,L])
    ax.set_ylim([(Hf+Hs)/3,Hf+Hs+5])
    u_plot = Function(V)
    u_plot.vector()[:] = disp[n][:]
    return plot(u_plot, mode = 'displacement', cmap = 'Reds'),
# blit=True re-draws only the parts that have changed.
anim = animation.FuncAnimation(fig, drawframe, frames=len(lmbda), interval=40, blit=True)
anim

# %%



