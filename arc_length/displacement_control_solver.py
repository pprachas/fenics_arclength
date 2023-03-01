from dolfin import *
import numpy as np
from math import sqrt

'''
An FEniCS based arc-length with prescribed non-zero displacement. The predictor-corrector step is heavily based on the paper:

Kadapa, Chennakesava. "A simple extrapolated predictor for overcoming the starting and tracking issues in the arc-length method for nonlinear structural mechanics." Engineering Structures 234 (2021): 111755.
 
and the displacement-based formulation is heavily based on:
Verhoosel, Clemens V., Joris JC Remmers, and Miguel A. Gutiérrez. "A dissipation‐based arc‐length method for robust simulation of brittle and ductile failure." International journal for numerical methods in engineering 77.9 (2009): 1290-1321.

'''

class displacement_control: 
    ''' The arc-length displacement control solver of this library
    
    Args:
        psi: the scalar arc-length parameter. When psi = 1, the method becomes the spherical arc-length method and when psi = 0 the method becomes the cylindrical arc-length method
        tol : tolerance for the linear solver
        lmbda0 : the initial load parameter
        max_iter : maximum number of iterations for the linear solver
        u : the solution function
        F_int : First variation of strain energy (internal nodal forces)
        F_ext : Externally applied load (external applied force)
        J : The Jacobian of the residual with respect to the deformation (tangential stiffness matrix)
        displacement_factor : The incremental load factor
        solver : (optional): type of linear solver for the FEniCS linear solve function -- default FEniCS linear solver is used if no argument is used.
    '''

    def __init__(self, psi, tol, lmbda0, max_iter, u, F_int, F_ext, bcs, J, displacement_factor, solver='default'):
        # Initialize Variables
        self.psi = psi
        self.tol = tol
        self.lmbda = lmbda0
        self.max_iter = max_iter
        self.F_int = F_int
        self.F_ext = F_ext
        self.u = u
        self.J = J
        self.bcs = bcs
        self.displacement_factor = displacement_factor
        self.residual = F_int-F_ext
        self.solver = solver
        self.counter = 0
        self.converged = True
    
    def __update_nodal_values(self, u_new):
        '''
        Function to update solution (i.e. displacement) vector after each solver iteration
        
        Args:
            u_new: updated solution

        '''

        # Function to update displacements
        u_nodal_values = u_new.get_local()
        self.u.vector().set_local(u_nodal_values)
    
    def __initial_step(self):
        '''
        Inital step of the arc-length method. 
        For the displacement control formulation, this function constructs the constraint matrix and the initial arc-length step size.
        '''    

        ii=0
        print('Starting initial Displacement Control Control with Newton Method:')
        
        # Find DoFs for homogenous and nonhomogenous Dirichlet BCs:
        self.displacement_factor.t = 1.0 #set nonhomogenous DoFs to 1 
        
        self.dofs_hom = [] # list of DoFs with homogenous DoFs
        self.dofs_nonhom = [] # list of DoFs with nonhomogenous DoFs
        for bc in self.bcs:
            for (dof, value) in bc.get_boundary_values().items():
                if value == 0:
                    self.dofs_hom.append(dof)
                else:
                    self.dofs_nonhom.append(dof)
        
        self.dofs_free = np.arange(0,len(self.u.vector()))
                              
        self.dofs_free = np.delete(self.dofs_free,(self.dofs_hom+self.dofs_nonhom)) # vector of free DoF
        
        # set up DoF vector of Dirchlet BCs:
        self.u_p = Vector()
        self.u_p.init(self.u.vector().size())
        u_p_temp = [0]*self.u.vector().size()
        
        for dof in self.dofs_nonhom:
            u_p_temp[dof] = 1
        
        self.u_p.set_local(u_p_temp)
        
                                  
        
        # Construct Constraint Matrix C
        C = PETScMatrix()
        C_mat = C.mat() # use petsc interface
        C_mat.create()
        C_mat.setSizes((self.u.vector().size(), len(self.dofs_free)))
        C_mat.setType('aij')
        C_mat.setUp()
        for j, i in enumerate(self.dofs_free):
            C_mat.setValue(i, j, 1)
        C_mat.assemble()
        self.C = PETScMatrix(C_mat)

        self.u_f = Vector()
        self.C.init_vector(self.u_f, 1)
        
        # Initialize reduced residual vector
        R_star = Vector()
        self.C.init_vector(R_star, 1)
        
        # Initialize Q
        self.Q = Vector()
        self.C.init_vector(self.Q, 1)
        
        du = Vector()
        self.C.init_vector(du, 0)
        self.C_mat = self.C.mat()

        # Apply all Dirichlet (both homogenous and non-homogenous) BCs to u:
        self.displacement_factor.t = self.lmbda
        for bc in self.bcs:
            bc.apply(self.u.vector())

        while True:
            # Assemble K and R
            K = assemble(self.J)
            R = assemble(self.residual)

            # Get K*, F*, and R* (reduced matrices and vectors)
            K_mat = as_backend_type(K).mat()
            temp = self.C_mat.copy().transpose().matMult(K_mat)
            K_star_mat = temp.matMult(self.C_mat) # Reduced stiffness matrix
            K_star = Matrix(PETScMatrix(K_star_mat))
            PETScMatrix(temp).mult(-self.u_p, self.Q) # vector of Dirichlet BCs
            self.C.transpmult(R, R_star) # reduced residual vector


            norm = R_star.norm('l2')
            print(f'Iteration {ii}: \nResidual error: {norm:.4e}')
            if norm < self.tol:
                self.delta_s = sqrt(self.u_f.inner(self.u_f) + self.psi * self.lmbda**2 * self.Q.inner(self.Q))
                self.counter = 1
                break

            ii+=1
            assert ii <= self.max_iter, 'Newton Solver not converging'

            du_f = Vector()
            du = Vector()
            solve(K_star, du_f, R_star, self.solver) # solve reduced problem
            self.C.mult(du_f, du) # convert to global displacement vector
            self.__update_nodal_values(self.u.vector()-du) # update solution
            self.u_f -= du_f # update free DoFs for calculation of s
        
    def solve(self):
        '''
        Main function to increment through the arc-length scheme. 
        '''
        if self.counter == 0:
            print('Initializing solver parameters...')
            self.__initial_step()
        print('\nArc-Length Step', self.counter,':')
        # initialization
        u_update = Vector()
        u_update.init(self.u.vector().size())
        
        R_star = Vector()
        self.C.init_vector(R_star,1)
        
        if self.counter == 1:
            self.converged = False

            self.u_f_n, self.u_f_n_1 = Vector(), Vector()
            self.u_f_n.init(self.u_f.size())
            self.u_f_n_1.init(self.u_f.size())
            self.lmbda_n = 0
            self.lmbda_n_1 = 0

        # Predictor Step: 
        else:
            alpha = self.delta_s / self.delta_s_n
            self.u_f = (1+alpha) * self.u_f_n - alpha * self.u_f_n_1 # update the free nodes
            self.C.mult(self.u_f, u_update) # free nodes to global displacement vector
            self.__update_nodal_values(u_update)            
            self.lmbda = (1+alpha) * self.lmbda_n - alpha * self.lmbda_n_1 # update displacement factor
            
        # Apply boundary conditions (both homogenous and nonhomogenous)
        self.displacement_factor.t = self.lmbda
        for bc in self.bcs:
            bc.apply(self.u.vector())
                    
        delta_u_f = self.u_f.copy() - self.u_f_n
        delta_lmbda = self.lmbda - self.lmbda_n

        self.converged_prev = self.converged
        self.converged = False
        
        # Corrector Step(i.e. arc-length solver):
        solver_iter = 0
        norm = 1
        while norm > self.tol and solver_iter < self.max_iter:
            solver_iter += 1
            
            # Assemble K and R
            K = assemble(self.J)
            R = assemble(self.residual)

            #Get K*, F*, and R* (reduced matrices and vectors)
            K_mat = as_backend_type(K).mat()
            temp = self.C_mat.copy().transpose().matMult(K_mat)
            K_star_mat = temp.matMult(self.C_mat) 
            K_star = Matrix(PETScMatrix(K_star_mat)) # reduced stiffness matrix
            PETScMatrix(temp).mult(-self.u_p, self.Q) # vector of Dirichlet BCs
            self.C.transpmult(R, R_star) # reduced residual vector
            
            QQ = self.Q.inner(self.Q)

            # solve for d_lmbda, d_u:
            a = 2 * delta_u_f
            b = 2 * self.psi * delta_lmbda * QQ

            A = delta_u_f.inner(delta_u_f) + self.psi * delta_lmbda**2 * QQ - self.delta_s**2
            R_star_norm = R_star.norm('l2')
            norm = sqrt(R_star_norm**2 + A**2)
            print(f'Iteration: {solver_iter} \n|Total Norm: {norm:.4e} |Residual Norm: {R_star_norm:.4e} |A: {A:.4e}|')
            if norm < self.tol:
                self.converged = True
                break

            du_f_1 = Vector()
            du_f_2 = Vector()

            solve(K_star, du_f_1, self.Q, self.solver)
            solve(K_star, du_f_2, R_star, self.solver)

            dlmbda = (a.inner(du_f_2) - A) / (b + a.inner(du_f_1))
            du_f = -du_f_2 + dlmbda * du_f_1

            # update delta_u, delta_lmbda, u, lmbda
            delta_lmbda += dlmbda
            self.lmbda += dlmbda
            delta_u_f += du_f
            self.u_f += du_f
            self.C.mult(self.u_f, u_update)
            self.__update_nodal_values(u_update)
            self.displacement_factor.t = self.lmbda
            for bc in self.bcs:
                bc.apply(self.u.vector())
            

         # Solution Update
        if self.converged:
            if self.counter == 1:
                self.delta_s_max = self.delta_s
                self.delta_s_min = self.delta_s / 1024.0
            
            self.delta_s_n = self.delta_s
                
            self.u_f_n_1 = self.u_f_n
            self.lmbda_n_1 = self.lmbda_n
        
            self.u_f_n = self.u_f.copy()
            self.lmbda_n = self.lmbda
                
            self.counter +=1 
            
            if self.converged_prev:
                self.delta_s = min(max(2*self.delta_s, self.delta_s_min), self.delta_s_max)
                
                
        else:
            if self.converged_prev:
                self.delta_s = max(self.delta_s / 2, self.delta_s_min)
            else:
                self.delta_s = max(self.delta_s / 4, self.delta_s_min)