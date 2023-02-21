from dolfin import *
import numpy as np
from math import sqrt

'''
An FEniCS based arc-length with prescribed traction. The code is heavily based on the paper:
Kadapa, Chennakesava. "A simple extrapolated predictor for overcoming the starting and tracking issues in the arc-length method for nonlinear structural mechanics." Engineering Structures 234 (2021): 111755.

'''
class force_control:
    ''' The arc-length displacement control solver of this library
    
    Args:
        psi: the scalar arc-length parameter. When psi = 1, the method becomes the shperical arc-length method and when psi = 0 the method becomes the cylindrical arc-length method
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

    def __init__(self, psi, tol, lmbda0, max_iter, u, F_int, F_ext, bcs, J, load_factor, solver='default'):
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
        self.load_factor = load_factor
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
        '''        
        
        ii=0
        print('Starting initial Force Control with Newton Method:')
        
        #------------Calculate inner(F,F)---------------#
        self.load_factor.t = 1 # for construction of FF
        self.F_ext_vec = assemble(self.F_ext)
        
        for bc in self.bcs:
            bc.apply(self.F_ext_vec)
        
        self.FF = self.F_ext_vec.inner(self.F_ext_vec)
        
        self.load_factor.t = self.lmbda
        
        while True:
            # Assemble and apply BCs
            K,R = assemble_system(self.J, self.residual, self.bcs)
            norm = R.norm('l2')
            print(f'Iteration {ii}: \nResidual error: {norm:.4e}')
            if norm < self.tol:
                self.delta_s = sqrt(self.u.vector().inner(self.u.vector()) + self.psi * self.lmbda**2 * self.FF)
                self.counter = 1
                break
            
            ii += 1
            assert ii <= self.max_iter, 'Newton Solver not converging'

            du = Vector()
            solve(K, du, R, self.solver)
            self.__update_nodal_values(self.u.vector()-du)
        
    def solve(self):
        '''
        Main function to increment through the arc-length scheme. 
        '''
        if self.counter == 0:
            print('Initializing solver parameters...')
            self.__initial_step()

        
        print('\nArc-Length Step', self.counter,':')
        # initialization
        if self.counter == 1:
            self.converged = False

            self.u_n, self.u_n_1 = Vector(), Vector()
            self.u_n.init(self.u.vector().size())
            self.u_n_1.init(self.u.vector().size())
            self.lmbda_n = 0
            self.lmbda_n_1 = 0

        # Predictor Step:
        else:
            alpha = self.delta_s / self.delta_s_n
            self.__update_nodal_values((1+alpha) * self.u_n - alpha * self.u_n_1)
            self.lmbda = (1+alpha) * self.lmbda_n - alpha * self.lmbda_n_1
         
        self.load_factor.t = self.lmbda
        delta_u = self.u.vector().copy() - self.u_n
        delta_lmbda = self.lmbda - self.lmbda_n
           
        self.converged_prev = self.converged
        self.converged = False
        
        # Corrector Step(i.e. arc-length solver):
        solver_iter = 0
        norm = 1
        while norm > self.tol and solver_iter < self.max_iter:
            solver_iter += 1
            
            # Assemble K and R
            K,R = assemble_system(self.J, self.residual, self.bcs) # assemble system

            # solve for d_lmbda, d_u:
            a = 2 * delta_u
            b = 2 * self.psi * delta_lmbda * self.FF

            A = delta_u.inner(delta_u) + self.psi * delta_lmbda**2 * self.FF - self.delta_s**2
            R_norm = R.norm('l2')
            norm = sqrt(R_norm**2 + A**2)
            print(f'Iteration: {solver_iter} \n|Total Norm: {norm:.4e} |Residual Norm: {R_norm:.4e} |A: {A:.4e}|')
            if norm < self.tol:
                self.converged = True
                break

            du_1 = Vector()
            du_2 = Vector()

            solve(K, du_1, self.F_ext_vec, self.solver)
            solve(K, du_2, R, self.solver)

            dlmbda = (a.inner(du_2) - A) / (b + a.inner(du_1))
            du = -du_2 + dlmbda * du_1

            # update delta_u, delta_lmbda, u, lmbda
            delta_u += du
            delta_lmbda += dlmbda
            self.__update_nodal_values(self.u.vector() + du)

            self.lmbda += dlmbda
            self.load_factor.t = self.lmbda

         # Solution Update
        if self.converged:
            if self.counter == 1:
                self.delta_s_max = self.delta_s
                self.delta_s_min = self.delta_s / 1024.0
            
            self.delta_s_n = self.delta_s
                
            self.u_n_1 = self.u_n
            self.lmbda_n_1 = self.lmbda_n
        
            self.u_n = self.u.vector().copy()
            self.lmbda_n = self.lmbda
                
            self.counter +=1 
            
            if self.converged_prev:
                # Predictor update rule if solution converges
                self.delta_s = min(max(2*self.delta_s, self.delta_s_min), self.delta_s_max)
                
                
        else:
            # Predictor update rule id solution doesn't converge
            if self.converged_prev:
                # Rule if previous step converged
                self.delta_s = max(self.delta_s / 2, self.delta_s_min)
            else:
                # Rule if previous step did not converge
                self.delta_s = max(self.delta_s / 4, self.delta_s_min)