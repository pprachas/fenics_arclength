.. default-role:: math

Math Preliminaries
==================

Here is the outline for the basic theory of solving nonlinear finite elements and our implementation of the arc-length solver.

Nonlinear Finite Elements
-------------------------
A nonlinear finite element problem seeks to minimize the residual vector that comes from discretizing the weak form of the energy balance equation (e.g. continuum or beam balance equations). In general the residual cannot be solved exactly and must be approximated through linearization. A common method to solve nonlinear finite element problems uses the Newton-Raphson method:

.. math:: \mathcal{R}(\mathbf{u}_{n+1}) = \mathcal{R}(\mathbf{u}_{n})+\frac{\partial \mathcal{R}(\mathbf{u}_{n})}{\partial \mathbf{u}_{n}}\Delta u
 
where `\Delta \mathbf{u} = \mathbf{u}_{n+1}-\mathbf{u}_n`.

Newton's method is solved incrementally until the desired convergence criterion. The term `\frac{\partial \mathcal R(\mathbf u_n)}{\partial \mathbf u_n}`
is typically called the tangential stiffness matrix `K_T`. The first term `\mathcal R(\mathbf{u_n})` is the difference between the internal force of the previous step and applied force `F^{ext}`, while the second term
`\frac{\partial \mathcal R(\mathbf u_n)}{\partial \mathbf u_n}\Delta \mathbf{u}` is the correction of the internal force `F^{int}`. In general, the nonlinear problem is too difficult for the Newton solver to converge. As such, the external load is applied incrementally with the load factor `\lambda^k` where `k` is the increment. Putting it all together, the nonlinear problem can be written as:

.. math:: \mathcal{R}(\mathbf{u}_{n+1},\lambda_{n+1}) = F^{int}(\mathbf{u}_{n+1};\mathbf{u}_{n},\lambda_{n+1})-\lambda_{n+1} F^{ext}(\mathbf{u}_{n})

Conservative Loading
#####################

In most cases the external force does not depend on the solution `\mathbf{u}` (i.e. `F^{ext} (\mathbf{u}_n) = F^{ext}` ). These cases are called conservative loading. The problem than can be simplified to:

.. math:: \mathcal{R}(\mathbf{u}_{n+1},\lambda_{n+1}) = F^{int}(\mathbf{u}_{n+1};\mathbf{u}_{n})-\lambda^k F^{ext}

In this case the tangential stiffness matrix `K_T` can be constructed using just the internal energy (i.e. `K_T = \frac{\partial F^{int}(\mathbf u_{n+1};\mathbf u_n)}{\partial \mathbf{u}_n}`)

Non-conservative loading
########################
In the case where the external force depends on the solution `\mathbf{u}`, the above assumption cannot be made and the whole residual must be taken into account when constructing the tangential stiffness matrix `K_T`. As a result, `K_T` will be non-symmetric. Examples of these specific special cases are applied moments around a fixed axis, follower loads (i.e. loads that change direction based on the deformed configuration), pressure loads, etc.

The Arc-length method
---------------------

One of the main drawbacks of Newton's method is its inability to trace equilibrium paths with limit points. As a workaround, the load parameter `\lambda_n` is now also an unknown parameter at each increment, and additional arc-length constraint is added. In this repository, we implement both the arc-length method for force control (i.e. problems with force boundary conditions) and displacement control (i.e problems with non-homogenous displacement boundary conditions).

Force Control
#############

The additional arc-length constraint for force control is:

.. math:: \mathcal{A}(\mathbf{\mathbf{u}_{n+1}},\lambda_{n+1}) = \Delta\mathbf{u}^T\Delta\mathbf{u} + \psi\Delta\lambda^2 F_{ext}(\mathbf{u}_{n})^T F_{ext}(\mathbf{u}_{n})-(\Delta s)^2

where `\Delta s` determines how far to search for the next equilibrium point and `\psi` is the arc length parameter that gives you different arc-length solver schemes. When `\psi = 1` (as like the examples in this repository), the arc-length equation is also known as the *spherical arc-length method*, and when `\psi = 0` the *cylindrical arc-length* method is recovered.

Displacement Control
#####################

Sometimes instead of prescribing traction, the problem has a boundary condition with prescribed non-zero displacement (i.e. non-homogenous Dirichlet boundary conditions). In this case, similar to Ref.2, the problem is formulated similar to a multifreedom constraint and we construct a constraint matrix `C` such that: 

.. math:: \mathbf{u} = C\mathbf{u}_f+\lambda \mathbf{u}_p

where `\mathbf{u}_f` and `\mathbf{u}_p` are the free and prescribed displacement nodes respectively, and `\lambda` is the incremental displacement factor.


The arc length equation needs to be modified and now becomes:

.. math::\mathcal{A}(\mathbf{u}_f,\lambda) = \Delta\mathbf{u}_f^T\Delta\mathbf{u}_f + \psi\Delta\lambda^2Q^TQ-(\Delta s)^2

where:

.. math:: Q = C^TK_T\mathbf{u}_p

Predictor-Corrector Scheme
--------------------------

The predictor our arc-length implementation for both the force and displacement control scheme follows the implementation from Ref.3. The prediction step takes in the previous solution and extrapolates where:

.. math:: \mathbf{u}_{n+1}^{predicted} = [1+\alpha] \mathbf{u}_{n} -\alpha u_{n-1} \\ \lambda_{n+1}^{predicted} = [1+\alpha] \lambda_n -\alpha \lambda_{n-1}

where `\alpha=\frac{\Delta s_n}{\Delta s_{n-1}}` is the extrapolation parameter that depends on the arc-length parameter for the previous and current step. Using this extrapolation scheme both provides a good initial guess for the next equilibrium solution as well as identifies the correct direction for the next point in the equilibrium path. Our implementation allows for easy modification of this extrapolation scheme, as seen in the 3D beam example.

Following Ref.3 and Ref.4 (see below), force control corrector scheme solves the augmented matrix equation:

.. math:: 
    \begin{bmatrix} 
    K_T & -F^{ext} \\ 
    \frac{\partial \mathcal{A}}{\partial \mathbf{u}} & \frac{\partial \mathcal{A}}{\partial \lambda}
     \end{bmatrix} \begin{bmatrix} \delta \mathbf{u} \\ \delta \lambda \end{bmatrix} = \begin{bmatrix} \mathcal{R} \\ \mathcal{A} \end{bmatrix}

where 

.. math:: \Delta \mathbf{u}_{n+1} = \Delta \mathbf{u}_n + \delta \mathbf{u}  

.. math:: \Delta \lambda_{n+1} = \Delta \lambda_n + \delta \lambda

The displacement control corrector scheme modifies the above equation to:

.. math::
    \begin{bmatrix}
    C^\top K_T C & C^\top K \mathbf{u}_p \\
    \frac{\partial \mathcal{A}}{\partial \mathbf{u}_f} & \frac{\partial \mathcal{A}}{\partial \lambda}
    \end{bmatrix} 
    \begin{bmatrix}
    \delta \mathbf{u} \\ \delta \lambda
    \end{bmatrix}
    = 
    \begin{bmatrix}
    \mathcal{R} \\ \mathcal{A}
    \end{bmatrix}


Similar to Ref. 3 and Ref. 4, we solve the block system of equations by part using the Sherman–Morrison formula. For more details refer to the Ref 3 and Ref 4. 

Additional Resources
--------------------
More information on the arc-length method and the solution approach can be found in:

#. `Nonlinear Analysis of Structures: The Arc Length Method <https://scholar.harvard.edu/files/vasios/files/ArcLength.pdf>`_  
#. `Incremental displacement algorithms for nonlinear problems <https://onlinelibrary.wiley.com/doi/10.1002/nme.1620140811>`_  
#. `A simple extrapolated predictor for overcoming the starting and tracking issues in the arc-length method for nonlinear structural mechanics <https://arxiv.org/abs/2005.10192>`_  
#. `A dissipation-based arc-length method for robust simulation of brittle and ductile failure <https://onlinelibrary.wiley.com/doi/10.1002/nme.2447>`_
