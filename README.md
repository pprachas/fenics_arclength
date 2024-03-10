# FEniCS Implementation of the Arc-length Method
[![Documentation Status](https://readthedocs.org/projects/fenics-arclength/badge/?version=latest)](https://fenics-arclength.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/12Ss1RN7Smyrnn2671JloGphtiO1r_67J#scrollTo=AHwJEykSdkGp)
[![DOI](https://zenodo.org/badge/590268096.svg)](https://zenodo.org/doi/10.5281/zenodo.10563094)
[![status](https://joss.theoj.org/papers/c4791755683fc95dd5ccd9c8e8708a08/status.svg)](https://joss.theoj.org/papers/c4791755683fc95dd5ccd9c8e8708a08)


This repository contains the arc-length Riks solver written with FEniCS. Additional information, documentation, and example problems can be found our ReadTheDocs documentation:
  
<div align="center">

 |[ **Link to ReadTheDocs Documentation** ](https://fenics-arclength.readthedocs.io/en/latest/index.html)|
 |--------------------------------------------------------------------------|
 
<div align = "left">

**[The accompanying paper to this Github Repository can be found here.](https://joss.theoj.org/papers/10.21105/joss.05727)**

More information on the arc-length method and the solution approach can be found in:  
1. [Nonlinear Analysis of Structures: The Arc Length Method](https://scholar.harvard.edu/files/vasios/files/ArcLength.pdf)  
2. [Incremental displacement algorithms for nonlinear problems](https://onlinelibrary.wiley.com/doi/10.1002/nme.1620140811)  
3. [A simple extrapolated predictor for overcoming the starting and tracking issues in the arc-length method for nonlinear structural mechanics](https://arxiv.org/abs/2005.10192)  
4. [A dissipation-based arc-length method for robust simulation of brittle and ductile failure](https://onlinelibrary.wiley.com/doi/10.1002/nme.2447)

## Table of Contents
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Contents in this Repository](#contents)
* [Validation](#validation)
* [Theory](#theory)
* [Community Contributions](#community)


## Dependencies <a name="dependencies"></a>
This package relies on FEniCS 2019.1.0. (Note that this is the legacy version NOT FEniCSx). Brief installation instructions are outline below. For more information see the [official FEniCS installation instructions.](https://fenicsproject.org/download/archive/)

### FEniCS on Windows
The simplest way to install FEniCS on Windows 10 is to install [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) with Ubuntu distribution. Then you can follow the FEniCS installation instructions for a Linux machine.

### FEniCS on Ubuntu
 To install FEniCS on Ubuntu, run these commands:
 
        sudo apt-get install software-properties-common
        sudo add-apt-repository ppa:fenics-packages/fenics
        sudo apt-get update
        sudo apt-get install fenics

**Note that due to recent ufl updates and its issues with legacy dolfin there may be issues running fenics-arclength if you install FEniCS through PPA. See [here](https://fenicsproject.discourse.group/t/announcement-ufl-legacy-and-legacy-dolfin/11583) about known issues and solutions to ufl and legacy dolfin.**
        
### FEniCS on Anaconda (Linux and Mac only) *Recommended*:
 
        conda create -n fenicsproject -c conda-forge fenics
        conda activate fenicsproject

However, the conda FEniCS installation is missing a few important libraries (i.e. scipy, mshr). For convenience, we provide can environment.yml file that contains all the dependencies except for Jupyter notebooks. 
To use the environment.yml file, navigate the to root directory and use the following commands:
 
        conda create -n fenicsproject
        conda activate fenicsproject
        conda env update -f environment.yml

 For M1 macs there might be issues with installing FEniCS. As a workaround, you must first set the conda environment variable to osx-64. As such, the full command to install FEniCS and all the dependencies on an M1 Mac are:
 
        conda create -n fenicsproject
        conda activate fenicsproject
        conda config --env --set subdir osx-64
        conda env update -f environment.yml
 
While the validation scripts and the package can be used without Jupyter notebooks, Jupyter notebooks are required to run the examples in the examples directory. To install jupyter notebooks, the following command should be ruin after using the environment.ymk file:
 
        conda install -c conda-forge jupyter
 
### FEniCS on Docker (Windows, Mac, Linux)
First install [Docker Desktop](https://fenicsproject.org/download/archive/) then run the following command:

        curl -s https://get.fenicsproject.org | bash
You also can start the Docker container with the following command:

        docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:current
A more comprehensive and detailed instructions on Docker installation can be found here: [Docker Installation Instructions](https://fenics.readthedocs.io/projects/containers/en/latest/introduction.html).

**Note: For Docker installation, numpy, scipy, and matplotlib should be installed alongside FEniCS by default. To enable jupyter notebooks for FEniCS, please see this link: [Jupyter notebooks for Docker installation of FEniCS](https://fenics.readthedocs.io/projects/containers/en/latest/jupyter.html). However, the scripts in the ``validation`` can be run without Jupyter notebooks.**
 
## Installation
Once FEniCS has been installed, our package can easily be install through pip with the following command:

    pip install git+https://github.com/pprachas/fenics_arclength.git@master

Note that in the case of Docker installation, you might have to add the option --user to use pip installation (i.e. ``pip install --user git+https://github.com/pprachas/fenics_arclength.git@master``).
In cases where the whole github repository is needed, the github repo can first be cloned before installation:

    git clone https://github.com/pprachas/fenics_arclength.git
    cd fenics_arclength
    pip install .

In the case of developer's version, the last line can be replaced with ``pip install -e .``

## Running FEniCS-arclength on Google Colab
Our notebooks in the [examples directory](https://github.com/pprachas/fenics_arclength/tree/master/examples) can be run on Google Colab! 
Simply navigate to the notebook in the GitHub repository and change the domain from ``github.com`` to ``githubtocolab.com`` and run the notebook!

## Contents in this repository <a name="contents"></a>

| Folder| Contents of the folder|
|-------|--------|
|[arc_length](arc_length)| contains the code for our arc-length implementation; both force and displacement control solvers are contained there |
|[docs](docs)| the build files for our [readthedocs documentation](https://fenics-arclength.readthedocs.io/en/latest/) |
|[examples](examples) | contains Jupyter notebook examples to use our arc-length implementation. Note that Jupyter notebooks has to be installed in the FEniCS environment for the notebooks to run. | 
|[validation](validation)| contains python scripts to compare our solver with analytical solutions/solutions in literature. |

## Validation <a name="validation"></a>
To validate that our arc-length solver works we provide 3 validation examples. To run the examples go to the ``validation`` directory and run ``python3 validate_xx.py``. The output plots will be in the ``plots`` folder. The visualize the Paraview files, append ``-p`` at the end of the command line (i.e. ``python3 validate_xx.py -p``). The available python scripts are:
* ``validate_3Dbeamsmall.py``
    * This script solves a clamped cantilever beam with a small applied force and moment at the free end. The solution (i.e. the reaction shear, moment, and curvature) of from the arc-length solver is compared with linear beam theory.
    * *Outputs:* The outputs of the script is the percent differences between the analytical solution and arc-length solution for reaction shear, reaction moment, and beam curvature. If the solutions are within 1% difference, then the validation is complete.
* ``validate_3Dbeamlarge.py``
    * This script solves a clamped cantilever beam with a large applied force and moment at the free end. The curvature of the beam in the same direction of the applied moment is compared with the moment-curvature relation ($\\kappa = \frac{M}{EI}$) since the beam constitutive model is linear elastic and does not couple deformation modes. 
    * *Outputs:* The outputs of the script is the percent differences between the analytical solution and arc-length solution for beam curvature. If the solutions are within 1% difference, then the validation is complete.
* ``validate_leeframe.py``
    * This scripts solves Lee's frame, a popular benchmarking problem in nonlinear solid mechanics. The resulting equilibrium path and critical buckling load is compared with literature obtained [here](https://www.sciencedirect.com/science/article/pii/S014102962034356X).
    * *Outputs:* The outputs of the script are the Pearson correlation coefficient of the equilibrium paths between our arc-length solver and the solution from literature. Both of the equilibrium paths are also plotted and saved in ``valiation/validation_leeframe.png``.
 * ``validate_bilayer.py``
    * This scripts solve the bilayer wrinkling problem. The resulting wrinkling wavelength and critical buckling strain is compared with literature obtained [here](https://royalsocietypublishing.org/doi/epdf/10.1098/rsta.2016.0163) and [here](https://groups.seas.harvard.edu/hutchinson/papers/WrinklingPhenonmena-JAMF.pdf).
    * *Outputs:* The outputs of the script are the percent differences between the analytical solutions (critical strain and wavelength) and FEA solution. The comparison plots are also saved in ``plots/validation_bilayer_stresstrain.png`` and ``plots/validation_bilayer_wavelength.png``.
    
 **Note that the beam validation scripts should be fast to run ($\sim$ 5 secs for small deformation and Lee's frame, $\sim$ 1 min for large deformation). The bilayer wrinkling will take longer to run ($\sim$ 25 mins).**

## Theory <a name="theory"></a>
Here is outline the basic theory of solving nonlinear finite elements and our implementation of the arc-length solver.
### Nonlinear Finite Elements
A nonlinear finite element problem seeks to minimize the residual vector that comes from discretizing the weak form of the energy balance equation (e.g. continuum for beam balance equations). In general the residual cannot be solved exactly and must be approximated through linearization. A common method to solve nonlinear finite element problems uses the Newton-Raphson method:

 ```math
\mathcal{R}(\mathbf{u}_{n+1}) = \mathcal{R}(\mathbf{u}_{n})+\frac{\partial \mathcal{R}(\mathbf{u}_{n})}{\partial \mathbf{u}_{n}}\Delta \mathbf{u}
 ```
 
where $\Delta u = \mathbf{u}_{n+1}-\mathbf{u}_n$.

Newton's method is solved incrementally until the desired convergence criterion. The term $\frac{\partial \mathcal R(\mathbf u_n)}{\partial \mathbf u_n}$
is typically called the tangential stiffness matrix $K_T$. The first term $\mathcal R(\mathbf u_n)$ is the difference between the internal force of the previous step $F^{ext}$, while the second term
$\frac{\partial \mathcal R(\mathbf u_n)}{\partial \mathbf u_n}\Delta \mathbf u$ is the correction of the internal force $F^{int}$, In general, the nonlinear problem is too difficult for the Newton solver to converge. As such, the external load is applied incrementally with the load factor $\lambda^k$ where $k$ is the increment. Putting it all together, the nonlinear problem can be written as:

```math
\mathcal{R}(\mathbf{u}_{n+1},\lambda_{n+1}) = F^{int}(\mathbf{u}_{n+1};\mathbf{u}_{n},\lambda_{n+1})-\lambda_{n+1} F^{ext}(\mathbf{u}_{n})
 ```

#### Conservative Loading
In most cases the external force does not depend on the solution $u$ (i.e. $F^{ext} (u_n) = F^{ext}$ ). These cases are called conservative loading. The problem than can be simplified to:

```math
\mathcal{R}(\mathbf{u}_{n+1},\lambda_{n+1}) = F^{int}(\mathbf{u}_{n+1};\mathbf{u}_{n})-\lambda^k F^{ext}
```

In this case the tangential stiffness matrix $K_T$ can be constructed using just the internal energy (i.e. $K_T = \frac{\partial F^{int}(\mathbf u_{n+1};\mathbf u_n)}{\partial \mathbf{u}_n}$)

#### Non-conservative loading
In the case where the external force depends on the solution $u$, the above assumption cannot be made and the whole residual must be taken into account when constructing the tangential stiffness matrix $K_T$. As a result, $K_T$ will be non-symmetric. Examples of these specific special cases are applied moments around a fixed axis, follower loads (i.e. loads that change direction based on the deformed configuration), pressure loads, etc.

### The Arc-length method
One of the main drawbacks of Newton's method is its inability to trace equilibrium paths with limit points. As a workaround, the load parameter $\lambda_n$ is now also an unknown parameter at each increment, and additional arc-length constraint is added. In this repository, we implement both the arc-length method for force control (i.e. problems with force boundary conditions) and displacement control (i.e problems with non-homogenous displacement boundary conditions).
#### Force Control
The additional arc-length constraint for force control is:

```math
\mathcal{A}(\mathbf{\mathbf{u}_{n+1}},\lambda_{n+1}) = \Delta\mathbf{u}^T\Delta\mathbf{u} + \psi\Delta\lambda^2 F_{ext}(\mathbf{u}_{n})^T F_{ext}(\mathbf{u}_{n})-(\Delta s)^2
```

where $\Delta s$ determines how far to search for the next equilibrium point and $\psi$ is the arc length parameter that gives you different arc-length solver schemes. When $\psi = 1$ (as like the examples in this repository), the arc-length equation is also known as the *spherical arc-length method*, and when $\psi = 0$ the *cylindrical arc-length* method is recovered.

#### Displacement Control
Sometimes instead of prescribing traction, the problem has a boundary condition with prescribed non-zero displacement (i.e. non-homogenous Dirichlet boundary conditions). In this case, similar to Ref.2, the problem is formulated similar to a multifreedom constraint and we construct a constraint matrix $C$ such that: 

```math
\mathbf{u} = C\mathbf{u}_f+\lambda \mathbf{u}_p
```

where $u_f$ and $u_p$ are the free and prescribed displacement nodes respectively, and $\lambda$ is the incremental displacement factor.


The arc length equation needs to be modified and now becomes:

```math
 \mathcal{A}(\mathbf{u}_f,\lambda) = \Delta\mathbf{u}_f^T\Delta\mathbf{u}_f + \psi\Delta\lambda^2Q^TQ-(\Delta s)^2
 ```

where:

```math
 Q = C^TK_T\mathbf{u}_p
```

### Predictor-Corrector Scheme
The predictor our arc-length implementation for both the force and displacement control scheme follows the implementation from Ref.3. The prediction step takes in the previous solution and extrapolates where:

```math
\mathbf{u}_{n+1}^{predicted} = [1+\alpha] \mathbf{u}_{n} -\alpha u_{n-1} \\
\lambda_{n+1}^{predicted} = [1+\alpha] \lambda_n -\alpha \lambda_{n-1}
```

where $\alpha=\frac{\Delta s_n}{\Delta s_{n-1}}$ is the extrapolation parameter that depends on the arc-length parameter for the previous and current step. Using this extrapolation scheme both provides a good initial guess for the next equilibrium solution as well as identifies the correct direction for the next point in the equilibrium path. Our implementation allows for easy modification of this extrapolation scheme, as seen in the 3D beam example.

The following Ref.3 and Ref.4, force control corrector scheme solves the augmented matrix equation:

```math
\begin{bmatrix}
K_T & -F^{ext} \\
\frac{\partial \mathcal{A}}{\partial u} & \frac{\partial \mathcal{A}}{\partial \lambda}
\end{bmatrix} 
\begin{bmatrix}
\delta \mathbf{u} \\ \delta \lambda
\end{bmatrix}
= 
\begin{bmatrix}
\mathcal{R} \\ \mathcal{A}
\end{bmatrix}
```

where 

```math
\Delta \mathbf{u}_{n+1} = \Delta \mathbf{u}_n + \delta \mathbf{u}  
```

```math
\Delta \lambda_{n+1} = \Delta \lambda_n + \delta \lambda
```

The displacement control corrector scheme modifies the above equation to:

```math
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
```

Similar to Ref. 3 and Ref. 4, we solve the block system of equations by part using the Sherman–Morrison formula. For more details refer to the Ref 3 and Ref 4. 

## Community Contributions

If you find any issues, bugs or problems with FEniCS-arclength, please use report them using the [GitHub issue tracker](https://github.com/pprachas/fenics_arclength/issues). For efficiency, please include a minimum working example (to reproduce the error/bug) as well as the error itself.

We also welcome any questions and suggestions in the using [Github Discussions](https://github.com/pprachas/fenics_arclength/discussions)!

## Citation
If you find this repository helpful, we would appreciate it if you cite us!
```
@article{Prachaseree_FEniCS-arclength_A_numerical_2024,
author = {Prachaseree, Peerasait and Mohammadzadeh, Saeed and Dortdivanlioglu, Berkin and Lejeune, Emma},doi = {10.21105/joss.05727},
journal = {Journal of Open Source Software},
month = feb,
number = {94},
pages = {5727},
title = {{FEniCS-arclength: A numerical continuation package in FEniCS for nonlinear problems in solid mechanics}},
url = {https://joss.theoj.org/papers/10.21105/joss.05727},
volume = {9},
year = {2024}}
```
