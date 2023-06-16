---
title: "FEniCS-arclength: A numerical continuation package in FEniCS for nonlinear problems in solid mechanics"
tags:
  - Python
  - FEniCS
  - Finite Element Analysis
  - Solid Mechanics
authors:
  - name: Peerasait Prachaseree
    orcid: 0009-0000-3325-1410
    affiliation: 1 
  - name: Saeed Mohammadzadeh
    orcid: 0000-0001-9879-044X
    affiliation: 2
  - name: Berkin Dortdivanlioglu
    orcid: 0000-0001-7105-1452
    affiliation: "3, 4"
  - name: Emma Lejeune
    orcid: 0000-0001-8099-3468
    corresponding: true
    affiliation: 1
affiliations:
 - name: Department of Mechanical Engineering, Boston University, Massachusetts, the United States of America
   index: 1
 - name: Department of Systems Engineering, Boston University, Massachusetts, the United States of America
   index: 2
 - name: Department of Civil, Architectural and Environmental Engineering, The University of Texas at Austin, Austin, the United States of America
   index: 3
 - name: Oden Institute for Computational Engineering and Sciences, The University of Texas at Austin, Austin, the United States of America
   index: 4
date: 16 June 2023
bibliography: paper.bib
---

# Summary

FEniCS-arclength is a package that implements the numerical arcLength continuation method, following the pioneering work in developing the method [@riks1979incremental,crisfield1981fast], as a nonlinear solver in the open source finite element software framework FEniCS [@alnaes2015fenics, logg2012automated]. The arclength method is a robust nonlinear solver that is capable of capturing complex equilibrium paths including unstable states in solid mechanics problems with structural instabilities.
Such structural instabilities are induced by geometric nonlinearities, whereas material and contact nonlinearities can also be accommodated in the framework. 

FEniCS [@alnaes2015fenics, logg2012automated] is an open-source software that solves partial differential equations (PDEs) through finite element analysis (FEA). Building on the versatile FEniCS platform, researchers have created both detailed tutorials [@bleyer2018numericaltours] and software packages [@finsberg2019pulse, rodriguez2019fenics] to facilitate research in specific application domains. For example, there are multiple open-source software libraries designed to increase the existing capabilities of FEniCS [@farrell2016computation, Kamensky2019, Mitusch2019]. Building on this established precedent, we have implemented our arclength solver as an add-on library on top of FEniCS, and also demonstrated usage of our solver via multiple documented examples involving nonlinear continuum and beam finite elements at large deformations. Notably, our intent is to keep the solver usage similar to already established FEniCS solvers to enable off-the-shelf implementation within FEniCS workflows. To this end, our arclength solver can be integrated into existing FEniCS code that is based on a Newton solver just by changing the choice of solver. Furthermore, our code and documentation are organized such that more advanced users can make custom modifications as needed, for example pinpointing and branch switching algorithms~[@wriggersQuadraticallyConvergentProcedure1988]. 

In brief, our implementation follows a recently published method [@kadapa2021simple] adopting an extrapolated predictor scheme using the previously converged two solution steps, and solves the augmented arc-length matrix equation as the corrector scheme. The main advantage of the extrapolation prediction scheme is in its ability to determine the forward direction of the equilibrium path without doing extra computations such as evaluating the determinant of the stiffness matrix~[@crisfield1981fast], although calculating the sign of the determinant can be simple depending on the type of linear solver (i.e. the sign of the determinant can be extracted furing numerical factorization ). For displacement-controlled arclength problems, similar to a multifreedom constraint problem, we introduced a constraint matrix to augment the prescribed non-zero Dirichlet boundary conditions and the free degrees of freedoms to solve the arc-length equations~[@batoz1979incremental,verhoosel2009dissipation]. More details about the theory, implementation, and usage of our arc-length solver and general mathematical preliminaries can be found in our extensive ReadTheDocs documentation (<https://fenics-arclength.readthedocs.io/en/latest/>). 

# References