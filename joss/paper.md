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

FEniCS [@alnaes2015fenics, logg2012automated] is an open-source software that solves partial differential equations (PDEs) through finite element analysis (FEA). Building on the versatile FEniCS platform, researchers have created both detailed tutorials [@bleyer2018numericaltours] and software packages [@finsberg2019pulse, rodriguez2019fenics] to facilitate research in specific application domains. For example, there are multiple open-source software libraries designed to increase the existing capabilities of FEniCS [@farrell2016computation, Kamensky2019, Mitusch2019]. Building on this established precedent, we have implemented our arclength solver as an add-on library on top of FEniCS, and also demonstrated usage of our solver via multiple documented examples involving nonlinear continuum and beam finite elements at large deformations. Notably, our intent is to keep the solver usage similar to already established FEniCS solvers to enable off-the-shelf implementation within FEniCS workflows. To this end, our arclength solver can be integrated into existing FEniCS code that is based on a Newton solver just by changing the choice of solver. Furthermore, our code and documentation are organized such that more advanced users can make custom modifications as needed, for example pinpointing and branch switching algorithms [@wriggersQuadraticallyConvergentProcedure1988]. 

In brief, our implementation follows a recently published method [@kadapa2021simple] adopting an extrapolated predictor scheme using the previously converged two solution steps, and solves the augmented arc-length matrix equation as the corrector scheme. The main advantage of the extrapolation prediction scheme is in its ability to determine the forward direction of the equilibrium path without doing extra computations such as evaluating the determinant of the stiffness matrix~[@crisfield1981fast], although calculating the sign of the determinant can be simple depending on the type of linear solver (i.e. the sign of the determinant can be extracted furing numerical factorization ). For displacement-controlled arclength problems, similar to a multifreedom constraint problem, we introduced a constraint matrix to augment the prescribed non-zero Dirichlet boundary conditions and the free degrees of freedoms to solve the arc-length equations [@batoz1979incremental,verhoosel2009dissipation]. More details about the theory, implementation, and usage of our arc-length solver and general mathematical preliminaries can be found in our extensive ReadTheDocs documentation (<https://fenics-arclength.readthedocs.io/en/latest/>). 

# Statement of Need
Materials and structures that undergo large deformation often exhibit fascinating mechanical instabilities due to their geometrically nonlinear behavior. Recently, there has been a renewed interest in studying these mechanical instabilities in engineered architected materials [@medina2020navigating], and in leveraging our understanding of mechanics to interpret the non-linear behavior of natural materials [@lee2021elastic] and the elastocapillary behavior of soft solids [@dortdivanliogluPlateauRayleighInstability2022a]. In these applied mechanics problems, the finite element method is often used to computationally simulate these complex structures, particularly when geometry and/or material behavior are complex and consequently previously derived and validated analytical solutions are not applicable. In the solutions of geometrically nonlinear problems with finite element analysis, the equilibrium path of the structure is traced either in a force-controlled scheme (i.e. prescribed incremental Neumann boundary conditions) or displacement-controlled scheme (i.e. prescribed incremental non-zero Dirichlet boundary conditions). However, in extreme cases, such as in the case of snap-back or snap-through instabilities, standard Newton solver based force/displacement-controlled schemes fail. To address this shortcoming of standard Newton solvers, various nonlinear solvers have been proposed based on the common approach of adding an extra constraint equation on top of the Newton solver [@leon2011unified]. In particular, the arclength method addresses this challenge through an arclength constraint equation where an incremental load or displacement factor at each load iteration, an additional unknown, is modified. Various arclength-type solvers are available in commercial FEA software (e.g., see Riks method in ABAQUS [@abaqus]) to capture the full equilibrium path. The goal of this work is to provide the research community with an open-source and robust arclength solver that can readily be applied to a myriad of nonlinear solid mechanics problems.

# Examples of FEniCS Arclength

![A demonstration of our arclength solver via the snap-back instability behavior of Lee's frame, a popular benchmark problem in nonlinear structural mechanics: a) The problem definition for Lee's frame; b) The equilibrium path from our solver compared to the literature [@kadapa2021simple]; c) The final deformed configuration from our solver compared to the literature [@kadapa2021simple].\label{fig:leesframe}](joss_arclength.png)

In Fig. \autoref{fig:leesframe}, we show an example of our arclength implementation for solving the Lee's frame problem [@lee1968large], a popular benchmark problem in nonlinear structural mechanics. In this problem, the frame is discretized using 2D Simo-Reissner beams [@simo1986three] with mixed elements. The elements have linear Lagrange interpolation functions for both displacement and rotation and use reduced integration to avoid locking. Figure \autoref{fig:leesframe}a) gives an overview of the Lee's frame problem and illustrates the equilibrium path (Figure \autoref{fig:leesframe}b) and the deformed configuration (Figure \autoref{fig:leesframe}c) in comparison with the literature [@kadapa2021simple]. The Lee's frame example demonstrates that our solver is able to trace unstable states, where traditional Newton solvers would fail. Additional validation examples on problems in nonlinear structural mechanics with beam elements (compared to solutions from Kadapa et al. [@kadapa2021simple]), and bilayer wrinkling with continuum elements (compared to analytical solutions [@allensandwich, budday2017wrinkling, cao2012wrinkling]) can be found in the examples folder of the Github Repository (<https://github.com/pprachas/fenics_arclength/tree/master/examples>), or in the ReadTheDocs (<https://fenics-arclength.readthedocs.io/en/latest/index.html#notebook-examples>) documentation. Critically, we demonstrated that our solver is able to trace the complex equilibrium paths of popular benchmark problems such as $215^\circ$ and $180^\circ$ archs (see Table {#tab:examples} for all the provided examples in our documentation). In addition to validation problems, our documentation also contains the solver's API, example usage of our solver for continuum problems, and example usage of our displacement controlled arclength scheme.

::: {#tab:examples}
+--------------------------------------+
| **Force-controlled Scheme**          |
+--------------------------------------+
| *Continuum elements*                 |
|                                      |
| -   Lee's Frame (Continuum Elements) |
|                                      |
| *Beam elements*                      |
|                                      |
| -   Lee's Frame $^\S$                |
|                                      |
| -   $215^\circ$ arc $^\S$            |
|                                      |
| -   $180^\circ$ arc $^\S$            |
|                                      |
| -   3D Helical Beam                  |
+--------------------------------------+
| **Displacement-controlled Scheme**   |
+:=====================================+
| *Continuum Elements*                 |
|                                      |
| -   Bilayer Wrinkling $^\bigstar$    |
|                                      |
| *Beam Elements*                      |
|                                      |
| -   Random Fiber Network             |
+--------------------------------------+

: A list of all available examples in our GitHub Repository and
ReadTheDocs documentation. Examples marked with ($\S$) are compared with
solutions from the literature obtained from [@kadapa2021simple].
Examples with ($\bigstar$) are compared with solutions from
[@allensandwich; @budday2017wrinkling; @cao2012wrinkling] \label{fig:leesframe}
:::

[]{#tab:examples label="tab:examples"}

# Acknowledgements

This work was made possible with funding through the Boston University David R. Dalton Career Development Professorship, the Hariri Institute Junior Faculty Fellowship, the Haythornthwaite Foundation Research Initiation Grant, and the National Science Foundation Grant CMMI-2127864. This support is gratefully acknowledged.


# References
