.. fenics-arclength documentation master file, created by
   sphinx-quickstart on Mon Jan 23 12:24:26 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Fenics-arclength Documentation
============================================

**fenics-arclength** is a Python implementation of the arclength solver built on top of FEniCS. 
The Arc Length Method, sometimes referred to as the Riks Method, is a method used to solve solid mechanics problems with geometric nonlinearity
with complex equilibrium paths. 
The library aims to keep the usage as similar to `FEniCS <https://fenicsproject.org/download/archive/>`_ (version 2019.1.0) to allow
for off-the-shelf implementation and integration with other FEniCS workflows.


.. button-link:: https://github.com/pprachas/fenics_arclength
   :color: secondary
   :expand:

   Link to Github Repository


.. image:: imgs/fenics_project.png
   :align: center
   :width: 700

Usage
-----

Dependencies
############

This package solely depends on FEniCS 2019.1.0. It is not tested on other packages.
`Instructions to install FEniCS (legacy) can be found here. <https://fenicsproject.org/download/archive/>`_
Our `GitHub repository <https://github.com/pprachas/fenics_arclength>`_ also contains more detailed installation instructions. 
For conda installations (available only on Mac and Linux), we also provide an environment.yml file in our GitHub repository.

Package Installation
####################
To use the arc-length solver, you can install it using pip or by downloading and adding the repository to your Python path.

Install directly using pip
**************************

The recommended way to install the arc-length solver is to use pip. To do this, run the following command in your terminal:

.. code-block:: bash

   pip install git+https://github.com/pprachas/fenics_arclength.git@master

This will install the latest version of the arc-length solver from the master branch of the repository.

Install by cloning the GitHub Repository
****************************************

If you prefer to download and install the arc-length solver and run the notebook examples/validation scripts, you can do so by following these steps:

1. Clone the repository to your local machine:
   
   .. code-block:: bash

      git clone https://github.com/pprachas/fenics_arclength.git

2. Change directory to the fenics_arclength directory:
   
   .. code-block:: bash

      cd fenics_arclength

3. Install the arc-length solver:

   .. code-block:: bash

      pip install .

This will install the arc-length solver from the source code in the fenics_arclength directory.

In the case of developer's version, the last step can be replaced as ``pip install -e .``

Once the arc-length solver is installed, you can import it into your Python scripts and use it to solve arc-length problems.


Theory
--------------------------
The general concepts behind non-linear finite element analysis (FEA) are introduced here for both educational purposes and completeness. In addition, we briefly describe the implementation of both arclength solvers in our package. We also provide an example on modifying the arc-length predictor scheme for geometrically exact beams. The links below contain significant additional detail.

.. toctree::
   :maxdepth: 1

   math_prelim
   examples/force_control/beam/README

Arc-Length solvers in this package
----------------------------------

.. toctree::
   :maxdepth: 2

   modules


Notebook Examples
----------------------------------
To execute the Jupyter Notebooks in this repository, please download the
`repository <https://github.com/pprachas/fenics_arclength>`_ and use the notebooks in the ``examples`` folder. Note that downloading
the notebook directly from ReadTheDocs will not guarantee that the examples will run since the files will be
missing other components (i.e. meshes, other scripts).  

.. toctree::
   :maxdepth: 3

   nb_displacement
   
.. toctree::
   :maxdepth: 3 

   nb_force


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
