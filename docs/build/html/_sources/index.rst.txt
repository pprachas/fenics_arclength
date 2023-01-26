.. fenics-arclength documentation master file, created by
   sphinx-quickstart on Mon Jan 23 12:24:26 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Fenics-arclength Documentation
============================================

**fenics-arclength** is a Python implementation of the arclength solver built on top of FEniCS. 
The arclength method, or the Riks method, is a method to generally used to solve solid mechanics problems
with complex equilibrium paths. 
The library aims to keep the usage as similiar to `FEniCS <https://fenicsproject.org/download/archive/>`_ (version 2019.1.0) to allow
for off-the-shelf usage.

Usage
-----

Dependencies
############

This package solely depends on FEniCS 2019.1.0. It is not tested on other packages.
`Instructions to install FEniCS (legacy) can be found here <https://fenicsproject.org/download/archive/>`_

Package Installation
####################
To use our arc-length solver, download and append this `repository <https://github.com/pprachas/fenics_arclength>`_ to the python path. This can be done by:

* Add directory to ``PYTHONPATH``:

.. code-block:: bash

   export PYTHONPATH <path/to/fenics_arclength>:$PYTHONPATH

* Append directory to path in python script:

.. code-block:: python

   import sys
   sys.path.append('path/to/fenics_arclength')


Theory
--------------------------
The general idea of non-linear finite element analysis (FEA) is introduced here. We also 
briefly describe the implementation of arclength solvers in our package.

.. toctree::
   :maxdepth: 1

   math_prelim

Arc-Length solvers in this package
----------------------------------

.. toctree::
   :maxdepth: 2

   modules


Notebook Examples
----------------------------------
To execute the Jupyter Notebooks in this repository, please download the
`<repository>`_ and use the notebooks in the ``examples`` folder. Note that downloading
the notebook directly will not guarantee that the examples will run since the files will be
missing other components (i.e. meshes, other scripts).  

.. toctree::
   :maxdepth: 2

   nb_displacement
   
.. toctree::
   :maxdepth: 3 

   nb_force


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
