.. fenics-arclength documentation master file, created by
   sphinx-quickstart on Mon Jan 23 12:24:26 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Fenics-arclength Documentation
============================================

**fenics-arclength** is a Python implementation of the arclength solver built on top of FEniCS. 
The Arc Length Method, sometimes referred to as the Riks Method, is a method used to solve solid mechanics problems with geometric nonlinearity. 
with complex equilibrium paths. 
The library aims to keep the usage as similiar to `FEniCS <https://fenicsproject.org/download/archive/>`_ (version 2019.1.0) to allow
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
`Instructions to install FEniCS (legacy) can be found here <https://fenicsproject.org/download/archive/>`_

Package Installation
####################
To use our arc-length solver, download and append our `repository <https://github.com/pprachas/fenics_arclength>`_ to the python path. Common methods to do this are:

* Add directory to ``PYTHONPATH``:

.. code-block:: bash

   export PYTHONPATH <path/to/fenics_arclength>:$PYTHONPATH

* Append directory to path in python script:

.. code-block:: python

   import sys
   sys.path.append('path/to/fenics_arclength')


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
