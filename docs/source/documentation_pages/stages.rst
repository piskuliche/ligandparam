Available Parametrization Stages
=================================

The parametrization process is broken down into a series of stages, each of which is a subclass of the :doc:`./stages/abstractstage` class. 
Each stage is responsible for a specific part of the parametrization process, and can be easily added or removed from the pipleine. 
New stages can easily be added by subclassing the :doc:`./stages/abstractstage` class and implementing the required methods.



.. toctree::
   :maxdepth: 2
   :caption: Available Stages:


   ./stages/initialize.rst
   ./stages/charge.rst
   ./stages/gaussian.rst
   ./stages/resp.rst
   ./stages/parmchk.rst
   ./stages/leap.rst
   ./stages/build_system.rst
   ./stages/typematching.rst


.. toctree::
   :maxdepth: 2
   :caption: Development Stages:

   ./stages/abstractstage.rst
   ./stages/teststage.rst







