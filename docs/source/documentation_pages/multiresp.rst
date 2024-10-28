Multi-State RESP Fitting
=========================

The MultiResp subpackage contains a series of useful classes for handling Amber parm7 files, and for calculating multi-state RESP fitting to obtain
charges for MD simulation forcefields. The primary use case of this package is to be a combination of :class:`ligandparam.stages.gaussian.StageGaussianRotation` 
stage and the :class:`ligandparam.stages.resp.StageMultiRespFit` stage.


.. toctree::
   :maxdepth: 2
   :caption: Available Stages:


   ./multiresp/mdinutils.rst
   ./multiresp/parmhelper.rst
   ./multiresp/residueresp.rst
   ./multiresp/endstate.rst
   ./multiresp/intermolequiv.rst
   ./multiresp/functions.rst
   ./multiresp/respfunctions.rst




Attribution
-----------
This mini-core package is adapted from the (unreleased) parmutils package originally developed by Prof. Tim Giese at Rutgers University. These functions
and classes have been adapted to work with the ligand_param package, but have otherwise been left unchanged (except for documentation additions).





   






