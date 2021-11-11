.. _Usage:

############
Using ARKODE
############

We present usage for ARKODE in five separate subsections.  The first three
focus on the C and C++ interfaces to each of ARKODE's time
stepping modules: :ref:`ARKStep <Usage.ARKStep>`, :ref:`ERKStep <Usage.ERKStep>`,
and :ref:`MRIStep <Usage.MRIStep>`.  Following these, we describe the
:ref:`modern Fortran 2003 interface <Usage.Fortran>` to these
time-stepping modules, and discuss the :ref:`use of GPU accelerators <Usage.GPU>`
for ARKODE simulations.


.. toctree::
   :maxdepth: 1

   ARKStep_c_interface/index.rst
   ERKStep_c_interface/index.rst
   MRIStep_c_interface/index.rst
   ARKode_f_interface/index.rst
   GPU
