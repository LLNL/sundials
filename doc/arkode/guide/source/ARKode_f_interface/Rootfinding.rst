..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _FInterface.RootFinding:

Usage of the FARKROOT interface to rootfinding
-----------------------------------------------

The FARKROOT interface package allows programs written in Fortran to
use the rootfinding feature of the ARKStep solver module. The
user-callable functions in FARKROOT, with the corresponding ARKStep
functions, are as follows: 

* :f:func:`FARKROOTINIT()` interfaces to :c:func:`ARKStepRootInit()`,

* :f:func:`FARKROOTINFO()` interfaces to
  :c:func:`ARKStepGetRootInfo()`, and 

* :f:func:`FARKROOTFREE()` interfaces to :c:func:`ARKStepRootInit()`,
  freeing memory by calling the initializer with no root functions.

Note that at this time, FARKROOT does not provide support to specify
the direction of zero-crossing that is to be monitored.  Instead, all
roots are considered.  However, the actual direction of zero-crossing
may be captured by the user through monitoring the sign of any
non-zero elements in the array *INFO* returned by
:f:func:`FARKROOTINFO()`. 

In order to use the rootfinding feature of ARKStep, after calling :f:func:`FARKMALLOC()` but prior to
calling :f:func:`FARKODE()`, the user must call
:f:func:`FARKROOTINIT()` to allocate and initialize memory for the FARKROOT module: 

.. f:subroutine:: FARKROOTINIT(NRTFN, IER)
   
   Initializes the Fortran interface to the FARKROOT module.
      
   **Arguments:** 
      * *NRTFN* (``int``, input) -- total number of root functions.
      * *IER* (``int``, output) -- return flag (0 success, -1 if
	ARKStep memory is ``NULL``, and -11 if a memory allocation
	error occurred).
      

If rootfinding is enabled, the user must specify the functions whose
roots are to be found.  These rootfinding functions should be
implemented in the user-supplied :f:func:`FARKROOTFN()` subroutine:

.. f:subroutine:: FARKROOTFN(T, Y, G, IPAR, RPAR, IER)
   
   User supplied function implementing the vector-valued function
   :math:`g(t,y)` such that the roots of the *NRTFN* components
   :math:`g_i(t,y)=0` are sought.
      
   **Arguments:** 
      * *T* (``realtype``, input) -- independent variable value :math:`t`.  
      * *Y* (``realtype``, input) -- dependent variable array :math:`y`. 
      * *G* (``realtype``, output) -- function value array :math:`g(t,y)`.  
      * *IPAR* (``long int``, input/output) -- integer user data
	array, the same as the array passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input/output) -- real-valued user data
	array, the same as the array passed to :f:func:`FARKMALLOC()`. 
      * *IER*  (``int``, output) -- return flag (0 success, :math:`<0`
	if error).
      

When making calls to :f:func:`FARKODE()` to solve the ODE system, the
occurrence of a root is flagged by the return value *IER = 2*.  In
that case, if *NRTFN > 1*, the functions :math:`g_i(t,y)` which were
found to have a root can be identified by calling the routine
:f:func:`FARKROOTINFO()`:

.. f:subroutine:: FARKROOTINFO(NRTFN, INFO, IER)
   
   Initializes the Fortran interface to the FARKROOT module.
      
   **Arguments:** 
      * *NRTFN* (``int``, input) -- total number of root functions.

      * *INFO* (``int``, input/output) -- array of length *NRTFN* with
	root information (must be allocated by the user).  For each
	index, *i = 1, ..., NRTFN*:
	
	* *INFO(i) = 1*  if :math:`g_i(t,y)` was found to have a root,
	  and :math:`g_i` is increasing.

	* *INFO(i) = -1*  if :math:`g_i(t,y)` was found to have a root,
	  and :math:`g_i` is decreasing.

	* *INFO(i) = 0*  otherwise.

      * *IER* (``int``, output) -- return flag (0 success, :math:`<0`
	if error).
      

The total number of calls made to the root function
:f:func:`FARKROOTFN()`, denoted *NGE*, can be obtained from
*IOUT(12)*.  If the FARKODE/ARKStep memory block is reinitialized to
solve a different problem via a call to :f:func:`FARKREINIT()`, then
the counter *NGE* is reset to zero. 

Lastly, to free the memory resources allocated by a prior call to
:f:func:`FARKROOTINIT()`, the user must make a call to
:f:func:`FARKROOTFREE()`:


.. f:subroutine:: FARKROOTFREE()
   
   Frees memory associated with the FARKODE rootfinding module.
      

