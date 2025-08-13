..
   Mustafa Aggul @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDomEigEst.API:

The SUNDomEigEstimator API
=============================

The SUNDomEigEst API defines several dominant eigenvalue estimation operations that enable
SUNDIALS packages to utilize this API.  These functions can be divided into three categories.
The first are the core dominant eigenvalue estimation functions.  The second consist of "set"
routines to supply the dominant eigenvalue estimator with functions provided by the SUNDIALS
packages and to modify estimator parameters.  The final group consists of "get" routines for
retrieving dominant eigenvalue estimation statistics.  All of these functions are defined in
the header file ``sundials/sundials_domeigestimator.h``.


.. _SUNDomEigEst.CoreFn:

SUNDomEigEstimator core functions
-----------------------------------------------------

The SUNDomEigEstimator base class provides two **utility** routines for implementers, 
:c:func:`SUNDomEigEst_NewEmpty` and :c:func:`SUNDomEigEst_FreeEmpty`.

Implementations of SUNDomEigEstimators must include a set of **required** functions: 
:c:func:`SUNDomEigEst_SetATimes` provides a :c:type:`SUNATimesFn` function pointer,
as well as a ``void*`` pointer to a data structure used by this routine,
:c:func:`SUNDomEig_Estimate` estimates the dominant eigenvalue,
and :c:func:`SUNDomEigEst_Destroy` destroys an estimator object.

.. c:function:: SUNDomEigEstimator SUNDomEigEst_NewEmpty(SUNContext sunctx)

   This function allocates a new ``SUNDomEigEstimator`` object and
   initializes its content pointer and the function pointers in the operations
   structure to ``NULL``.

   **Arguments:**

      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**

      If successful, this function returns a ``SUNDomEigEstimator`` object.
      If an error occurs when allocating the object, then this routine will
      return ``NULL``.


.. c:function:: SUNErrCode SUNDomEigEst_Initialize(SUNDomEigEstimator DEE)

   This *optional* function performs dominant eigenvalue estimator initialization (assuming that all  
   estimator-specific options have been set).

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. c:function:: SUNErrCode SUNDomEig_Estimate(SUNDomEigEstimator DEE, sunrealtype* lambdaR, sunrealtype* lambdaI)

   This *required* function estimates the dominant eigenvalue,
   :math:`\lambda_{\max} = \lambda_{R} + \lambda_{I}i` such that
   :math:`|\lambda| = \max\{|\lambda_i| : A \vec{v_i} = \lambda_i \vec{v_i}, \ \vec{v_i} \neq \vec{0} \}`.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *lambdaR* -- The real part of the dominant eigenvalue
      * *lambdaI* -- The imaginary part of the dominant eigenvalue

   **Return value:**

      `SUN_SUCCESS` for a successful call, or a relevant error code from
      :numref:`SUNDomEigEst.ErrorCodes` upon failure.


.. c:function:: SUNErrCode SUNDomEigEst_FreeEmpty(SUNDomEigEstimator DEE)

   This routine frees the ``SUNDomEigEstimator`` object, under the
   assumption that any implementation-specific data that was allocated
   within the underlying content structure has already been freed.
   It will additionally test whether the ops pointer is ``NULL``,
   and, if it is not, it will free it as well.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object

   **Return value:**

      A :c:type:`SUNErrCode`.


.. c:function:: SUNErrCode SUNDomEigEst_Destroy(SUNDomEigEstimator* DEEptr)

   Frees memory allocated by the dominant eigenvalue estimator.

   **Arguments:**

      * *DEEptr* -- a SUNDomEigEstimator object pointer.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEst_Destroy(&DEE);


.. _SUNDomEigEst.SetFn:

SUNDomEigEstimator "set" functions
-------------------------------------

The following functions supply dominant eigenvalue estimator modules with
functions defined by the SUNDIALS packages and modify estimator parameters.
Only the routine for setting the matrix-vector product routine is required.
Otherwise, all other set functions are optional.  SUNDomEigEst implementations
that do not provide the functionality for any optional routine should leave the corresponding
function pointer ``NULL`` instead of supplying a dummy routine.


.. c:function:: SUNErrCode SUNDomEigEst_SetATimes(SUNDomEigEstimator DEE, void* A_data, SUNATimesFn ATimes)

   This *required* function provides a :c:type:`SUNATimesFn` function pointer, as well as a ``void*`` pointer to a
   data structure used by this routine, to the dominant eigenvalue estimator object
   `DEE``.  SUNDIALS packages call this function to set the matrix-vector product function
   to either an estimator-provided difference-quotient via vector operations or a user-supplied
   estimator-specific routine.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *A_data* -- pointer to structure for ``ATimes``,
      * *ATimes* -- function pointer to perform :math:`Av` product.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. c:function:: SUNErrCode SUNDomEigEst_SetNumPreProcess(SUNDomEigEstimator DEE, int numpreprocess)

   This *optional* routine should set the number of "warm-up" matrix-vector multiplications,
   which is executed by :c:func:`SUNDomEig_Estimate` before each estimate.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *numpreprocess* -- the number of preprocessing iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.

      
   .. note:: When the DEE is used within LSRKStep (see :c:func:`LSRKStepSetDomEigEstimator`), 
      this number of warmup iterations will be overwritten after the first call to 
      `SUNDomEig_Estimate` (see :c:func:`LSRKSetNumDomEigEstPreprocessIters`).


.. c:function:: SUNErrCode SUNDomEigEst_SetRelTol(SUNDomEigEstimator DEE, sunrealtype rel_tol)

   This *optional* routine sets the estimator's :ref:`relative tolerance <pi_rel_tol>`.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *rel_tol* -- the requested eigenvalue accuracy.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. c:function:: SUNErrCode SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE, long int max_iters)

   This *optional* routine sets the maximum number of iterations.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *max_iters* -- the maximum number of iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. _SUNDomEigEst.GetFn:

SUNDomEigEstimator "get" functions
----------------------------------

The following functions allow SUNDIALS packages to retrieve results from a
dominant eigenvalue estimator.  *All routines are optional.*

.. c:function:: SUNErrCode SUNDomEigEst_GetRes(SUNDomEigEstimator DEE, sunrealtype* cur_res)

   This *optional* routine should return the final residual from
   the most-recent call to :c:func:`:SUNDomEigEst_Estimate`. 

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *cur_res* -- the current residual

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         sunrealtype cur_res;
         retval = SUNDomEigEst_GetRes(DEE, &cur_res);


.. c:function:: SUNErrCode SUNDomEigEst_GetNumIters(SUNDomEigEstimator DEE, long int* num_iters)

   This *optional* routine should return the number of estimator
   iterations performed in the most-recent call to :c:func:`:SUNDomEigEst_Estimate`.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *num_iters* -- the number of iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         long int num_iters;
         retval = SUNDomEigEst_GetNumIters(DEE, &num_iters);


.. c:function:: SUNErrCode SUNDomEigEst_GetMaxNumIters(SUNDomEigEstimator DEE, long int* max_niter)

   This *optional* routine should return the maximum number of iterations
   performed in any previous "estimator" call so far.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *max_niter* -- the maximum number of iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         long int max_niter;
         retval = SUNDomEigEst_GetMaxNumIters(DEE, &max_niter);


.. c:function:: SUNErrCode SUNDomEigEst_GetMinNumIters(SUNDomEigEstimator DEE, long int* min_niter)

   This *optional* routine should return the minimum number of iterations
   performed in any previous "estimator" call so far.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *min_niter* -- the minimum number of iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         long int min_niter;
         retval = SUNDomEigEst_GetMinNumIters(DEE, &min_niter);

.. c:function:: SUNErrCode SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE, long int* num_ATimes)

   This *optional* routine should return the number of calls to the :c:type:`SUNATimesFn` function.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *num_ATimes* -- the number of calls to the ``Atimes`` function.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         long int num_ATimes;
         retval = SUNDomEigEst_GetNumATimesCalls(DEE, &num_ATimes);


.. c:function:: SUNErrCode SUNDomEigEst_Write(SUNDomEigEstimator DEE, FILE* outfile)

   This *optional* routine prints the dominant eigenvalue estimator statistics
   to the output stream *outfile*.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object,
      * *outfile* -- the output stream.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. _SUNDomEigEst.SUNSuppliedFn:

Functions provided by SUNDIALS packages
---------------------------------------------

To interface with SUNDomEigEst modules, the SUNDIALS packages supply a routine
:c:type:`SUNATimesFn` for evaluating the matrix-vector product.  This package-provided
routine translates between the user-supplied ODE or DAE systems and the generic
dominant eigenvalue estimator API.  The function types
for these routines are defined in the header file ``sundials/sundials_iterative.h``.

.. _SUNDomEigEst.ReturnCodes:

SUNDomEigEstimator return codes
------------------------------------

The functions provided to SUNDomEigEst modules by each SUNDIALS package,
and functions within the SUNDIALS-provided SUNDomEigEst implementations,
utilize a common set of return codes, listed in :numref:`SUNDomEigEst.ErrorCodes`.


.. _SUNDomEigEst.ErrorCodes:
.. table:: SUNDomEigEst error codes
   :align: center

   +------------------------------------+-------+---------------------------------------------------+
   | Error code                         | Value | Meaning                                           |
   +====================================+=======+===================================================+
   | ``SUN_SUCCESS``                    | 0     | successful call or converged estimate             |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_BAD_NVECTOR``        | -9973 | bad NVector                                       |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_NULL_ATIMES``        | -9972 | the ``Atimes`` function ptr is ``NULL``           |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_ATIMES_FAIL_REC``    | -9971 | an unrecoverable failure occurred in the          |
   |                                    |       | ``ATimes`` routine                                |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_ATIMES_FAIL_UNREC``  | -9970 | a recoverable failure occurred in the             |
   |                                    |       | ``ATimes`` routine                                |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_NULL_HES``           | -9969 | the Hessenberg matrix is ``NULL``                 |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_NULL_MEM``           | -9968 | the DEE memory is ``NULL``                        |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_NULL_CONTENT``       | -9967 | the DEE content is ``NULL``                       |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_LAPACK_FAIL``        | -9966 | LAPACK ``_dgeev/_sgeev`` function failure         |
   |                                    |       |                                                   |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_NULL_ESTIMATE``      | -9965 | estimate function ptr is ``NULL``                 |
   |                                    |       |                                                   |
   +------------------------------------+-------+---------------------------------------------------+
   | ``SUN_ERR_DEE_NULL_FREE``          | -9964 | free function ptr is ``NULL``                     |
   |                                    |       |                                                   |
   +------------------------------------+-------+---------------------------------------------------+


.. _SUNDomEigEst.Generic:

The generic SUNDomEigEstimator module
-----------------------------------------

SUNDIALS packages interact with dominant eigenvalue estimator implementations through the
:c:type:`SUNDomEigEstimator` class.  A :c:type:`SUNDomEigEstimator` is a pointer to the
:c:struct:`_generic_SUNDomEigEstimator` structure:

.. c:type:: struct SUNDomEigEstimator_ *SUNDomEigEstimator

.. c:struct:: _generic_SUNDomEigEstimator

   The structure defining the SUNDIALS dominant eigenvalue estimator class.

   .. c:member:: void *content

      Pointer to the dominant eigenvalue estimator-specific member data

   .. c:member:: SUNDomEigEstimator_Ops ops

      A virtual table of dominant eigenvalue estimator operations provided by a specific
      implementation

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context

The virtual table structure is defined as

.. c:type:: struct SUNDomEigEstimator_Ops_ *SUNDomEigEstimator_Ops

.. c:struct:: _generic_SUNDomEigEstimator_Ops

   The structure defining :c:type:`SUNDomEigEstimator` operations.

   .. c:member:: SUNErrCode (*setatimes)(SUNDomEigEstimator, void*, SUNATimesFn)

      The function implementing :c:func:`SUNDomEigEst_SetATimes`

   .. c:member:: SUNErrCode (*setmaxiters)(SUNDomEigEstimator, int)

      The function implementing :c:func:`SUNDomEigEst_SetMaxIters`

   .. c:member:: SUNErrCode (*setnumpreprocess)(SUNDomEigEstimator, int)

      The function implementing :c:func:`SUNDomEigEst_SetNumPreProcess`

   .. c:member:: SUNErrCode (*settol)(SUNDomEigEstimator, sunrealtype)

      The function implementing :c:func:`SUNDomEigEst_SetRelTol`

   .. c:member:: SUNErrCode (*initialize)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEst_Initialize`

   .. c:member:: SUNErrCode (*estimate)(SUNDomEigEstimator, sunrealtype*, sunrealtype*)

      The function implementing :c:func:`SUNDomEig_Estimate`

   .. c:member:: sunrealtype (*getcurres)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEst_GetRes`

   .. c:member:: int (*getcurniters)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEst_GetNumIters`

   .. c:member:: int (*getmaxniters)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEst_GetMaxNumIters`

   .. c:member:: int (*getminniters)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEst_GetMinNumIters`

   .. c:member:: long int (*getnumatimescalls)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEst_GetNumATimesCalls`
      
   .. c:member:: SUNErrCode (*write)(SUNDomEigEstimator, FILE*)

      The function implementing :c:func:`SUNDomEigEst_Write`
      
   .. c:member:: SUNErrCode (*destroy)(SUNDomEigEstimator*)

      The function implementing :c:func:`SUNDomEigEst_Destroy`

The generic SUNDomEigEst class defines and implements the dominant eigenvalue estimator
operations defined in :numref:`SUNDomEigEst.CoreFn` -- :numref:`SUNDomEigEst.GetFn`.
These routines are in fact only wrappers to the dominant eigenvalue estimator operations
defined by a particular SUNDomEigEst implementation, which are accessed through
the *ops* field of the ``SUNDomEigEstimator`` structure.  To illustrate this
point we show below the implementation of a typical dominant eigenvalue estimator operation
from the ``SUNDomEigEstimator`` base class, namely :c:func:`SUNDomEigEst_Initialize`,
that initializes a ``SUNDomEigEstimator`` object for use after it has been
created and configured, and returns a flag denoting a successful or failed
operation:

.. code-block:: c

   SUNErrCode SUNDomEigEst_Initialize(SUNDomEigEstimator DEE)
   {
     return (DEE->ops->initialize(DEE));
   }


Additionally, a ``SUNDomEigEstimator`` implementation *may* do the following:

* Define and implement additional user-callable "set" routines
  acting on the ``SUNDomEigEstimator``, e.g., for setting various
  configuration options to tune the dominant eigenvalue estimator
  for a particular problem.

* Provide additional user-callable "get" routines acting on the
  ``SUNDomEigEstimator`` object, e.g., for returning various estimator
  statistics.
  

.. _SUNDomEigEst.Intended:

SUNDIALS modules SUNDomEigEstimator interface
==============================================

In :numref:`SUNDomEigEst.Intended.Usage`, we list the SUNDomEigEst module functions used
within SUNDIALS modules.  We emphasize that the user does not need to know
detailed usage of dominant eigenvalue estimator functions by a SUNDIALS module
in order to use a module.  The information is presented as an implementation detail for
the interested reader.

.. _SUNDomEigEst.Intended.Usage:
.. table:: List of SUNDomEigEst functions called by a SUNDIALS module dominant eigenvalue
           estimator interface.  Functions marked with "X" are required;
           functions marked with "O" are only called if they are non-``NULL`` and
           functions marked with "N/A" are not applicable in the ``SUNDomEigEstimator``
           implementation that is being used.
   :align: center

   +----------------------------------------------------+---------------------+---------------------+
   | Routine                                            |   POWER ITERATION   |  ARNOLDI ITERATION  |
   |                                                    |                     |                     |
   +====================================================+=====================+=====================+
   | :c:func:`SUNDomEigEst_SetATimes`                   |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetMaxIters`\ :sup:`1`       |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetNumPreProcess`            |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetRelTol`\ :sup:`1`         |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_Initialize`                  |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEig_Estimate`                       |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_GetRes`\ :sup:`2`            |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_GetNumIters`\ :sup:`3`       |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_GetMaxNumIters`\ :sup:`3`    |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_GetMinNumIters`\ :sup:`3`    |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_GetNumATimesCalls`           |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_Write`                       |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_Destroy`\ :sup:`4`           |                     |                     |
   +----------------------------------------------------+---------------------+---------------------+


Notes:

1. :c:func:`SUNDomEigEst_SetMaxIters()` and :c:func:`SUNDomEigEst_SetRelTol()` might or 
   might not be required depending on ``SUNDomEigEstimator`` implementation that is being used. 
   These flags must be left ``NULL`` if it is not applicable for an estimator.

2. Although :c:func:`SUNDomEigEst_GetRes()` is optional, if it is not
   implemented by the ``SUNDomEigEstimator`` then the interface will consider all
   estimates a being *exact*.

3. :c:func:`SUNDomEigEst_GetNumIters()`, :c:func:`SUNDomEigEst_GetMaxNumIters()`
   and :c:func:`SUNDomEigEst_GetMinNumIters()` are optional, if they are not
   implemented by the ``SUNDomEigEstimator`` then the interface will consider all
   estimates as requiring zero iterations.

4. Although the interface does not call :c:func:`SUNDomEigEst_Destroy()`
   directly, this routine should be available for users to call when
   cleaning up from a simulation.