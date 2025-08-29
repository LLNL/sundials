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

Implementations of SUNDomEigEstimators must include a **required**
:c:func:`SUNDomEig_Estimate` function to estimate the dominant eigenvalue.

.. c:function:: SUNDomEigEstimator SUNDomEigEst_NewEmpty(SUNContext sunctx)

   This function allocates a new ``SUNDomEigEstimator`` object and
   initializes its content pointer and the function pointers in the operations
   structure to ``NULL``.

   **Arguments:**

      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`).

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
      * *lambdaR* -- The real part of the dominant eigenvalue.
      * *lambdaI* -- The imaginary part of the dominant eigenvalue.

   **Return value:**

      `SUN_SUCCESS` for a successful call, or a relevant error code from
      :c:type:`SUNErrCode` upon failure.


.. c:function:: SUNErrCode SUNDomEigEst_FreeEmpty(SUNDomEigEstimator DEE)

   This routine frees the ``SUNDomEigEstimator`` object, under the
   assumption that any implementation-specific data that was allocated
   within the underlying content structure has already been freed.
   It will additionally test whether the ops pointer is ``NULL``,
   and, if it is not, it will free it as well.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. c:function:: SUNErrCode SUNDomEigEst_Destroy(SUNDomEigEstimator* DEEptr)

   Frees memory allocated by the dominant eigenvalue estimator.

   **Arguments:**

      * *DEEptr* -- a SUNDomEigEstimator object pointer.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. _SUNDomEigEst.SetFn:

SUNDomEigEstimator "set" functions
-------------------------------------

The following functions supply dominant eigenvalue estimator modules with
functions defined by the SUNDIALS packages and modify estimator parameters.
When using the matrix-vector product routine provided by a SUNDIALS integration,
the ``SetATimes`` is required. Otherwise, all set functions are optional.
SUNDomEigEst implementations that do not provide the functionality for any
optional routine should leave the corresponding function pointer ``NULL``
instead of supplying a dummy routine.


.. c:function:: SUNErrCode SUNDomEigEst_SetATimes(SUNDomEigEstimator DEE, void* A_data, SUNATimesFn ATimes)

   This function provides a :c:type:`SUNATimesFn` function for performing
   matrix-vector products, as well as a ``void*`` pointer to a data structure
   used by this routine, to the dominant eigenvalue estimator. This function is
   *required* when using the matrix-vector product function provided by a
   SUNDIALS integrator, otherwise the function is *optional*.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *A_data* -- pointer to structure for ``ATimes``.
      * *ATimes* -- function pointer to perform :math:`Av` product.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. c:function:: SUNErrCode SUNDomEigEst_SetNumPreprocessIters(SUNDomEigEstimator DEE, int num_iters)

   This *optional* routine should set the number of preprocessing matrix-vector
   multiplications, performed at the beginning of each
   :c:func:`SUNDomEig_Estimate` evaluation.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *num_iters* -- the number of preprocessing iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.

   .. note::

      Prior to computing the dominant eigenvalue in :c:func:`SUNDomEig_Estimate`
      an implementation may perform ``num_iters`` power iterations on ``q`` to
      generate an improved initial guess.  Preprocessing iterations can help
      reduce some computational overhead, and may be useful if the initial guess
      ``q`` is not a good approximation of the dominant eigenvector.

      When the estimator is used in a time-dependent context, it is likely that
      the most-recent ``q`` will provide a suitable initial guess for subsequent
      calls to :c:func:`SUNDomEig_Estimate`. Thus, when the estimator is used
      with LSRKStep (see :c:func:`LSRKStepSetDomEigEstimator`), the initial
      value of ``num_iters`` should be set with
      :c:func:`LSRKStepSetNumDomEigEstInitPreprocessIters` while the number of
      preprocessing iterations for subsequent calls should be set with
      :c:func:`LSRKStepSetNumDomEigEstPreprocessIters`.

      Both the Arnodli and Power implementations provided with SUNDIALS use a
      default value of 100. This default value is particularly chosen to
      minimize the memory footprint by lowering the required ``kry_dim`` in the
      Arnoldi iteration, or reducing computational overhead when estimating with
      the power iteration. With either implementation, supplying a ``num_iters``
      argument that is :math:`< 0`, it will reset the value to the default.

.. c:function:: SUNErrCode SUNDomEigEst_SetRelTol(SUNDomEigEstimator DEE, sunrealtype rel_tol)

   This *optional* routine sets the estimator's :ref:`relative tolerance <pi_rel_tol>`.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *rel_tol* -- the requested eigenvalue accuracy.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. c:function:: SUNErrCode SUNDomEigEst_SetMaxIters(SUNDomEigEstimator DEE, long int max_iters)

   This *optional* routine sets the maximum number of iterations.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
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
   the most-recent call to :c:func:`SUNDomEig_Estimate`.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *cur_res* -- the residual.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         sunrealtype cur_res;
         retval = SUNDomEigEst_GetRes(DEE, &cur_res);


.. c:function:: SUNErrCode SUNDomEigEst_GetNumIters(SUNDomEigEstimator DEE, long int* num_iters)

   This *optional* routine should return the number of estimator
   iterations performed in the most-recent call to :c:func:`SUNDomEig_Estimate`.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *num_iters* -- the number of iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         long int num_iters;
         retval = SUNDomEigEst_GetNumIters(DEE, &num_iters);


.. c:function:: SUNErrCode SUNDomEigEst_GetNumATimesCalls(SUNDomEigEstimator DEE, long int* num_ATimes)

   This *optional* routine should return the number of calls to the :c:type:`SUNATimesFn` function.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *num_ATimes* -- the number of calls to the ``Atimes`` function.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         long int num_ATimes;
         retval = SUNDomEigEst_GetNumATimesCalls(DEE, &num_ATimes);


.. c:function:: SUNErrCode SUNDomEigEst_Write(SUNDomEigEstimator DEE, FILE* outfile)

   This *optional* routine prints the dominant eigenvalue estimator settings to
   the file pointer.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object.
      * *outfile* -- the output stream.

   **Return value:**

      A :c:type:`SUNErrCode`.


.. _SUNDomEigEst.SUNSuppliedFn:

Functions provided by SUNDIALS packages
---------------------------------------------

To interface with SUNDomEigEst modules, the SUNDIALS packages supply a
:c:type:`SUNATimesFn` function for evaluating the matrix-vector product. This
package-provided routine translates between the user-supplied ODE or DAE systems
and the generic dominant eigenvalue estimator API. The function types for these
routines are defined in the header file ``sundials/sundials_iterative.h``.


.. _SUNDomEigEst.Generic:

The generic SUNDomEigEstimator module
-----------------------------------------

SUNDIALS packages interact with dominant eigenvalue estimator implementations through the
:c:type:`SUNDomEigEstimator` class.  A :c:type:`SUNDomEigEstimator` is a pointer to the
:c:struct:`SUNDomEigEstimator_` structure:

.. c:type:: struct SUNDomEigEstimator_ *SUNDomEigEstimator

.. c:struct:: SUNDomEigEstimator_

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

.. c:struct:: SUNDomEigEstimator_Ops_

   The structure defining :c:type:`SUNDomEigEstimator` operations.

   .. c:member:: SUNErrCode (*setatimes)(SUNDomEigEstimator, void*, SUNATimesFn)

      The function implementing :c:func:`SUNDomEigEst_SetATimes`

   .. c:member:: SUNErrCode (*setmaxiters)(SUNDomEigEstimator, int)

      The function implementing :c:func:`SUNDomEigEst_SetMaxIters`

   .. c:member:: SUNErrCode (*setnumpreprocess)(SUNDomEigEstimator, int)

      The function implementing :c:func:`SUNDomEigEst_SetNumPreprocessIters`

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

In :numref:`SUNDomEigEst.Intended.Usage`, we list the SUNDomEigEst module
functions used within SUNDIALS packages. We emphasize that the user does not
need to know detailed usage of dominant eigenvalue estimator functions by a
SUNDIALS package in order to use it. The information is presented as an
implementation detail for the interested reader.

.. _SUNDomEigEst.Intended.Usage:
.. table:: List of SUNDomEigEst functions called by a SUNDIALS module dominant eigenvalue
           estimator interface.  Functions marked with "X" are required;
           functions marked with "O" are only called if they are non-``NULL`` and
           functions marked with "N/A" are not applicable in the ``SUNDomEigEstimator``
           implementation that is being used.
   :align: center

   +----------------------------------------------------+---------------------+---------------------+
   | Routine                                            |   Power Iteration   |  Arnoldi Iteration  |
   |                                                    |                     |                     |
   +====================================================+=====================+=====================+
   | :c:func:`SUNDomEigEst_SetATimes`                   |          X          |          X          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetMaxIters`\ :sup:`1`       |          O          |         N/A         |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_SetNumPreprocessIters`       |          O          |          O          |
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
   | :c:func:`SUNDomEigEst_GetNumATimesCalls`           |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_Write`                       |          O          |          O          |
   +----------------------------------------------------+---------------------+---------------------+
   | :c:func:`SUNDomEigEst_Destroy`\ :sup:`4`           |                     |                     |
   +----------------------------------------------------+---------------------+---------------------+


Notes:

1. :c:func:`SUNDomEigEst_SetMaxIters` and :c:func:`SUNDomEigEst_SetRelTol` might
   or might not be required depending on ``SUNDomEigEstimator`` implementation
   that is being used. These operations should be left as ``NULL`` if it is not
   applicable for an estimator.

2. Although :c:func:`SUNDomEigEst_GetRes` is optional, if it is not implemented
   by the ``SUNDomEigEstimator`` then the interface will consider all estimates
   a being *exact*.

3. :c:func:`SUNDomEigEst_GetNumIters` is optional, if it is not implemented by
   the ``SUNDomEigEstimator`` then the interface will consider all estimates as
   requiring zero iterations.

4. Although the interface does not call :c:func:`SUNDomEigEst_Destroy` directly,
   this routine should be available for users to call when cleaning up from a
   simulation.
