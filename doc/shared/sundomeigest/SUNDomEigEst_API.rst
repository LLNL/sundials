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
SUNDIALS packages to utilize this API. These functions can be divided into
three categories. The first are the core dominant eigenvalue estimation functions. The
second consist of "set" routines to supply the dominant eigenvalue estimatitor with functions
provided by the SUNDIALS packages and to modify estimator parameters. The
final group consists of "get" routines for retrieving dominant eigenvalue estimation
statistics. All of these functions are defined in the header file
``sundials/sundials_domeigestimator.h``.


.. _SUNDomEigEst.CoreFn:

SUNDomEigEstimator core functions
-----------------------------------------------------

The core dominant eigenvalue estimatimator functions consist of several **required**
functions: :c:func:`SUNDomEigEstSetATimes` provides a :c:type:`SUNATimesFn` function pointer,
as well as a ``void*`` pointer to a data structure used by this routine,
:c:func:`SUNDomEigEstInitialize` initializes the estimator object once
all estimator-specific options have been set, and :c:func:`SUNDomEigEstimate` estimates
the dominant eigenvalue.

The remaining **optional** functions returns the estimator ID (:c:func:`SUNDomEigEstGetID`),
preprocess the estimator object to "warm-up" the estimator for a more appropriate initial vector
for iterations (:c:func:`SUNDomEigEstPreProcess`), computes Hessenberg matrix (when necessary)
(:c:func:`SUNDomEigEstComputeHess`), and destroy an estimator object (:c:func:`SUNDomEigEstFree`).


.. c:function:: SUNErrCode SUNDomEigEstInitialize(SUNDomEigEstimator DEE)

   *Required.*

   Performs dominant eigenvalue estimatimator initialization (assuming that all
   estimator-specific options have been set).

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstInitialize(DEE);


.. c:function:: SUNErrCode SUNDomEigEstComputeHess(SUNDomEigEstimator DEE)

   *Required* for some estimators (ARNOLDI) and *not applicable* for others (POWER)

   Performs Hessenberg matrix computation (assuming that the estimator is
   already initialized).

   **Return value:**

      Zero for a successful call, a positive value for a recoverable failure,
      and a negative value for an unrecoverable failure.  Ideally this should
      return one of the generic error codes listed in
      :numref:`SUNDomEigEst.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstComputeHess(DEE);


.. c:function:: SUNErrCode SUNDomEigEstimate(SUNDomEigEstimator DEE, suncomplextype* dom_eig)

   This *required* function estimates the dominant eigenvalue,
   :math:`\lambda_{\max} = \lambda` such that
   :math:`|\lambda| = \max\{|\lambda_i| : A \vec{v_i} = \lambda_i \vec{v_i}, \ \vec{v_i} \neq \vec{0} \}`.

   **Arguments:**

      * *DEE* -- a SUNDomEigEst object.
      * *dom_eig* -- a ``SUNMatrix`` object.

   **Return value:**

      Zero for a successful call, a positive value for a recoverable failure,
      and a negative value for an unrecoverable failure.  Ideally this should
      return one of the generic error codes listed in
      :numref:`SUNDomEigEst.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstimate(DEE, dom_eig);



.. c:function:: SUNErrCode SUNDomEigEstFree(SUNDomEigEstimator DEE)

   Frees memory allocated by the dominant eigenvalue estimatimator.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstFree(DEE);


.. _SUNDomEigEst.SetFn:

SUNDomEigEstimator "set" functions
-------------------------------------

The following functions supply dominant eigenvalue estimatimator modules with
functions defined by the SUNDIALS packages and modify estimator parameters.
Only the routine for setting the matrix-vector product routine is required.
Otherwise, all other set functions are optional. SUNDomEigEst implementations
that do not provide the functionality for any optional routine should leave the corresponding
function pointer ``NULL`` instead of supplying a dummy routine.


.. c:function:: SUNErrCode SUNDomEigEstSetATimes(SUNDomEigEstimator DEE, void* A_data, SUNATimesFn ATimes)

   *Required.*

   Provides a :c:type:`SUNATimesFn` function pointer, as well as a ``void*``
   pointer to a data structure used by this routine, to the dominant eigenvalue estimatimator object
   *DEE*.  SUNDIALS packages call this function to set the matrix-vector product function to either
   an estimator-provided difference-quotient via vector operations or a user-supplied
   estimator-specific routine.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstSetATimes(DEE, A_data, ATimes);


.. c:function:: SUNErrCode SUNDomEigEstSetMaxPowerIter(SUNDomEigEstimator DEE, sunindextype max_powiter)

   This *optional* routine sets the number of max power iterations.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstSetMaxPowerIter(DEE, max_powiter);


.. c:function:: SUNErrCode SUNDomEigEstSetTol(SUNDomEigEstimator DEE, sunrealtype tol)

   This *optional* routine sets the estimator tolerance.

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstSetTol(DEE, tol);


.. c:function:: SUNErrCode SUNDomEigEstPreProcess(SUNDomEigEstimator DEE)

   This *optional* routine executes the "warm-up" matrix-vector multiplications,
   whose number is set by (:c:func:`SUNDomEigEstSetNumPreProcess`).

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstPreProcess(DEE);


.. c:function:: SUNErrCode SUNDomEigEstSetNumPreProcess(SUNDomEigEstimator DEE, sunindextype numofperprocess)

   This *optional* routine should set the number of "warm-up" matrix-vector multiplications,
   which then should be executed by (:c:func:`SUNDomEigEstPreProcess`).

   **Return value:**

      A :c:type:`SUNErrCode`.

   **Usage:**

      .. code-block:: c

         retval = SUNDomEigEstSetNumPreProcess(DEE, numofperprocess);


.. _SUNDomEigEst.GetFn:

SUNDomEigEstimator "get" functions
----------------------------------

The following functions allow SUNDIALS packages to retrieve results from a
dominant eigenvalue estimatimate.  *All routines are optional.*


.. c:function:: SUNDomEigEstimator_ID SUNDomEigEstGetID(SUNDomEigEstimator DEE)

   This *optional* routine returns a non-negative estimator identifier (of type ``int``)
   for the dominant eigenvalue estimator *DEE*.

   **Return value:**

      Non-negative estimator identifier (of type ``int``), defined by the
      enumeration ``SUNDomEigEstimator_ID``, with values shown in
      :numref:`SUNDomEigEst.API.IDs` and defined in the ``sundials_domeigestimator.h``
      header file.

   **Usage:**

      .. code-block:: c

         id = SUNDomEigEstGetID(DEE);

   .. note::

      It is recommended that a user-supplied ``SUNDomEigEstimator`` return the
      ``SUNDSOMEIGESTIMATOR_CUSTOM`` identifier.


.. c:function:: int SUNDomEigEstNumIters(SUNDomEigEstimator DEE)

   This *optional* routine should return the number of estimator
   iterations performed in the most-recent "estimator" call.

   **Usage:**

      .. code-block:: c

         its = SUNDomEigEstNumIters(DEE);


.. c:function:: sunrealtype SUNDomEigEstRes(SUNDomEigEstimator DEE)

   This *optional* routine should return the final residual from
   the most-recent "estimator" call.

   **Usage:**

      .. code-block:: c

         res = SUNDomEigEstRes(DEE);


.. _SUNDomEigEst.SUNSuppliedFn:

Functions provided by SUNDIALS packages
---------------------------------------------

To interface with SUNDomEigEst modules, the SUNDIALS packages supply
a routine for evaluating the matrix-vector product.  This package-provided
routine translate between the user-supplied ODE, DAE, or linear and nonlinear
systems and the generic dominant eigenvalue estimatimator API. The
function types for these routines are defined in the header file
``sundials/sundials_iterative.h``, and are described below.


.. c:type:: int (*SUNATimesFn)(void *A_data, N_Vector v, N_Vector z)

   Computes the action of a matrix on a vector, performing the
   operation :math:`z \gets Av`.  Memory for *z* will already be
   allocated prior to calling this function.  The parameter
   *A_data* is a pointer to any information about :math:`A` which
   the function needs in order to do its job. The vector :math:`v`
   should be left unchanged.

   **Return value:**

      Zero for a successful call, and non-zero upon failure.


.. _SUNDomEigEst.ReturnCodes:

SUNDomEigEstimator return codes
------------------------------------

The functions provided to SUNDomEigEst modules by each SUNDIALS package,
and functions within the SUNDIALS-provided SUNDomEigEst implementations,
utilize a common set of return codes, listed in
:numref:`SUNDomEigEst.ErrorCodes`.  These adhere to a common pattern:

* 0 indicates success
* a positive value corresponds to a recoverable failure, and
* a negative value indicates a non-recoverable failure.

Aside from this pattern, the actual values of each error code
provide additional information to the user in case of an estimatitor
failure.

TODO:Add the right list here!

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
   | ``SUN_ERR_DEE_NULL_ATIMES``        | -9972 | the ``Atimes`` function is ``NULL``               |
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
   | ``SUN_ERR_DEE_LAPACK_FAIL``        | -9966 | LAPACK ``_dgeev`` function failure                |
   |                                    |       |                                                   |
   +------------------------------------+-------+---------------------------------------------------+



.. _SUNDomEigEst.Generic:

The generic SUNDomEigEstimator module
-----------------------------------------

SUNDIALS packages interact with dominant eigenvalue estimatitor implementations through the
:c:type:`SUNDomEigEstimator` class. A :c:type:`SUNDomEigEstimator` is a pointer to the
:c:struct:`_generic_SUNDomEigEstimator` structure:

.. c:type:: struct _generic_SUNDomEigEstimator *SUNDomEigEstimator

.. c:struct:: _generic_SUNDomEigEstimator

   The structure defining the SUNDIALS dominant eigenvalue estimatitor class.

   .. c:member:: void *content

      Pointer to the dominant eigenvalue estimatitor-specific member data

   .. c:member:: SUNDomEigEstimator_Ops ops

      A virtual table of dominant eigenvalue estimatitor operations provided by a specific
      implementation

   .. c:member:: SUNContext sunctx

      The SUNDIALS simulation context

The virtual table structure is defined as

.. c:type:: struct _generic_SUNDomEigEstimator_Ops *SUNDomEigEstimator_Ops

.. c:struct:: _generic_SUNDomEigEstimator_Ops

   The structure defining :c:type:`SUNDomEigEstimator` operations.

   .. c:member:: SUNDomEigEstimator_ID (*getid)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEstGetID`

   .. c:member:: SUNErrCode (*setatimes)(SUNDomEigEstimator, void*, SUNATimesFn)

      The function implementing :c:func:`SUNDomEigEstSetATimes`

   .. c:member:: SUNErrCode (*setmaxpoweriter)(SUNDomEigEstimator, sunindextype)

      The function implementing :c:func:`SUNDomEigEstSetMaxPowerIter`

   .. c:member:: SUNErrCode (*settol)(SUNDomEigEstimator, sunrealtype)

      The function implementing :c:func:`SUNDomEigEstSetTol`

   .. c:member:: SUNErrCode (*setnumofperprocess)(SUNDomEigEstimator, sunindextype)

      The function implementing :c:func:`SUNDomEigEstSetNumPreProcess`

   .. c:member:: SUNErrCode (*initialize)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEstInitialize`

   .. c:member:: SUNErrCode (*preprocess)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEstPreProcess`

   .. c:member:: SUNErrCode (*computehess)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEstComputeHess`

   .. c:member:: SUNErrCode (*estimate)(SUNDomEigEstimator, suncomplextype*)

      The function implementing :c:func:`SUNDomEigEstimate`

   .. c:member:: sunindextype (*getnumofiters)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEstNumIters`

   .. c:member:: sunrealtype (*res)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEstRes`

   .. c:member:: SUNErrCode (*free)(SUNDomEigEstimator)

      The function implementing :c:func:`SUNDomEigEstFree`

The generic SUNDomEigEst class defines and implements the dominant eigenvalue estimatitor
operations defined in :numref:`SUNDomEigEst.CoreFn` -- :numref:`SUNDomEigEst.GetFn`.
These routines are in fact only wrappers to the dominant eigenvalue estimatitor operations
defined by a particular SUNDomEigEst implementation, which are accessed through
the *ops* field of the ``SUNDomEigEstimator`` structure.  To illustrate this
point we show below the implementation of a typical dominant eigenvalue estimatitor operation
from the ``SUNDomEigEstimator`` base class, namely :c:func:`SUNDomEigEstInitialize`,
that initializes a ``SUNDomEigEstimator`` object for use after it has been
created and configured, and returns a flag denoting a successful or failed
operation:

.. code-block:: c

   SUNErrCode SUNDomEigEstInitialize(SUNDomEigEstimator DEE)
   {
     return ((SUNErrCode) DEE->ops->initialize(DEE));
   }


.. _SUNDomEigEst.API.Custom:

Implementing a custom SUNDomEigEstimator module
--------------------------------------------------

A particular implementation of the ``SUNDomEigEstimator`` module must:

* Specify the *content* field of the SUNDomEigEst module.

* Define and implement the required dominant eigenvalue estimatitor operations.

  .. note::

     The names of these routines should be unique to that
     implementation in order to permit using more than one
     SUNDomEigEst module (each with different ``SUNDomEigEstimator``
     internal data representations) in the same code.

* Define and implement user-callable constructor and destructor
  routines to create and free a ``SUNDomEigEstimator`` with
  the new *content* field and with *ops* pointing to the
  new dominant eigenvalue estimatitor operations.

We note that the function pointers for all unsupported optional
routines should be set to ``NULL`` in the *ops* structure.  This
allows the SUNDIALS package that is using the SUNDomEigEst object
to know whether the associated functionality is supported.

To aid in the creation of custom ``SUNDomEigEstimator`` modules the generic
``SUNDomEigEstimator`` module provides the utility function
:c:func:`SUNDomEigEstNewEmpty`. When used in custom ``SUNDomEigEstimator``
constructors this function will ease the introduction of any new optional dominant
eigenvalue estimatitor operations to the ``SUNDomEigEstimator`` API by ensuring that only required
operations need to be set.

.. c:function:: SUNDomEigEstimator SUNDomEigEstNewEmpty(SUNContext sunctx)

   This function allocates a new generic ``SUNDomEigEstimator`` object and
   initializes its content pointer and the function pointers in the operations
   structure to ``NULL``.

   **Return value:**

      If successful, this function returns a ``SUNDomEigEstimator`` object.
      If an error occurs when allocating the object, then this routine will
      return ``NULL``.

.. c:function:: void SUNDomEigEstFreeEmpty(SUNDomEigEstimator DEE)

   This routine frees the generic ``SUNDomEigEstimator`` object, under the
   assumption that any implementation-specific data that was allocated
   within the underlying content structure has already been freed.
   It will additionally test whether the ops pointer is ``NULL``,
   and, if it is not, it will free it as well.

   **Arguments:**

      * *DEE* -- a SUNDomEigEstimator object


Additionally, a ``SUNDomEigEstimator`` implementation *may* do the following:

* Define and implement additional user-callable "set" routines
  acting on the ``SUNDomEigEstimator``, e.g., for setting various
  configuration options to tune the dominant eigenvalue estimatitor
  for a particular problem.

* Provide additional user-callable "get" routines acting on the
  ``SUNDomEigEstimator`` object, e.g., for returning various estimator
  statistics.


.. c:enum:: SUNDomEigEstimator_ID

   Each SUNDomEigEst implementation included in SUNDIALS has a unique identifier
   specified in enumeration and shown in :numref:`SUNDomEigEst.API.IDs`. It is
   recommended that a user-supplied SUNDomEigEst implementation use the
   ``SUNDSOMEIGESTIMATOR_CUSTOM`` identifier.

.. _SUNDomEigEst.API.IDs:
.. table:: Identifiers associated with :c:type:`SUNDomEigEstimator`
           modules supplied with SUNDIALS
   :align: center

   ==================================  =======================================================  ============
   SUNDomEigEst ID                        Dominant eigenvalue estimatitor type                    ID Value
   ==================================  =======================================================  ============
   SUNDSOMEIGESTIMATOR_POWER           Power Iteration (internal)                                 0
   SUNDSOMEIGESTIMATOR_ARNOLDI         Arnoldi Iteration (internal)                               1
   SUNDSOMEIGESTIMATOR_CUSTOM          User-provided custom dominant eigenvalue estimatitor      15
   ==================================  =======================================================  ============