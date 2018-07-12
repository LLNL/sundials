..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _ARKStep_CInterface.UserSupplied:

User-supplied functions
=============================

The user-supplied functions for ARKStep consist of:

* at least one function defining the ODE (required),

* a function that handles error and warning messages (optional),

* a function that provides the error weight vector (optional),

* a function that provides the residual weight vector (optional),

* a function that handles adaptive time step error control (optional),

* a function that handles explicit time step stability (optional),

* a function that defines the root-finding problem(s) to solve
  (optional),

* one or two functions that provide Jacobian-related information for
  the linear solver, if a Newton-based nonlinear iteration is chosen
  (optional),

* one or two functions that define the preconditioner for use in any
  of the Krylov iterative algorithms, if a Newton-based nonlinear
  iteration and iterative linear solver are chosen (optional), and

* if the problem involves a non-identity mass matrix :math:`M\ne I`:

  * one or two functions that provide mass-matrix-related information
    for the linear and mass matrix solvers (required),

  * one or two functions that define the mass matrix preconditioner
    for use in an iterative mass matrix solver is chosen (optional), and

* a function that handles vector resizing operations, if the
  underlying vector structure supports resizing (as opposed to
  deletion/recreation), and if the user plans to call
  :c:func:`ARKStepResize()` (optional).




.. _ARKStep_CInterface.ODERHS:

ODE right-hand side
-----------------------------

The user must supply at least one function of type :c:type:`ARKRhsFn` to
specify the explicit and/or implicit portions of the ODE system:


.. c:type:: typedef int (*ARKRhsFn)(realtype t, N_Vector y, N_Vector ydot, void* user_data)

   These functions compute the ODE right-hand side for a given
   value of the independent variable :math:`t` and state vector :math:`y`.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector, :math:`y(t)`.
      * *ydot* -- the output vector that forms a portion of the ODE RHS :math:`f_E(t,y) + f_I(t,y)`.
      * *user_data* -- the `user_data` pointer that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   An *ARKRhsFn* should return 0 if successful, a positive value if a
   recoverable error occurred (in which case ARKStep will attempt to
   correct), or a negative value if it failed unrecoverably (in which
   case the integration is halted and *ARK_RHSFUNC_FAIL* is returned).

   **Notes:** Allocation of memory for `ydot` is handled within the
   ARKStep module.  A recoverable failure error return from the *ARKRhsFn* is
   typically used to flag a value of the dependent variable :math:`y`
   that is "illegal" in some way (e.g., negative where only a
   nonnegative value is physically meaningful).  If such a return is
   made, ARKStep will attempt to recover (possibly repeating the
   nonlinear iteration, or reducing the step size) in order to avoid
   this recoverable error return.  There are some situations in which
   recovery is not possible even if the right-hand side function
   returns a recoverable error flag.  One is when this occurs at the
   very first call to the *ARKRhsFn* (in which case
   ARKStep returns *ARK_FIRST_RHSFUNC_ERR*).  Another is when a
   recoverable error is reported by *ARKRhsFn* after the integrator
   completes a successful stage, in which case ARKStep returns
   *ARK_UNREC_RHSFUNC_ERR*).




.. _ARKStep_CInterface.ErrorHandler:

Error message handler function
--------------------------------------

As an alternative to the default behavior of directing error and
warning messages to the file pointed to by `errfp` (see
:c:func:`ARKStepSetErrFile()`), the user may provide a function of type
:c:type:`ARKErrHandlerFn` to process any such messages.



.. c:type:: typedef void (*ARKErrHandlerFn)(int error_code, const char* module, const char* function, char* msg, void* user_data)

   This function processes error and warning messages from
   ARKStep and is sub-modules.

   **Arguments:**
      * *error_code* -- the error code.
      * *module* -- the name of the ARKStep module reporting the error.
      * *function* -- the name of the function in which the error occurred.
      * *msg* -- the error message.
      * *user_data* -- a pointer to user data, the same as the
        *eh_data* parameter that was passed to :c:func:`ARKStepSetErrHandlerFn()`.

   **Return value:**
   An *ARKErrHandlerFn* function has no return value.

   **Notes:** *error_code* is negative for errors and positive
   (*ARK_WARNING*) for warnings.  If a function that returns a
   pointer to memory encounters an error, it sets *error_code* to
   0.




.. _ARKStep_CInterface.ErrorWeight:

Error weight function
--------------------------------------

As an alternative to providing the relative and absolute tolerances,
the user may provide a function of type :c:type:`ARKEwtFn` to compute a
vector *ewt* containing the weights in the WRMS norm
:math:`\|v\|_{WRMS} = \left(\frac{1}{n} \sum_{i=1}^n \left(ewt_i\; v_i\right)^2
\right)^{1/2}`.  These weights will be used in place of those defined
in the section :ref:`Mathematics.Error.Norm`.



.. c:type:: typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void* user_data)

   This function computes the WRMS error weights for the vector
   :math:`y`.

   **Arguments:**
      * *y* -- the dependent variable vector at which the
        weight vector is to be computed.
      * *ewt* -- the output vector containing the error weights.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   An *ARKEwtFn* function must return 0 if it
   successfully set the error weights, and -1 otherwise.

   **Notes:** Allocation of memory for *ewt* is handled within ARKStep.

   The error weight vector must have all components positive.  It is
   the user's responsibility to perform this test and return -1 if it
   is not satisfied.



.. _ARKStep_CInterface.ResidualWeight:

Residual weight function
--------------------------------------

As an alternative to providing the scalar or vector absolute residual
tolerances (when the IVP units differ from the solution units), the
user may provide a function of type :c:type:`ARKRwtFn` to compute a
vector *rwt* containing the weights in the WRMS norm
:math:`\|v\|_{WRMS} = \left(\frac{1}{n} \sum_{i=1}^n \left(rwt_i\; v_i\right)^2
\right)^{1/2}`.  These weights will be used in place of those defined
in the section :ref:`Mathematics.Error.Norm`.



.. c:type:: typedef int (*ARKRwtFn)(N_Vector y, N_Vector rwt, void* user_data)

   This function computes the WRMS residual weights for the vector
   :math:`y`.

   **Arguments:**
      * *y* -- the dependent variable vector at which the
        weight vector is to be computed.
      * *rwt* -- the output vector containing the residual weights.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   An *ARKRwtFn* function must return 0 if it
   successfully set the residual weights, and -1 otherwise.

   **Notes:** Allocation of memory for *rwt* is handled within ARKStep.

   The residual weight vector must have all components positive.  It is
   the user's responsibility to perform this test and return -1 if it
   is not satisfied.



.. _ARKStep_CInterface.AdaptivityFn:

Time step adaptivity function
--------------------------------------

As an alternative to using one of the built-in time step adaptivity
methods for controlling solution error, the user may provide a
function of type :c:type:`ARKAdaptFn` to compute a target step size
:math:`h` for the next integration step.  These steps should be chosen
as the maximum value such that the error estimates remain below 1.



.. c:type:: typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h1, realtype h2, realtype h3, realtype e1, realtype e2, realtype e3, int q, realtype* hnew, void* user_data)

   This function implements a time step adaptivity algorithm
   that chooses :math:`h` satisfying the error tolerances.

   **Arguments:**
      * *y* -- the current value of the dependent variable vector, :math:`y(t)`.
      * *t* -- the current value of the independent variable.
      * *h1* -- the current step size, :math:`t_m - t_{m-1}`.
      * *h2* -- the previous step size, :math:`t_{m-1} - t_{m-2}`.
      * *h3* -- the step size :math:`t_{m-2}-t_{m-3}`.
      * *e1* -- the error estimate from the current step, :math:`m`.
      * *e2* -- the error estimate from the previous step, :math:`m-1`.
      * *e3* -- the error estimate from the step :math:`m-2`.
      * *q* -- the global order of accuracy for the method.
      * *hnew* -- the output value of the next step size.
      * *user_data* -- a pointer to user data, the same as the
        *h_data* parameter that was passed to :c:func:`ARKStepSetAdaptivityFn()`.

   **Return value:**
   An *ARKAdaptFn* function should return 0 if it
   successfuly set the next step size, and a non-zero value otherwise.




.. _ARKStep_CInterface.StabilityFn:

Explicit stability function
--------------------------------------

A user may supply a function to predict the maximum stable step size
for the explicit portion of the ImEx system, :math:`f_E(t,y)`.  While
the accuracy-based time step adaptivity algorithms may be sufficient
for retaining a stable solution to the ODE system, these may be
inefficient if :math:`f_E(t,y)` contains moderately stiff terms.  In
this scenario, a user may provide a function of type :c:type:`ARKExpStabFn`
to provide this stability information to ARKStep.  This function
must set the scalar step size satisfying the stability restriction for
the upcoming time step.  This value will subsequently be bounded by
the user-supplied values for the minimum and maximum allowed time
step, and the accuracy-based time step.



.. c:type:: typedef int (*ARKExpStabFn)(N_Vector y, realtype t, realtype* hstab, void* user_data)

   This function predicts the maximum stable step size for the
   explicit portions of the ImEx ODE system.

   **Arguments:**
      * *y* -- the current value of the dependent variable vector, :math:`y(t)`.
      * *t* -- the current value of the independent variable
      * *hstab* -- the output value with the absolute value of the
 	maximum stable step size.
      * *user_data* -- a pointer to user data, the same as the
        *estab_data* parameter that was passed to :c:func:`ARKStepSetStabilityFn()`.

   **Return value:**
   An *ARKExpStabFn* function should return 0 if it
   successfully set the upcoming stable step size, and a non-zero
   value otherwise.

   **Notes:**  If this function is not supplied, or if it returns
   *hstab* :math:`\le 0.0`, then ARKStep will assume that there is no explicit
   stability restriction on the time step size.



.. _ARKStep_CInterface.RootfindingFn:

Rootfinding function
--------------------------------------

If a rootfinding problem is to be solved during the integration of the
ODE system, the user must supply a function of type :c:type:`ARKRootFn`.



.. c:type:: typedef int (*ARKRootFn)(realtype t, N_Vector y, realtype* gout, void* user_data)

   This function implements a vector-valued function
   :math:`g(t,y)` such that the roots of the *nrtfn* components
   :math:`g_i(t,y)` are sought.

   **Arguments:**
      * *t* -- the current value of the independent variable
      * *y* -- the current value of the dependent variable vector, :math:`y(t)`.
      * *gout* -- the output array, of length *nrtfn*, with components :math:`g_i(t,y)`.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   An *ARKRootFn* function should return 0 if successful
   or a non-zero value if an error occurred (in which case the
   integration is halted and ARKStep returns *ARK_RTFUNC_FAIL*).

   **Notes:** Allocation of memory for *gout* is handled within ARKStep.



.. _ARKStep_CInterface.JacobianFn:

Jacobian information (direct method Jacobian)
--------------------------------------------------------------

If the direct linear solver interface is used (i.e.,
:c:func:`ARKDlsSetLinearSolver()` is called in the section
:ref:`ARKStep_CInterface.Skeleton`), the user may provide a function of type
:c:type:`ARKDlsJacFn` to provide the Jacobian approximation.



.. c:type:: typedef int (*ARKDlsJacFn)(realtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the Jacobian matrix :math:`J =
   \frac{\partial f_I}{\partial y}` (or an approximation to it).

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector, namely
        the predicted value of :math:`y(t)`.
      * *fy* -- the current value of the vector :math:`f_I(t,y)`.
      * *Jac* -- the output Jacobian matrix.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.
      * *tmp1*, *tmp2*, *tmp3* -- pointers to memory allocated to
        variables of type ``N_Vector`` which can be used by an
        ARKDlsacFn as temporary storage or work space.

   **Return value:**
   An *ARKDlsJacFn* function should return 0 if successful, a positive
   value if a recoverable error occurred (in which case ARKStep will
   attempt to correct, while ARKDLS sets *last_flag* to
   *ARKDLS_JACFUNC_RECVR*), or a negative value if it failed
   unrecoverably (in which case the integration is halted,
   :c:func:`ARKStepEvolve()` returns *ARK_LSETUP_FAIL* and ARKDLS sets
   *last_flag* to *ARKDLS_JACFUNC_UNRECVR*).

   **Notes:** Information regarding the structure of the specific
   ``SUNMatrix`` structure (e.g.~number of rows, upper/lower
   bandwidth, sparsity type) may be obtained through using the
   implementation-specific ``SUNMatrix`` interface functions
   (see the section :ref:`SUNMatrix` for details).

   Prior to calling the user-supplied Jacobian function, the Jacobian
   matrix :math:`J(t,y)` is zeroed out, so only nonzero elements need
   to be loaded into *Jac*.

   If the user's :c:type:`ARKDlsJacFn` function uses difference
   quotient approximations, then it may need to access quantities not
   in the argument list.  These include the current step size, the
   error weights, etc.  To obtain these, the user will need to add a
   pointer to the ``ark_mem`` structure to their ``user_data``, and
   then use the ARKStepGet* functions listed in
   :ref:`ARKStep_CInterface.OptionalOutputs`. The unit roundoff can be
   accessed as ``UNIT_ROUNDOFF``, which is defined in the header
   file ``sundials_types.h``.

   **dense**:

   A user-supplied dense Jacobian function must load the
   *N* by *N* dense matrix *Jac* with an approximation to the Jacobian
   matrix :math:`J(t,y)` at the point :math:`(t,y)`. The accessor
   macros ``SM_ELEMENT_D`` and ``SM_COLUMN_D`` allow the user to read
   and write dense matrix elements without making explicit references
   to the underlying representation of the SUNMATRIX_DENSE type.
   ``SM_ELEMENT_D(J, i, j)`` references the ``(i,j)``-th element of
   the dense matrix ``J`` (for ``i``, ``j`` between 0 and
   N-1). This macro is meant for small problems for which
   efficiency of access is not a major concern. Thus, in terms of the
   indices :math:`m` and :math:`n` ranging from 1 to *N*, the
   Jacobian element :math:`J_{m,n}` can be set using the statement
   ``SM_ELEMENT_D(J, m-1, n-1) =`` :math:`J_{m,n}`.  Alternatively,
   ``SM_COLUMN_D(J, j)`` returns a pointer to the first element of the
   ``j``-th column of ``J`` (for ``j`` ranging from 0 to `N`-1),
   and the elements of the ``j``-th column can then be accessed using
   ordinary array indexing. Consequently, :math:`J_{m,n}` can be
   loaded using the statements
   ``col_n = SM_COLUMN_D(J, n-1); col_n[m-1] =`` :math:`J_{m,n}`.
   For large problems, it is more efficient to use ``SM_COLUMN_D``
   than to use ``SM_ELEMENT_D``.  Note that both of these macros
   number rows and columns starting from 0.  The SUNMATRIX_DENSE type
   and accessor macros are documented in section
   :ref:`SUNMatrix_Dense`.

   **band**:

   A user-supplied banded Jacobian function must load the band
   matrix *Jac* with the elements of the Jacobian
   :math:`J(t,y)` at the point :math:`(t,y)`. The accessor macros
   ``SM_ELEMENT_B``, ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B``
   allow the user to read and write band matrix elements without
   making specific references to the underlying representation of the
   SUNMATRIX_BAND type.  ``SM_ELEMENT_B(J, i, j)`` references the
   ``(i,j)``-th element of the band matrix ``J``, counting
   from 0. This macro is meant for use in small problems for
   which efficiency of access is not a major concern. Thus, in terms
   of the indices :math:`m` and :math:`n` ranging from 1 to *N* with
   :math:`(m, n)` within the band defined by *mupper* and
   *mlower*, the Jacobian element :math:`J_{m,n}` can be loaded
   using the statement ``SM_ELEMENT_B(J, m-1, n-1)`` :math:`=
   J_{m,n}`. The elements within the band are those with *-mupper*
   :math:`\le m-n \le` *mlower*.  Alternatively, ``SM_COLUMN_B(J, j)``
   returns a pointer to the diagonal element of the ``j``-th column of
   ``J``, and if we assign this address to ``realtype *col_j``, then
   the ``i``-th element of the ``j``-th column is given by
   ``SM_COLUMN_ELEMENT_B(col_j, i, j)``, counting from 0. Thus, for
   :math:`(m,n)` within the band, :math:`J_{m,n}` can be loaded by
   setting ``col_n = SM_COLUMN_B(J, n-1); SM_COLUMN_ELEMENT_B(col_n, m-1,
   n-1)`` :math:`= J_{m,n}` . The elements of the ``j``-th column can
   also be accessed via ordinary array indexing, but this approach
   requires knowledge of the underlying storage for a band matrix of
   type SUNMATRIX_BAND. The array ``col_n`` can be indexed from
   *-mupper* to *mlower*. For large problems, it is more efficient
   to use ``SM_COLUMN_B`` and ``SM_COLUMN_ELEMENT_B`` than to use the
   ``SM_ELEMENT_B`` macro. As in the dense case, these macros all
   number rows and columns starting from 0. The SUNMATRIX_BAND type
   and accessor macros are documented in section :ref:`SUNMatrix_Band`.

   **sparse**:

   A user-supplied sparse Jacobian function must load the
   compressed-sparse-column (CSC) or compressed-sparse-row (CSR)
   matrix *Jac* with an approximation to the Jacobian matrix
   :math:`J(t,y)` at the point :math:`(t,y)`.  Storage for *Jac*
   already exists on entry to this function, although the user should
   ensure that sufficient space is allocated in *Jac* to hold the
   nonzero values to be set; if the existing space is insufficient the
   user may reallocate the data and index arrays as needed.  The
   amount of allocated space in a SUNMATRIX_SPARSE object may be
   accessed using the macro ``SM_NNZ_S`` or the routine
   :c:func:`SUNSparseMatrix_NNZ()`.  The SUNMATRIX_SPARSE type is
   further documented in the section :ref:`SUNMatrix_Sparse`.




.. _ARKStep_CInterface.JTimesFn:

Jacobian information (matrix-vector product)
--------------------------------------------------------------

If the ARKSPILS solver interface is selected
(i.e. :c:func:`ARKSpilsSetLinearSolver()` is called in the
section :ref:`ARKStep_CInterface.Skeleton`), the user may provide a function
of type :c:type:`ARKSpilsJacTimesVecFn` in the following form, to
compute matrix-vector products :math:`Jv`. If such a function is not
supplied, the default is a difference quotient approximation to these
products.


.. c:type:: typedef int (*ARKSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void* user_data, N_Vector tmp)

   This function computes the product :math:`Jv =
   \left(\frac{\partial f_I}{\partial y}\right)v` (or an approximation to it).

   **Arguments:**
      * *v* -- the vector to multiply.
      * *Jv* -- the output vector computed.
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *fy* -- the current value of the vector :math:`f_I(t,y)`.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.
      * *tmp* -- pointer to memory allocated to a variable of type
        ``N_Vector`` which can be used as temporary storage or work space.

   **Return value:**
   The value to be returned by the Jacobian-vector product
   function should be 0 if successful. Any other return value will
   result in an unrecoverable error of the generic Krylov solver,
   in which case the integration is halted.

   **Notes:** If the user's :c:type:`ARKSpilsJacTimesVecFn` function
   uses difference quotient approximations, it may need to access
   quantities not in the argument list.  These include the current
   step size, the error weights, etc.  To obtain these, the
   user will need to add a pointer to the ``ark_mem`` structure to
   their ``user_data``, and then use the ARKStepGet* functions listed
   in :ref:`ARKStep_CInterface.OptionalOutputs`. The unit roundoff can be
   accessed as ``UNIT_ROUNDOFF``, which is defined in the header
   file ``sundials_types.h``.




.. _ARKStep_CInterface.JTSetupFn:

Jacobian information (matrix-vector setup)
--------------------------------------------------------------

If the user's Jacobian-times-vector requires that any Jacobian-related data
be preprocessed or evaluated, then this needs to be done in a
user-supplied function of type :c:type:`ARKSpilsJacTimesSetupFn`,
defined as follows:


.. c:type:: typedef int (*ARKSpilsJacTimesSetupFn)(realtype t, N_Vector y, N_Vector fy, void* user_data)

   This function preprocesses and/or evaluates any Jacobian-related
   data needed by the Jacobian-times-vector routine.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *fy* -- the current value of the vector :math:`f_I(t,y)`.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   The value to be returned by the Jacobian-vector setup
   function should be 0 if successful, positive for a recoverable
   error (in which case the step will be retried), or negative for an
   unrecoverable error (in which case the integration is halted).

   **Notes:**    Each call to the Jacobian-vector setup function is
   preceded by a call to the implicit :c:type:`ARKRhsFn` user
   function with the same :math:`(t,y)` arguments.  Thus, the setup
   function can use any auxiliary data that is computed and saved
   during the evaluation of the implicit ODE right-hand side.

   If the user's :c:type:`ARKSpilsJacTimesSetupFn` function uses
   difference quotient approximations, it may need to access
   quantities not in the argument list.  These include the current
   step size, the error weights, etc.  To obtain these, the
   user will need to add a pointer to the ``ark_mem`` structure to
   their ``user_data``, and then use the ARKStepGet* functions listed
   in :ref:`ARKStep_CInterface.OptionalOutputs`. The unit roundoff can be
   accessed as ``UNIT_ROUNDOFF``, which is defined in the header
   file ``sundials_types.h``.




.. _ARKStep_CInterface.PrecSolveFn:

Preconditioning (linear system solution)
--------------------------------------------------------------

If preconditioning is used with the ARKSPILS solver interface, then
the user must provide a function of type
:c:type:`ARKSpilsPrecSolveFn` to solve the linear system
:math:`Pz=r`, where :math:`P` may be either a left or right
preconditioning matrix.  Here :math:`P` should approximate (at least
crudely) the Newton matrix :math:`A=M-\gamma J`, where :math:`M` is
the mass matrix (typically :math:`M=I` unless working in a
finite-element setting) and :math:`J = \frac{\partial f_I}{\partial
y}`  If preconditioning is done on both sides, the product of the two
preconditioner matrices should approximate :math:`A`.



.. c:type:: typedef int (*ARKSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr, void* user_data)

   This function solves the preconditioner system :math:`Pz=r`.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *fy* -- the current value of the vector :math:`f_I(t,y)`.
      * *r* -- the right-hand side vector of the linear system.
      * *z* -- the computed output solution vector.
      * *gamma* -- the scalar :math:`\gamma` appearing in the Newton
        matrix given by :math:`A=M-\gamma J`.
      * *delta* -- an input tolerance to be used if an iterative method
        is employed in the solution.  In that case, the resdual vector
        :math:`Res = r-Pz` of the system should be made to be less than *delta*
        in the weighted :math:`l_2` norm, i.e. :math:`\left(\sum_{i=1}^n
        \left(Res_i * ewt_i\right)^2 \right)^{1/2} < \delta`, where :math:`\delta =`
        `delta`.  To obtain the ``N_Vector`` *ewt*, call
        :c:func:`ARKStepGetErrWeights()`.
      * *lr* -- an input flag indicating whether the preconditioner
        solve is to use the left preconditioner (*lr* = 1) or the right
        preconditioner (*lr* = 2).
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   The value to be returned by the preconditioner solve
   function is a flag indicating whether it was successful. This value
   should be 0 if successful, positive for a recoverable error (in
   which case the step will be retried), or negative for an
   unrecoverable error (in which case the integration is halted).




.. _ARKStep_CInterface.PrecSetupFn:

Preconditioning (Jacobian data)
--------------------------------------------------------------

If the user's preconditioner requires that any data be preprocessed or
evaluated, then these actions need to occur within a user-supplied
function of type :c:type:`ARKSpilsPrecSetupFn`.


.. c:type:: typedef int (*ARKSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy, booleantype jok, booleantype* jcurPtr, realtype gamma, void* user_data)

   This function preprocesses and/or evaluates Jacobian-related
   data needed by the preconditioner.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *y* -- the current value of the dependent variable vector.
      * *fy* -- the current value of the vector :math:`f_I(t,y)`.
      * *jok* -- is an input flag indicating whether the Jacobian-related
        data needs to be updated. The *jok* argument provides for the
        reuse of Jacobian data in the preconditioner solve function. When
        *jok* = ``SUNFALSE``, the Jacobian-related data should be recomputed
        from scratch. When *jok* = ``SUNTRUE`` the Jacobian data, if saved from the
        previous call to this function, can be reused (with the current
        value of *gamma*). A call with *jok* = ``SUNTRUE`` can only occur
        after a call with *jok* = ``SUNFALSE``.
      * *jcurPtr* -- is a pointer to a flag which should be set to
        ``SUNTRUE`` if Jacobian data was recomputed, or set to ``SUNFALSE`` if
        Jacobian data was not recomputed, but saved data was still reused.
      * *gamma* -- the scalar :math:`\gamma` appearing in the Newton
        matrix given by :math:`A=M-\gamma J`.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   The value to be returned by the preconditioner setup
   function is a flag indicating whether it was successful. This value
   should be 0 if successful, positive for a recoverable error (in
   which case the step will be retried), or negative for an
   unrecoverable error (in which case the integration is halted).

   **Notes:**  The operations performed by this function might include
   forming a crude approximate Jacobian, and performing an LU
   factorization of the resulting approximation to :math:`A = M -
   \gamma J`.

   Each call to the preconditioner setup function is preceded by a
   call to the implicit :c:type:`ARKRhsFn` user function with the
   same :math:`(t,y)` arguments.  Thus, the preconditioner setup
   function can use any auxiliary data that is computed and saved
   during the evaluation of the ODE right-hand side.

   This function is not called in advance of every call to the
   preconditioner solve function, but rather is called only as often
   as needed to achieve convergence in the Newton iteration.

   If the user's :c:type:`ARKSpilsPrecSetupFn` function uses
   difference quotient approximations, it may need to access
   quantities not in the call list. These include the current step
   size, the error weights, etc.  To obtain these, the user will need
   to add a pointer to the ``ark_mem`` structure to their
   ``user_data``, and then use the ARKStepGet* functions listed in
   :ref:`ARKStep_CInterface.OptionalOutputs`. The unit roundoff can be
   accessed as ``UNIT_ROUNDOFF``, which is defined in the header
   file ``sundials_types.h``.



.. _ARKStep_CInterface.MassFn:

Mass matrix information (direct method mass matrix)
---------------------------------------------------------------

If the direct mass matrix linear solver interface is
used (i.e., :c:func:`ARKDlsSetMassLinearSolver()` is called in the
section :ref:`ARKStep_CInterface.Skeleton`), the user must provide a function
of type :c:type:`ARKDlsMassFn` to provide the mass matrix
approximation.



.. c:type:: typedef int (*ARKDlsMassFn)(realtype t, SUNMatrix M, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the mass matrix :math:`M` (or an approximation to it).

   **Arguments:**
      * *N* -- the size of the ODE system.
      * *t* -- the current value of the independent variable.
      * *M* -- the output mass matrix.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.
      * *tmp1*, *tmp2*, *tmp3* -- pointers to memory allocated to
        variables of type ``N_Vector`` which can be used by an
        ARKDlsDenseMassFn as temporary storage or work space.

   **Return value:**
   An *ARKDlsMassFn* function should return 0 if successful, or a
   negative value if it failed unrecoverably (in which case the
   integration is halted, :c:func:`ARKStepEvolve()` returns
   *ARK_MASSSETUP_FAIL* and ARKDLS sets *last_flag* to
   *ARKDLS_MASSFUNC_UNRECVR*).

   **Notes:** Information regarding the structure of the specific
   ``SUNMatrix`` structure (e.g.~number of rows, upper/lower
   bandwidth, sparsity type) may be obtained through using the
   implementation-specific ``SUNMatrix`` interface functions
   (see the section :ref:`SUNMatrix` for details).

   ..
      Prior to calling the user-supplied mass matrix function, the mass
      matrix :math:`M(t)` is zeroed out, so only nonzero elements need to
      be loaded into *M*.

   Prior to calling the user-supplied mass matrix function, the mass
   matrix :math:`M` is zeroed out, so only nonzero elements need to
   be loaded into *M*.

   **dense**:

   ..
      A user-supplied dense mass matrix function must load the *N* by *N*
      dense matrix *M* with an approximation to the mass matrix
      :math:`M(t)`. As discussed above in section :ref:`ARKStep_CInterface.JacobianFn`,
      the accessor macros ``SM_ELEMENT_D`` and ``SM_COLUMN_D`` allow the user
      to read and write dense matrix elements without making explicit
      references to the underlying representation of the SUNMATRIX_DENSE
      type. Similarly, the SUNMATRIX_DENSE type and accessor macros
      ``SM_ELEMENT_D`` and ``SM_COLUMN_D`` are documented in the section
      :ref:`SUNMatrix_Dense`.

   A user-supplied dense mass matrix function must load the *N* by *N*
   dense matrix *M* with an approximation to the mass matrix
   :math:`M`. As discussed above in section :ref:`ARKStep_CInterface.JacobianFn`,
   the accessor macros ``SM_ELEMENT_D`` and ``SM_COLUMN_D`` allow the user
   to read and write dense matrix elements without making explicit
   references to the underlying representation of the SUNMATRIX_DENSE
   type. Similarly, the SUNMATRIX_DENSE type and accessor macros
   ``SM_ELEMENT_D`` and ``SM_COLUMN_D`` are documented in the section
   :ref:`SUNMatrix_Dense`.

   **band**:

   ..
      A user-supplied banded mass matrix function must load
      the band matrix *M* with the elements of the mass matrix
      :math:`M(t)`. As discussed above in section
      :ref:`ARKStep_CInterface.JacobianFn`, the accessor macros ``SM_ELEMENT_B``,
      ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` allow the user to read
      and write band matrix elements without making specific references
      to the underlying representation of the SUNMATRIX_BAND type.
      Similarly, the SUNMATRIX_BAND type and the accessor macros ``SM_ELEMENT_B``,
      ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` are documented in the section
      :ref:`SUNMatrix_Band`.

   A user-supplied banded mass matrix function must load
   the band matrix *M* with the elements of the mass matrix
   :math:`M`. As discussed above in section
   :ref:`ARKStep_CInterface.JacobianFn`, the accessor macros ``SM_ELEMENT_B``,
   ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` allow the user to read
   and write band matrix elements without making specific references
   to the underlying representation of the SUNMATRIX_BAND type.
   Similarly, the SUNMATRIX_BAND type and the accessor macros ``SM_ELEMENT_B``,
   ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` are documented in the section
   :ref:`SUNMatrix_Band`.

   **sparse**:

   ..
      A user-supplied sparse mass matrix function must load the
      compressed-sparse-column (CSR) or compressed-sparse-row (CSR)
      matrix *M* with an approximation to the mass matrix :math:`M(t)`.
      Storage for *M* already exists on entry to this function, although
      the user should ensure that sufficient space is allocated in *M*
      to hold the nonzero values to be set; if the existing space is
      insufficient the user may reallocate the data and row index arrays
      as needed.  The type of *M* is SUNMATRIX_SPARSE, and the amount of
      allocated space in a SUNMATRIX_SPARSE object may be
      accessed using the macro ``SM_NNZ_S`` or the routine
      :c:func:`SUNSparseMatrix_NNZ()`.  The SUNMATRIX_SPARSE type is
      further documented in the section :ref:`SUNMatrix_Sparse`.

   A user-supplied sparse mass matrix function must load the
   compressed-sparse-column (CSR) or compressed-sparse-row (CSR)
   matrix *M* with an approximation to the mass matrix :math:`M`.
   Storage for *M* already exists on entry to this function, although
   the user should ensure that sufficient space is allocated in *M*
   to hold the nonzero values to be set; if the existing space is
   insufficient the user may reallocate the data and row index arrays
   as needed.  The type of *M* is SUNMATRIX_SPARSE, and the amount of
   allocated space in a SUNMATRIX_SPARSE object may be
   accessed using the macro ``SM_NNZ_S`` or the routine
   :c:func:`SUNSparseMatrix_NNZ()`.  The SUNMATRIX_SPARSE type is
   further documented in the section :ref:`SUNMatrix_Sparse`.



.. _ARKStep_CInterface.MTimesFn:

Mass matrix information (matrix-vector product)
--------------------------------------------------------------

If the ARKSPILS solver interface is selected
(i.e. :c:func:`ARKSpilsSetMassLinearSolver()` is called in the
section :ref:`ARKStep_CInterface.Skeleton`), the user must provide a function
of type :c:type:`ARKSpilsMassTimesVecFn` in the following form, to
compute matrix-vector products :math:`Mv`.



.. c:type:: typedef int (*ARKSpilsMassTimesVecFn)(N_Vector v, N_Vector Mv, realtype t, void* mtimes_data)

   This function computes the product :math:`M*v` (or an approximation to it).

   **Arguments:**
      * *v* -- the vector to multiply.
      * *Mv* -- the output vector computed.
      * *t* -- the current value of the independent variable.
      * *mtimes_data* -- a pointer to user data, the same as the
        *mtimes_data* parameter that was passed to :c:func:`ARKSpilsSetMassTimes()`.

   **Return value:**
   The value to be returned by the mass-matrix-vector product
   function should be 0 if successful. Any other return value will
   result in an unrecoverable error of the generic Krylov solver,
   in which case the integration is halted.



.. _ARKStep_CInterface.MTSetupFn:

Mass matrix information (matrix-vector setup)
--------------------------------------------------------------

If the user's mass-matrix-times-vector requires that any mass
matrix-related data be preprocessed or evaluated, then this needs to
be done in a user-supplied function of type
:c:type:`ARKSpilsMassTimesSetupFn`, defined as follows:



.. c:type:: typedef int (*ARKSpilsMassTimesSetupFn)(realtype t, void* mtimes_data)

   This function preprocesses and/or evaluates any mass-matrix-related
   data needed by the mass-matrix-times-vector routine.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *mtimes_data* -- a pointer to user data, the same as the
        *mtimes_data* parameter that was passed to :c:func:`ARKSpilsSetMassTimes()`.

   **Return value:**
   The value to be returned by the mass-matrix-vector setup
   function should be 0 if successful. Any other return value will
   result in an unrecoverable error of the ARKSPILS mass matrix solver
   interface, in which case the integration is halted.



.. _ARKStep_CInterface.MassPrecSolveFn:

Mass matrix preconditioning (linear system solution)
--------------------------------------------------------------

If preconditioning is used with the ARKSPILS mass matrix solver
interface, then the user must provide a function of type
:c:type:`ARKSpilsMassPrecSolveFn` to solve the linear system
:math:`Pz=r`, where :math:`P` may be either a left or right
preconditioning matrix.  Here :math:`P` should approximate (at least
crudely) the mass matrix :math:`M`.  If preconditioning is done on
both sides, the product of the two preconditioner matrices should
approximate :math:`M`.


.. c:type:: typedef int (*ARKSpilsMassPrecSolveFn)(realtype t, N_Vector r, N_Vector z, realtype delta, int lr, void* user_data)

   This function solves the preconditioner system :math:`Pz=r`.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *r* -- the right-hand side vector of the linear system.
      * *z* -- the computed output solution vector.
      * *delta* -- an input tolerance to be used if an iterative method
        is employed in the solution.  In that case, the resdual vector
        :math:`Res = r-Pz` of the system should be made to be less than *delta*
        in the weighted :math:`l_2` norm, i.e. :math:`\left(\sum_{i=1}^n
        \left(Res_i * ewt_i\right)^2 \right)^{1/2} < \delta`, where :math:`\delta =`
        *delta*.  To obtain the ``N_Vector`` *ewt*, call
        :c:func:`ARKStepGetErrWeights()`.
      * *lr* -- an input flag indicating whether the preconditioner
        solve is to use the left preconditioner (*lr* = 1) or the right
        preconditioner (*lr* = 2).
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   The value to be returned by the preconditioner solve
   function is a flag indicating whether it was successful. This value
   should be 0 if successful, positive for a recoverable error (in
   which case the step will be retried), or negative for an
   unrecoverable error (in which case the integration is halted).




.. _ARKStep_CInterface.MasPrecSetupFn:

Mass matrix preconditioning (mass matrix data)
--------------------------------------------------------------

If the user's mass matrix preconditioner requires that any problem
data be preprocessed or evaluated, then these actions need to occur
within a user-supplied function of type
:c:type:`ARKSpilsMassPrecSetupFn`.



.. c:type:: typedef int (*ARKSpilsMassPrecSetupFn)(realtype t, void* user_data)

   This function preprocesses and/or evaluates mass-matrix-related
   data needed by the preconditioner.

   **Arguments:**
      * *t* -- the current value of the independent variable.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter that was passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
   The value to be returned by the mass matrix preconditioner setup
   function is a flag indicating whether it was successful. This value
   should be 0 if successful, positive for a recoverable error (in
   which case the step will be retried), or negative for an
   unrecoverable error (in which case the integration is halted).

   **Notes:**  The operations performed by this function might include
   forming a mass matrix and performing an incomplete
   factorization of the result.  Although such operations would
   typically be performed only once at the beginning of a simulation,
   these may be required if the mass matrix can change as a function
   of time.

   If both this function and a :c:type:`ARKSpilsMassTimesSetupFn` are
   supplied, all calls to this function will be preceded by a call to
   the :c:type:`ARKSpilsMassTimesSetupFn`, so any setup performed
   there may be reused.


.. _ARKStep_CInterface.VecResizeFn:

Vector resize function
--------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatial adaptivity in a
PDE simulation), the ARKStep integrator may be "resized" between
integration steps, through calls to the :c:func:`ARKStepResize()`
function. Typically, when performing adaptive simulations the solution
is stored in a customized user-supplied data structure, to enable
adaptivity without repeated allocation/deallocation of memory.  In
these scenarios, it is recommended that the user supply a customized
vector kernel to interface between SUNDIALS and their problem-specific
data structure.  If this vector kernel includes a function of type
:c:type:`ARKVecResizeFn` to resize a given vector implementation, then
this function may be supplied to :c:func:`ARKStepResize()` so that all
internal ARKStep vectors may be resized, instead of deleting and
re-creating them at each call.  This resize function should have the
following form:


.. c:type:: typedef int (*ARKVecResizeFn)(N_Vector y, N_Vector ytemplate, void* user_data)

   This function resizes the vector *y* to match the dimensions of the
   supplied vector, *ytemplate*.

   **Arguments:**
      * *y* -- the vector to resize.
      * *ytemplate* -- a vector of the desired size.
      * *user_data* -- a pointer to user data, the same as the
        *resize_data* parameter that was passed to :c:func:`ARKStepResize()`.

   **Return value:**
   An *ARKVecResizeFn* function should return 0 if it successfully
   resizes the vector *y*, and a non-zero value otherwise.

   **Notes:**  If this function is not supplied, then ARKStep will
   instead destroy the vector *y* and clone a new vector *y* off of
   *ytemplate*.
