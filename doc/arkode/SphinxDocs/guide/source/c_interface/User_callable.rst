..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _CInterface.UserCallable:

User-callable functions
=============================

This section describes the ARKode functions that are called by the
user to setup and then solve an IVP. Some of these are
required. However, starting with the section
:ref:`CInterface.OptionalInputs`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be
skipped for a casual use of ARKode. In any 
case, refer to the preceding section, :ref:`CInterface.Skeleton`, for
the correct order of these calls. 

On an error, each user-callable function returns a negative value and
sends an error message to the error handler routine, which prints the
message on ``stderr`` by default. However, the user can set a file as
error output or can provide her own error handler function
(see the section :ref:`CInterface.OptionalInputs` for details).



.. _CInterface.Initialization:

ARKode initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* ARKodeCreate()

   This function creates an internal memory block for a problem to be
   solved by ARKode. 

   **Arguments:**  None

   **Return value:**  If successful, a pointer to initialized problem memory
   of type ``void*``, to be passed to all user-facing ARKode routines
   listed below.  If unsuccessful, a ``NULL`` pointer will be
   returned, and an error message will be printed to ``stderr``.



.. c:function:: int ARKStepCreate(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, realtype t0, N_Vector y0)

   This function allocates and initializes memory for a problem to
   be solved using the ARK timestepping module in ARKode. 

   ..
      **Arguments:**
         * *arkode_mem* -- pointer to the ARKode memory block
           (that was returned by :c:func:`ARKodeCreate()`)
         * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
           defining the explicit portion of the right-hand side function in 
           :math:`M(t)\, \dot{y} = f_E(t,y) + f_I(t,y)` 
         * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
           defining the implicit portion of the right-hand side function in 
           :math:`M(t)\, \dot{y} = f_E(t,y) + f_I(t,y)`
         * *t0* -- the initial value of :math:`t`
         * *y0* -- the initial condition vector :math:`y(t_0)`

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block
        (that was returned by :c:func:`ARKodeCreate()`)
      * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit portion of the right-hand side function in 
        :math:`M\, \dot{y} = f_E(t,y) + f_I(t,y)` 
      * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit portion of the right-hand side function in 
        :math:`M\, \dot{y} = f_E(t,y) + f_I(t,y)`
      * *t0* -- the initial value of :math:`t`
      * *y0* -- the initial condition vector :math:`y(t_0)`

   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:**  *One and only one* of :c:func:`ARKStepCreate()` and
   :c:func:`ERKStepCreate()` should be called to define and initialize
   memory for solving the problem.
    

.. c:function:: int ERKStepCreate(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0)

   This function allocates and initializes memory for a problem to
   be solved using the ERK timestepping module in ARKode. 

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block
        (that was returned by :c:func:`ARKodeCreate()`)
      * *f* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function in 
        :math:`\dot{y} = f(t,y)` 
      * *t0* -- the initial value of :math:`t`
      * *y0* -- the initial condition vector :math:`y(t_0)`

   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:**  *One and only one* of :c:func:`ARKStepCreate()` and
   :c:func:`ERKStepCreate()` should be called to define and initialize
   memory for solving the problem.

   
.. c:function:: void ARKodeFree(void* arkode_mem)

   This function frees the problem memory *arkode_mem* created by
   :c:func:`ARKodeCreate()`, :c:func:`ARKStepCreate()` and
   :c:func:`ERKStepCreate()`.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:**  None



.. _CInterface.Tolerances:

ARKode tolerance specification functions
------------------------------------------------------

These functions specify the integration tolerances. One of them
**should** be called before the first call to :c:func:`ARKode()`; otherwise
default values of ``reltol = 1e-4`` and ``abstol = 1e-9`` will be
used, which may be entirely incorrect for a specific problem.

The integration tolerances ``reltol`` and ``abstol`` define a vector
of error weights, ``ewt``.  In the case of
:c:func:`ARKodeSStolerances()`, this vector has components 

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol);

whereas in the case of :c:func:`ARKodeSVtolerances()` the vector components
are given by 

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol[i]);

This vector is used in all error and convergence tests, which use a
weighted RMS norm on all error-like vectors v:

.. math::
    \|v\|_{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N (v_i\; ewt_i)^2 \right)^{1/2},

where :math:`N` is the problem dimension.

Alternatively, the user may supply a custom function to supply the
``ewt`` vector, through a call to :c:func:`ARKodeWFtolerances()`.



.. c:function:: int ARKodeSStolerances(void* arkode_mem, realtype reltol, realtype abstol)

   This function specifies scalar relative and absolute tolerances.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *reltol* -- scalar relative tolerance
      * *abstol* -- scalar absolute tolerance
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ARKode memory was not allocated by the timestepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int ARKodeSVtolerances(void* arkode_mem, realtype reltol, N_Vector abstol)

   This function specifies a scalar relative tolerance and a vector
   absolute tolerance (a potentially different absolute tolerance for
   each vector component).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *reltol* -- scalar relative tolerance
      * *abstol* -- vector containing the absolute tolerances for each
        solution component
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ARKode memory was not allocated by the timestepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int ARKodeWFtolerances(void* arkode_mem, ARKEwtFn efun)

   This function specifies a user-supplied function *efun* to compute
   the error weight vector ``ewt``. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *efun* -- the name of the function (of type :c:func:`ARKEwtFn()`)
        that implements the error weight vector computation.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ARKode memory was not allocated by the timestepping module



..
   Moreover, for problems involving a non-identity mass matrix
   :math:`M(t) \ne I`, the units of the solution vector :math:`y` may differ
   from the units of the IVP, posed for the vector :math:`M(t)\,y`.  When this
Moreover, for problems involving a non-identity mass matrix
:math:`M \ne I`, the units of the solution vector :math:`y` may differ
from the units of the IVP, posed for the vector :math:`My`.  When this
occurs, iterative solvers for the Newton linear systems and the mass
matrix linear systems may require a different set of tolerances.
Since the relative tolerance is dimensionless, but the absolute
tolerance encodes a measure of what is "small" in the units of the
respective quantity, a user may optionally define absolute tolerances
in the equation units.  In this case, ARKode defines a vector of residual
weights, ``rwt`` for measuring convergence of these iterative solvers.
In the case of :c:func:`ARKodeResStolerance()`, this vector has components 

.. code-block:: c

   rwt[i] = 1.0/(reltol*abs(My[i]) + rabstol);

whereas in the case of :c:func:`ARKodeResVtolerance()` the vector components
are given by 

.. code-block:: c

   rwt[i] = 1.0/(reltol*abs(My[i]) + rabstol[i]);

This residual weight vector is used in all iterative solver
convergence tests, which similarly use a weighted RMS norm on all
residual-like vectors v: 

.. math::
    \|v\|_{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N (v_i\; rwt_i)^2 \right)^{1/2},

where :math:`N` is the problem dimension.

As with the error weight vector, the user may supply a custom function
to supply the ``rwt`` vector, through a call to
:c:func:`ARKodeResFtolerance()`.  Further information on all three of
these functions is provided below.



.. c:function:: int ARKodeResStolerance(void* arkode_mem, realtype abstol)

   This function specifies a scalar absolute residual tolerance.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *rabstol* -- scalar absolute residual tolerance
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ARKode memory was not allocated by the timestepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int ARKodeResVtolerance(void* arkode_mem, N_Vector rabstol)

   This function specifies a vector of absolute residual tolerances.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *rabstol* -- vector containing the absolute residual
	tolerances for each solution component
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ARKode memory was not allocated by the timestepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int ARKodeResFtolerance(void* arkode_mem, ARKRwtFn rfun)

   This function specifies a user-supplied function *rfun* to compute
   the residual weight vector ``rwt``. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *rfun* -- the name of the function (of type :c:func:`ARKRwtFn()`)
        that implements the residual weight vector computation.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ARKode memory was not allocated by the timestepping module



General advice on the choice of tolerances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many users, the appropriate choices for tolerance values in
``reltol``, ``abstol`` and ``rabstol`` are a concern. The following pieces
of advice are relevant. 

(1) The scalar relative tolerance ``reltol`` is to be set to control
    relative errors. So a value of :math:`10^{-4}` means that errors
    are controlled to .01%. We do not recommend using ``reltol`` larger
    than :math:`10^{-3}`. On the other hand, ``reltol`` should not be so
    small that it is comparable to the unit roundoff of the machine
    arithmetic (generally around :math:`10^{-15}` for double-precision). 

(2) The absolute tolerances ``abstol`` (whether scalar or vector) need
    to be set to control absolute errors when any components of the
    solution vector :math:`y` may be so small that pure relative error
    control is meaningless.  For example, if :math:`y_i` starts at some
    nonzero value, but in time decays to zero, then pure relative
    error control on :math:`y_i` makes no sense (and is overly costly)
    after :math:`y_i` is below some noise level. Then ``abstol`` (if
    scalar) or ``abstol[i]`` (if a vector) needs to be set to that
    noise level. If the different components have different noise
    levels, then ``abstol`` should be a vector.  For example, see the
    example problem ``ark_robertson.c``, and the discussion
    of it in the ARKode Examples Documentation [R2013]_.  In that
    problem, the three components vary betwen 0 and 1, and have
    different noise levels; hence the ``atols`` vector therein. It is
    impossible to give any general advice on ``abstol`` values,
    because the appropriate noise levels are completely 
    problem-dependent. The user or modeler hopefully has some idea as
    to what those noise levels are. 

(3) The residual absolute tolerances ``rabstol`` (whether scalar or
    vector) follow a similar explanation as for ``abstol``, except
    that these should be set to the noise level of the equation
    components, i.e. the noise level of :math:`My`.  For problems in
    which :math:`M=I`, it is recommended that ``rabstol`` be left
    unset, which will default to the already-supplied ``abstol``
    values. 

(4) Finally, it is important to pick all the tolerance values
    conservatively, because they control the error committed on each
    individual step. The final (global) errors are an accumulation of
    those per-step errors, where that accumulation factor is
    problem-dependent.  A general rule of thumb is to reduce the
    tolerances by a factor of 10 from the actual desired limits on
    errors.  I.e. if you want .01% relative accuracy (globally), a good 
    choice for ``reltol`` is :math:`10^{-5}`.  But in any case, it is
    a good idea to do a few experiments with the tolerances to see how
    the computed solution values vary as tolerances are reduced.



Advice on controlling unphysical negative values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many applications, some components in the true solution are always
positive or non-negative, though at times very small.  In the
numerical solution, however, small negative (unphysical) values
can then occur. In most cases, these values are harmless, and simply
need to be controlled, not eliminated, but in other cases any value
that violates a constraint may cause a simulation to halt. For both of
these scenarios the following pieces of advice are relevant. 

(1) The best way to control the size of unwanted negative computed
    values is with tighter absolute tolerances.  Again this requires
    some knowledge of the noise level of these components, which may
    or may not be different for different components. Some
    experimentation may be needed. 

(2) If output plots or tables are being generated, and it is important
    to avoid having negative numbers appear there (for the sake of
    avoiding a long explanation of them, if nothing else), then
    eliminate them, but only in the context of the output medium. Then
    the internal values carried by the solver are unaffected. Remember
    that a small negative value in :math:`y` returned by ARKode, with
    magnitude comparable to ``abstol`` or less, is equivalent to zero
    as far as the computation is concerned. 

(3) The user's right-hand side routines :math:`f_E` and :math:`f_I`
    should never change a negative value in the solution vector :math:`y`
    to a non-negative value in attempt to "fix" this problem,
    since this can lead to numerical instability.  If the :math:`f_E`
    or :math:`f_I` routines cannot tolerate a zero or negative value
    (e.g. because there is a square root or log), then the offending
    value should be changed to zero or a tiny positive number in a
    temporary variable (not in the input :math:`y` vector) for the
    purposes of computing :math:`f_E(t, y)` or :math:`f_I(t, y)`. 

(4) Positivity and non-negativity constraints on components can be
    enforced by use of the recoverable error return feature in the
    user-supplied right-hand side functions, :math:`f_E` and
    :math:`f_I`. When a recoverable error is encountered, ARKode will
    retry the step with a smaller step size, which typically
    alleviates the problem.  However, because this option involves
    some additional overhead cost, it should only be exercised if the
    use of absolute tolerances to control the computed values is
    unsuccessful.

    

.. _CInterface.LinearSolvers:

Linear solver interface functions
-------------------------------------------

As previously explained, the Newton iterations used in solving
implicit systems within ARKode requires the solution of linear
systems of the form 

.. math::
   {\mathcal A}\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right)

where 

.. math::
   {\mathcal A} \approx M - \gamma J, \qquad J = \frac{\partial f_I}{\partial y}.

There are two ARKode linear solver interfaces currently available for this
task: ARKDLS and ARKSPILS.

The first corresponds to the use of Direct Linear Solvers, and
utilizes ``SUNMatrix`` objects to store the approximate Jacobian
:math:`J`, the Newton matrix :math:`{\mathcal A}`, the mass matrix
:math:`M`, and factorizations used throughout the solution process.

The second corresponds to the use of Scaled, Preconditioned, Iterative
Linear Solvers, utilizing matrix-free Krylov methods to solve the
Newton systems of equations.  With most of these methods,
preconditioning can be done on the left only, on the right only, on
both the left and the right, or not at all.  The exceptions to this
rule are SPFGMR that supports right preconditioning only and PCG
that performs symmetric preconditioning.  For the specification
of a preconditioner, see the iterative linear solver portions of the
sections :ref:`CInterface.OptionalInputs` and
:ref:`CInterface.UserSupplied`. 

If preconditioning is done, user-supplied functions should be used to
define left and right preconditioner matrices :math:`P_1` and 
:math:`P_2` (either of which could be the identity matrix), such that
the product :math:`P_{1}P_{2}` approximates the Newton matrix
:math:`{\mathcal A} = M - \gamma J`.  

To specify a generic linear solver for ARKode to use for the Newton
systems, after the call to :c:func:`ARKodeCreate()` but before any
calls to :c:func:`ARKode()`, the user's program must create the
appropriate ``SUNLinearSolver`` object and call either of the
functions :c:func:`ARKDlsSetLinearSolver()` or
:c:func:`ARKSpilsSetLinearSolver()`, as documented below.  The first
argument passed to these functions is the ARKode memory pointer 
returned by :c:func:`ARKodeCreate()`; the second argument passed to
these functions is the desired ``SUNLinearSolver`` object to use
for solving Newton systems.  A call to one of these functions initializes the
appropriate ARKode linear solver interface, linking this to the main
ARKode integrator, and allows the user to specify parameters which are
specific to a particular solver interface. 

The use of each of the generic linear solvers involves certain
constants and possibly some macros, that are likely to be needed in
the user code.  These are available in the corresponding header file
associated with the specific ``SUNMatrix`` or
``SUNLinearSolver`` module in question, as described in the
sections :ref:`SUNMatrix` and :ref:`SUNLinSol`.


.. c:function:: int ARKDlsSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix J)

   This function specifies the direct ``SUNLinearSolver`` object
   that ARKode should use, as well as a template Jacobian
   ``SUNMatrix`` object.  Its use requires inclusion of the
   header file  ``arkode/arkode_direct.h``.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *LS* -- the ``SUNLinearSolver`` object to use.
      * *J* -- the template Jacobian ``SUNMatrix`` object to use.
   
   **Return value:** 
      * *ARKDLS_SUCCESS*   if successful
      * *ARKDLS_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARKDLS_MEM_FAIL*  if there was a memory allocation failure
      * *ARKDLS_ILL_INPUT* if ARKDLS is incompatible with the
        provided *LS* or *J* input objects, or the current
        ``N_Vector`` module.
   
   **Notes:**  The template Jacobian matrix *J* will be used in the
   solve process, so if additional storage is required within the
   ``SUNMatrix`` object (e.g. for factorization of a banded
   matrix), ensure that the input object is allocated with sufficient
   size.

   The ARKDLS linear solver interface is not compatible
   with all implementations of the ``SUNLinearSolver`` and
   ``N_Vector`` modules.  Specifically, ARKDLS requires use of a
   *direct* ``SUNLinearSolver`` object and a serial or threaded
   ``N_Vector`` module.  Additional compatibility limitations
   for each ``SUNLinearSolver`` object (i.e. ``SUNMatrix``
   and ``N_Vector`` object compatibility) are described in the
   section :ref:`SUNLinSol`.  


.. c:function:: int ARKSpilsSetLinearSolver(void* arkode_mem, SUNLinearSolver LS)

   This function specifies the iterative ``SUNLinearSolver`` object
   that ARKode should use, initializing the ARKSPILS scaled,
   preconditioned, iterative linear solver interface.  Its use
   requires inclusion of the header file  ``arkode/arkode_spils.h``.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *LS* -- the ``SUNLinearSolver`` object to use.
   
   **Return value:** 
      * *ARKSPILS_SUCCESS*   if successful
      * *ARKSPILS_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARKSPILS_MEM_FAIL*  if there was a memory allocation failure
      * *ARKSPILS_ILL_INPUT* if ARKSPILS is incompatible with the
        provided *LS* input objects, or the current ``N_Vector``
        module. 
   
   **Notes:**  The ARKSPILS linear solver interface is not compatible
   with all implementations of the ``SUNLinearSolver`` and
   ``N_Vector`` modules.  Specifically, ARKSPILS requires use of an
   *iterative* ``SUNLinearSolver`` object, and a minimum
   required set of vector operations must be provided by the current
   ``N_Vector`` module.  Additional compatibility limitations 
   for each ``SUNLinearSolver`` object (i.e. required
   ``N_Vector`` routines) are described in the
   section :ref:`SUNLinSol`.  






.. _CInterface.MassMatrixSolvers:

Mass matrix solver specification functions
-------------------------------------------

As discussed in section :ref:`Mathematics.MassSolve`, if the ODE
system involves a non-identity mass matrix :math:`M\ne I`, then ARKode
must solve linear systems of the form 

.. math::
    M x = b.

The same solver interfaces listed above in the section
:ref:`CInterface.LinearSolvers` may be used for this purpose: ARKDLS
and ARKSPILS.  With the ARKSPILS interface preconditioning can be
applied.  For the specification of a preconditioner, see the iterative
linear solver portions of the sections
:ref:`CInterface.OptionalInputs` and :ref:`CInterface.UserSupplied`. 
If preconditioning is to be performed, user-supplied functions should
be used to define left and right preconditioner matrices :math:`P_1` and 
:math:`P_2` (either of which could be the identity matrix), such that
the product :math:`P_{1}P_{2}` approximates the mass matrix :math:`M`.  

To specify a generic linear solver for ARKode to use for mass matrix
systems, after the call to :c:func:`ARKodeCreate()` but before any 
calls to :c:func:`ARKode()`, the user's program must create the
appropriate ``SUNLinearSolver`` object and call either of the
functions :c:func:`ARKDlsSetMassLinearSolver()` or
:c:func:`ARKSpilsSetMassLinearSolver()`, as documented below.  The
first argument passed to these functions is the ARKode memory pointer
returned by :c:func:`ARKodeCreate()`; the second argument passed to
these functions is the desired ``SUNLinearSolver`` object to use
for solving mass matrix systems.  A call to one of these functions
initializes the appropriate ARKode mass matrix linear solver
interface, linking this to the main ARKode integrator, and allows the
user to specify parameters which are specific to a particular solver
interface.  

The use of each of the generic linear solvers involves certain
constants and possibly some macros, that are likely to be needed in
the user code.  These are available in the corresponding header file
associated with the specific ``SUNMatrix`` or
``SUNLinearSolver`` module in question, as described in the
sections :ref:`SUNMatrix` and :ref:`SUNLinSol`.


Note: if the user program includes linear solvers for *both* the
Newton and mass matrix systems, these must have the same type:

* If both are *direct*, then they must utilize the same
  ``SUNMatrix`` type.  In this case, both the Newton and mass
  matrix linear solver interfaces can use the same
  ``SUNLinearSolver`` object, although different objects
  (e.g. with different solver parameters) are also allowed.

* If both are *iterative*, then the Newton and mass matrix
  ``SUNLinearSolver`` objects must be different.  These may even
  use different solver algorithms (SPGMR, SPBCGS, etc.), if desired.
  For example, if the mass matrix is symmetric but the Jacobian is not,
  then PCG may be used for the mass matrix systems and SPGMR for the
  Newton systems.
     
As with the Newton system solvers, the mass matrix linear system
solvers listed below are all built on top of generic SUNDIALS solver
modules. 

.. c:function:: int ARKDlsSetMassLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix M, booleantype time_dep)

   This function specifies the direct ``SUNLinearSolver`` object
   that ARKode should use for mass matrix systems, as well as a
   template ``SUNMatrix`` object.  Its use requires inclusion of the
   header file  ``arkode/arkode_direct.h``.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *LS* -- the ``SUNLinearSolver`` object to use.
      * *M* -- the template mass ``SUNMatrix`` object to use.
      * *time_dep* -- flag denoting whether the mass matrix depends on
        the independent variable (:math:`M = M(t)`) or not (:math:`M
        \ne M(t)`).  Use ``SUNTRUE`` to indicate time-dependence of the
        mass matrix.
        *Currently, only values of "SUNFALSE" are supported*
   
   **Return value:** 
      * *ARKDLS_SUCCESS*   if successful
      * *ARKDLS_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARKDLS_MEM_FAIL*  if there was a memory allocation failure
      * *ARKDLS_ILL_INPUT* if ARKDLS is incompatible with the
        provided *LS* or *M* input objects, or the current
        ``N_Vector`` module.
   
   **Notes:**  The template mass matrix *M* will be used in the
   solve process, so if additional storage is required within the
   ``SUNMatrix`` object (e.g. for factorization of a banded
   matrix), ensure that the input object is allocated with sufficient
   size.

   If called with *time_dep* set to ``SUNFALSE``, then the mass matrix is
   only computed and factored once, with the results reused throughout
   the entire ARKode simulation.

   Unlike the system Jacobian, the system mass matrix cannot be
   approximated using finite-differences of any functions provided to
   ARKode.  Hence, use of the ARKDLS mass matrix solver interface
   requires the user to provide a mass-matrix constructor routine
   (see :c:type:`ARKDlsMassFn` and :c:func:`ARKDlsSetMassFn()`).
   
   The ARKDLS linear solver interface is not compatible
   with all implementations of the ``SUNLinearSolver`` and
   ``N_Vector`` modules.  Specifically, ARKDLS requires use of a
   *direct* ``SUNLinearSolver`` object and a serial or threaded
   ``N_Vector`` module.  Additional compatibility limitations
   for each ``SUNLinearSolver`` object (i.e. ``SUNMatrix``
   and ``N_Vector`` object compatibility) are described in the
   section :ref:`SUNLinSol`.

   

.. c:function:: int ARKSpilsSetMassLinearSolver(void* arkode_mem, SUNLinearSolver LS, booleantype time_dep)

   This function specifies the iterative ``SUNLinearSolver`` object
   that ARKode should use for mass matrix systems, initializing the
   ARKSPILS scaled, preconditioned, iterative mass matrix linear
   solver interface.  Its use requires inclusion of the header file
   ``arkode/arkode_spils.h``. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *LS* -- the ``SUNLinearSolver`` object to use.
      * *time_dep* -- flag denoting whether the mass matrix depends on
        the independent variable (:math:`M = M(t)`) or not (:math:`M
        \ne M(t)`).  Use ``SUNTRUE`` to indicate time-dependence of the
        mass matrix.
        *Currently, only values of "SUNFALSE" are supported*
   
   **Return value:** 
      * *ARKSPILS_SUCCESS*   if successful
      * *ARKSPILS_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARKSPILS_MEM_FAIL*  if there was a memory allocation failure
      * *ARKSPILS_ILL_INPUT* if ARKSPILS is incompatible with the
        provided *LS* input objects, or the current ``N_Vector``
        module. 
   
   **Notes:**  If called with *time_dep* set to ``SUNFALSE``, then the
   mass matrix-vector-product (if supplied) is only set up once, and
   the mass matrix preconditioner (if supplied) is only set up once,
   with the results reused throughout the entire ARKode simulation.

   Unlike the system Jacobian, the system mass matrix-vector-product
   cannot be approximated using finite-differences of any functions
   provided to ARKode.  Hence, use of the ARKSPILS mass matrix solver
   interface requires the user to provide a mass-matrix-times-vector
   product routine (see :c:type:`ARKSpilsMassTimesVecFn` and
   :c:func:`ARKSpilsSetMassTimes()`). 
   
   The ARKSPILS linear solver interface is not compatible
   with all implementations of the ``SUNLinearSolver`` and
   ``N_Vector`` modules.  Specifically, ARKSPILS requires use of an
   *iterative* ``SUNLinearSolver`` object, and a minimum
   required set of vector operations must be provided by the current
   ``N_Vector`` module.  Additional compatibility limitations 
   for each ``SUNLinearSolver`` object (i.e. required
   ``N_Vector`` routines) are described in the section :ref:`SUNLinSol`.  



.. _CInterface.RootFinding:

Rootfinding initialization function
--------------------------------------

As described in the section :ref:`Mathematics.Rootfinding`, while
solving the IVP ARKode has the capability to find the roots of a set
of user-defined functions.  To activate the root-finding algorithm,
call the following function.  This is normally called only once, prior
to the first call to :c:func:`ARKode()`, but if the rootfinding
problem is to be changed during the solution,
:c:func:`ARKodeRootInit()` can also be called prior to a continuation
call to :c:func:`ARKode()`. 


.. c:function:: int ARKodeRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved during the
   integration of the ODE system.  It must be called after
   :c:func:`ARKodeCreate()`, and before :c:func:`ARKode()`.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nrtfn* -- number of functions :math:`g_i`, an integer :math:`\ge` 0.
      * *g* -- name of user-supplied function, of type :c:func:`ARKRootFn()`,
        defining the functions :math:`g_i` whose roots are sought. 
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_MEM_FAIL*  if there was a memory allocation failure
      * *ARK_ILL_INPUT* if *nrtfn* is greater than zero but *g* = ``NULL``.
   
   **Notes:** To disable the rootfinding feature after it has already
   been initialized, or to free memory associated with ARKode's
   rootfinding module, call *ARKodeRootInit* with *nrtfn = 0*.  

   Similarly, if a new IVP is to be solved with a call to
   :c:func:`ARKodeReInit()`, where the new IVP has no rootfinding
   problem but the prior one did, then call *ARKodeRootInit* with
   *nrtfn = 0*.




.. _CInterface.Integration:

ARKode solver function
-------------------------

This is the central step in the solution process -- the call to perform
the integration of the IVP.  One of the input arguments (*itask*)
specifies one of two modes as to where ARKode is to return a
solution.  These modes are modified if the user has set a stop time
(with a call to the optional input function :c:func:`ARKodeSetStopTime()`) or
has requested rootfinding. 



.. c:function:: int ARKode(void* arkode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *tout* -- the next time at which a computed solution is desired
      * *yout* -- the computed solution vector
      * *tret* -- the time corresponding to *yout* (output)
      * *itask* -- a flag indicating the job of the solver for the next
        user step. 

	The *ARK_NORMAL* option causes the solver to take internal steps
	until it has reached or just passed the user-specified *tout*
	parameter. The solver then interpolates in order to return an
	approximate value of :math:`y(tout)`.  This interpolation may be
        slightly less accurate than the full time step solutions
	produced by the solver, since the interpolation uses a cubic
	Hermite polynomial even when the RK method is of higher order.

	To ensure that this returned value has full method accuracy,
	issue a call to :c:func:`ARKodeSetStopTime()` before the call
	to ARKode to specify a fixed stop time to end the time step
	and return to the user.  Once the integrator returns at a
	`tstop` time, any future testing for *tstop* is disabled (and
	can be reenabled only though a new call to
	:c:func:`ARKodeSetStopTime()`).  

	The *ARK_ONE_STEP* option tells the solver to take just one
	internal step and then return the solution at the point
	reached by that step.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_ROOT_RETURN* if ARKode succeeded, and found one or more roots.
        If *nrtfn* is greater than 1, call :c:func:`ARKodeGetRootInfo()` to see
        which :math:`g_i` were found to have a root at (*\*tret*). 
      * *ARK_TSTOP_RETURN* if ARKode succeeded and returned at *tstop*.
      * *ARK_MEM_NULL* if the *arkode_mem* argument was ``NULL``.
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if one of the inputs to ARKode is illegal, or
        some other input to the solver was either illegal or missing.
	Details will be provided in the error message.  Typical causes
	of this failure:
	
	(a) The tolerances have not been set. 

	(b) A component of the error weight vector became zero during
	    internal time-stepping. 

	(c) The linear solver initialization function (called by the
	    user after calling :c:func:`ARKodeCreate()`) failed to set
	    the linear solver-specific *lsolve* field in
	    *arkode_mem*. 

	(d) A root of one of the root functions was found both at a
	    point :math:`t` and also very near :math:`t`. 

      * *ARK_TOO_MUCH_WORK* if the solver took *mxstep* internal steps
        but could not reach *tout*.  The default value for *mxstep* is
        *MXSTEP_DEFAULT = 500*.
      * *ARK_TOO_MUCH_ACC* if the solver could not satisfy the accuracy
        demanded by the user for some internal step.
      * *ARK_ERR_FAILURE* if error test failures occurred either too many
        times (*ark_maxnef*) during one internal time step or occurred
        with :math:`|h| = h_{min}`. 
      * *ARK_CONV_FAILURE* if either convergence test failures occurred
        too many times (*ark_maxncf*) during one internal time step or
        occurred with :math:`|h| = h_{min}`. 
      * *ARK_LINIT_FAIL* if the linear solver's initialization function failed.
      * *ARK_LSETUP_FAIL* if the linear solver's setup routine failed in
        an unrecoverable manner.
      * *ARK_LSOLVE_FAIL* if the linear solver's solve routine failed in
        an unrecoverable manner.
      * *ARK_MASSINIT_FAIL* if the mass matrix solver's
	initialization function failed. 
      * *ARK_MASSSETUP_FAIL* if the mass matrix solver's setup routine
	failed.
      * *ARK_MASSSOLVE_FAIL* if the mass matrix solver's solve routine
	failed.
   
   **Notes:** The input vector *yout* can use the same memory as the
   vector *y0* of initial conditions that was passed to
   :c:func:`ARKStepCreate()` or :c:func:`ERKStepCreate()`.
   
   In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
   only to get the direction and a rough scale of the independent
   variable. All failure return values are negative and so testing the
   return argument for negative values will trap all ARKode failures.
   
   On any error return in which one or more internal steps were taken
   by ARKode, the returned values of *tret* and *yout* correspond to
   the farthest point reached in the integration.  On all other error
   returns, *tret* and *yout* are left unchanged from those provided
   to the routine.




.. _CInterface.OptionalInputs:

Optional input functions
-------------------------

There are numerous optional input parameters that control the behavior
of the ARKode solver, each of which may be modified from its default
value through calling an appropriate input function.  The following
tables list all optional input functions, grouped by which aspect of
ARKode they control.  Detailed information on the calling syntax and
arguments for each function are then provided following each table.  

The optional inputs are grouped into the following categories:

* General ARKode options (:ref:`CInterface.ARKodeInputTable`), 
* General ARK time-stepper module options (:ref:`CInterface.ARKStepInputTable`), 
* General ERK time-stepper module options (:ref:`CInterface.ERKStepInputTable`), 
* IVP method solver options (:ref:`CInterface.ARKodeMethodInputTable`), 
* Step adaptivity solver options (:ref:`CInterface.ARKodeAdaptivityInputTable`), 
* Implicit stage solver options (:ref:`CInterface.CInterface.ARKodeSolverInputTable`), 
* Direct linear solver interface options (:ref:`CInterface.ARKDlsInputs`),
* Iterative linear solver interface options (:ref:`CInterface.ARKSpilsInputs`).  

For the most casual use of ARKode, relying on the default set of
solver parameters, the reader can skip to the following section,
:ref:`CInterface.UserSupplied`.

Optional input routines may be similarly grouped into three
categories: 

* Options applicable to the main ARKode infrastructure -- these
  function names begin with ``ARKodeSet``

* Options specific to the ARK time-stepping module -- these function
  names begin with ``ARKStepSet``

* Options specific to the ERK time-stepping module -- these function
  names begin with ``ARKodeSet``

These routines are not generally listed separately, as the
aforementioned categories organized around functionality provide a
more intuitive grouping.

We note that, on an error return, all of the optional input functions
send an error message to the error handler function.  We also note
that all error return values are negative, so a test on the return
arguments for negative values will catch all errors. 



.. _CInterface.ARKodeInputTable:

Optional inputs for ARKode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==================================================  ====================================  ==============
Optional input                                      Function name                         Default
==================================================  ====================================  ==============
Return ARKode solver parameters to their defaults   :c:func:`ARKodeSetDefaults()`         internal
Set dense output order                              :c:func:`ARKodeSetDenseOrder()`       3
Supply a pointer to a diagnostics output file       :c:func:`ARKodeSetDiagnostics()`      ``NULL``
Supply a pointer to an error output file            :c:func:`ARKodeSetErrFile()`          ``stderr``
Supply a custom error handler function              :c:func:`ARKodeSetErrHandlerFn()`     internal fn
Disable time step adaptivity (fixed-step mode)      :c:func:`ARKodeSetFixedStep()`        disabled
Supply an initial step size to attempt              :c:func:`ARKodeSetInitStep()`         estimated
Maximum no. of warnings for :math:`t_n+h = t_n`     :c:func:`ARKodeSetMaxHnilWarns()`     10
Maximum no. of internal steps before *tout*         :c:func:`ARKodeSetMaxNumSteps()`      500
Maximum absolute step size                          :c:func:`ARKodeSetMaxStep()`          :math:`\infty`
Minimum absolute step size                          :c:func:`ARKodeSetMinStep()`          0.0
Set a value for :math:`t_{stop}`                    :c:func:`ARKodeSetStopTime()`         :math:`\infty`
Supply a pointer for user data                      :c:func:`ARKodeSetUserData()`         ``NULL``
==================================================  ====================================  ==============



.. c:function:: int ARKodeSetDefaults(void* arkode_mem)

   Resets all optional input parameters to ARKode's original 
   default values.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Does not change the *user_data* pointer or any
   parameters within the specified time-stepping module.
   
   Also leaves alone any data structures or options related to
   root-finding (those can be reset using :c:func:`ARKodeRootInit()`).



.. c:function:: int ARKodeSetDenseOrder(void* arkode_mem, int dord)

   Specifies the order of accuracy for the polynomial interpolant 
   used for dense output (i.e. interpolation of solution output values
   and implicit method predictors). 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *dord* -- requested polynomial order of accuracy
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Allowed values are between 0 and ``min(q,3)``, where ``q`` is
   the order of the overall integration method.



.. c:function:: int ARKodeSetDiagnostics(void* arkode_mem, FILE* diagfp)

   Specifies the file pointer for a diagnostics file where
   all ARKode step adaptivity and solver information is written.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *diagfp* -- pointer to the diagnostics output file
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This parameter can be ``stdout`` or ``stderr``, although the
   suggested approach is to specify a pointer to a unique file opened
   by the user and returned by ``fopen``.  If not called, or if called
   with a ``NULL`` file pointer, all diagnostics output is disabled.
   
   When run in parallel, only one process should set a non-NULL value
   for this pointer, since statistics from all processes would be
   identical.
   


.. c:function:: int ARKodeSetErrFile(void* arkode_mem, FILE* errfp)

   Specifies a pointer to the file where all ARKode warning and error
   messages will be written if the default internal error handling
   function is used.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *errfp* -- pointer to the output file. 
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value for *errfp* is ``stderr``.
    
   Passing a ``NULL`` value disables all future error message output
   (except for the case wherein the ARKode memory pointer is
   ``NULL``.  This use of the function is strongly discouraged.
   
   If used, this routine should be called before any other
   optional input functions, in order to take effect for subsequent
   error messages.



.. c:function:: int ARKodeSetErrHandlerFn(void* arkode_mem, ARKErrHandlerFn ehfun, void* eh_data)

   Specifies the optional user-defined function to be used
   in handling error messages.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *ehfun* -- name of user-supplied error handler function. 
      * *eh_data* -- pointer to user data passed to *ehfun* every time
        it is called
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Error messages indicating that the ARKode solver memory is
   ``NULL`` will always be directed to ``stderr``.




.. c:function:: int ARKodeSetFixedStep(void* arkode_mem, realtype hfixed)

   Disabled time step adaptivity within ARKode, and specifies the
   fixed time step size to use for all internal steps.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hfixed* -- value of the fixed step size to use
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Pass 0.0 to return ARKode to the default (adaptive-step) mode.  

   Use of this function is not recommended, since we it gives no
   assurance of the validity of the computed solutions.  It is
   primarily provided for code-to-code verification testing purposes.

   When using :c:func:`ARKodeSetFixedStep()`, any values provided to
   the functions 
   :c:func:`ARKodeSetInitStep()`,
   :c:func:`ARKStepSetAdaptivityFn()`,
   :c:func:`ARKStepSetMaxErrTestFails()`, 
   :c:func:`ARKStepSetAdaptivityMethod()`,
   :c:func:`ARKStepSetCFLFraction()`, 
   :c:func:`ARKStepSetErrorBias()`, 
   :c:func:`ARKStepSetFixedStepBounds()`,
   :c:func:`ARKStepSetMaxCFailGrowth()`,
   :c:func:`ARKStepSetMaxEFailGrowth()`,
   :c:func:`ARKStepSetMaxFirstGrowth()`,
   :c:func:`ARKStepSetMaxGrowth()`, 
   :c:func:`ARKStepSetSafetyFactor()`, 
   :c:func:`ARKStepSetSmallNumEFails()`,
   :c:func:`ARKStepSetStabilityFn()`,
   :c:func:`ERKStepSetAdaptivityFn()`,
   :c:func:`ERKStepSetMaxErrTestFails()`, 
   :c:func:`ERKStepSetAdaptivityMethod()`,
   :c:func:`ERKStepSetCFLFraction()`, 
   :c:func:`ERKStepSetErrorBias()`, 
   :c:func:`ERKStepSetFixedStepBounds()`,
   :c:func:`ERKStepSetMaxEFailGrowth()`,
   :c:func:`ERKStepSetMaxFirstGrowth()`,
   :c:func:`ERKStepSetMaxGrowth()`, 
   :c:func:`ERKStepSetSafetyFactor()`, 
   :c:func:`ERKStepSetSmallNumEFails()` and
   :c:func:`ERKStepSetStabilityFn()`,
   will be ignored, since temporal adaptivity is disabled. 

   If both :c:func:`ARKodeSetFixedStep()` and
   :c:func:`ARKodeSetStopTime()` are used, then the fixed step size
   will be used for all steps until the final step preceding the
   provided stop time (which may be shorter).  To resume use of the
   previous fixed step size, another call to
   :c:func:`ARKodeSetFixedStep()` must be made prior to calling
   :c:func:`ARKode()` to resume integration.

   It is *not* recommended that :c:func:`ARKodeSetFixedStep()` be used
   in concert with :c:func:`ARKodeSetMaxStep()` or
   :c:func:`ARKodeSetMinStep()`, since at best those latter two
   routines will provide no useful information to the solver, and at
   worst they may interfere with the desired fixed step size.





.. c:function:: int ARKodeSetInitStep(void* arkode_mem, realtype hin)

   Specifies the initial time step size ARKode should use after
   initialization or reinitialization.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hin* -- value of the initial step to be attempted :math:`(\ge 0)`
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Pass 0.0 to use the default value.  

   By default, ARKode estimates the initial step size to be the
   solution :math:`h` of the equation :math:`\left\| \frac{h^2
   \ddot{y}}{2}\right\| = 1`, where :math:`\ddot{y}` is an estimated
   value of the second derivative of the solution at *t0*.  




.. c:function:: int ARKodeSetMaxHnilWarns(void* arkode_mem, int mxhnil)

   Specifies the maximum number of messages issued by the
   solver to warn that :math:`t+h=t` on the next internal step, before
   ARKode will instead return with an error.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *mxhnil* -- maximum allowed number of warning messages (>0).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value is 10; set *mxhnil* to zero to specify
   this default.

   A negative value indicates that no warning messages should be issued.




.. c:function:: int ARKodeSetMaxNumSteps(void* arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time, before ARKode
   will return with an error.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *mxsteps* -- maximum allowed number of internal steps.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Passing *mxsteps* = 0 results in ARKode using the
   default value (500).
   
   Passing *mxsteps* < 0 disables the test (not recommended).



.. c:function:: int ARKodeSetMaxStep(void* arkode_mem, realtype hmax)

   Specifies the upper bound on the magnitude of the time step size.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hmax* -- maximum absolute value of the time step size :math:`(\ge 0)`
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Pass *hmax* :math:`\le 0.0` to set the default value of :math:`\infty`.  



.. c:function:: int ARKodeSetMinStep(void* arkode_mem, realtype hmin)

   Specifies the lower bound on the magnitude of the time step size.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hmin* -- minimum absolute value of the time step size :math:`(\ge 0)`
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Pass *hmin* :math:`\le 0.0` to set the default value of 0.



.. c:function:: int ARKodeSetStopTime(void* arkode_mem, realtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *tstop* -- stopping time for the integrator.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default is that no stop time is imposed.




.. c:function:: int ARKodeSetUserData(void* arkode_mem, void* user_data)

   Specifies the user data block *user_data* and
   attaches it to the main ARKode memory block.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *user_data* -- pointer to the user data.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** If specified, the pointer to *user_data* is passed to all
   user-supplied functions for which it is an argument; otherwise
   ``NULL`` is passed.
   
   If *user_data* is needed in user linear solver or preconditioner
   functions, the call to this function must be made *before* the call
   to specify the linear solver.




.. _CInterface.ARKStepInputTable:

Optional general inputs for ARKStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==================================================  =====================================  ==============
Optional input                                      Function name                          Default
==================================================  =====================================  ==============
Return ARKStep solver parameters to their defaults  :c:func:`ARKStepSetDefaults()`         internal
Maximum no. of ARKStep error test failures          :c:func:`ARKStepSetMaxErrTestFails()`  7
Set 'optimal' adaptivity params for a method        :c:func:`ARKStepSetOptimalParams()`    internal
==================================================  =====================================  ==============



.. c:function:: int ARKStepSetDefaults(void* arkode_mem)

   Resets all optional input parameters to ARKStep's original 
   default values.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Does not change problem-defining function pointers *fe*
   and *fi* or the *user_data* pointer.  

   
.. c:function:: int ARKStepSetMaxErrTestFails(void* arkode_mem, int maxnef)

   Specifies the maximum number of error test failures
   permitted in attempting one step, before returning with an error.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *maxnef* -- maximum allowed number of error test failures :math:`(>0)`
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value is 7; set *maxnef* :math:`\le 0`
   to specify this default.



.. c:function:: int ARKStepSetOptimalParams(void* arkode_mem)

   Sets all adaptivity and solver parameters to our 'best
   guess' values, for a given integration method (ERK, DIRK, ARK) and
   a given method order.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Should only be called after the method order and integration
   method have been set.  These values resulted from repeated testing
   of ARKode's solvers on a variety of training problems.  However,
   all problems are different, so these values may not be optimal for
   all users.




   
.. _CInterface.ERKStepInputTable:

Optional general inputs for ERKStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==================================================  =====================================  ==============
Optional input                                      Function name                          Default
==================================================  =====================================  ==============
Return ERKStep solver parameters to their defaults  :c:func:`ERKStepSetDefaults()`         internal
Maximum no. of ERKStep error test failures          :c:func:`ERKStepSetMaxErrTestFails()`  7
==================================================  =====================================  ==============



.. c:function:: int ERKStepSetDefaults(void* arkode_mem)

   Resets all optional input parameters to ERKStep's original 
   default values.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Does not change problem-defining function pointer *f*
   or the *user_data* pointer.  

   
.. c:function:: int ERKStepSetMaxErrTestFails(void* arkode_mem, int maxnef)

   Specifies the maximum number of error test failures
   permitted in attempting one step, before returning with an error.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *maxnef* -- maximum allowed number of error test failures :math:`(>0)`
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value is 7; set *maxnef* :math:`\le 0`
   to specify this default.







.. _CInterface.ARKodeMethodInputTable:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=================================  =================================  ==============
Optional input                     Function name                      Default
=================================  =================================  ==============
Set integrator method order        :c:func:`ARKStepSetOrder()`        4
                                   :c:func:`ERKStepSetOrder()`        4
Specify implicit/explicit problem  :c:func:`ARKStepSetImEx()`         ``SUNTRUE``
Specify explicit problem           :c:func:`ARKStepSetExplicit()`     ``SUNFALSE``
Specify implicit problem           :c:func:`ARKStepSetImplicit()`     ``SUNFALSE``
Set additive RK tables             :c:func:`ARKStepSetARKTables()`    internal
Set explicit RK table              :c:func:`ARKStepSetERKTable()`     internal
                                   :c:func:`ERKStepSetERKTable()`     internal
Set implicit RK table              :c:func:`ARKStepSetIRKTable()`     internal
Specify additive RK table numbers  :c:func:`ARKStepSetARKTableNum()`  internal
Specify explicit RK table number   :c:func:`ARKStepSetERKTableNum()`  internal
                                   :c:func:`ERKStepSetERKTableNum()`  internal
Specify implicit RK table number   :c:func:`ARKStepSetIRKTableNum()`  internal
=================================  =================================  ==============



.. c:function:: int ARKStepSetOrder(void* arkode_mem, int ord)

   Specifies the order of accuracy for the ARK/DIRK/ERK integration
   method in the ARK timestepper module.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *ord* -- requested order of accuracy.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** For explicit methods, the allowed values are :math:`2 \le`
   *ord* :math:`\le 8`.  For implicit methods, the allowed values are
   :math:`2\le` *ord* :math:`\le 5`, and for ImEx methods the allowed
   values are :math:`3 \le` *ord* :math:`\le 5`.  Any illegal input
   will result in the default value of 4. 
   
   Since *ord* affects the memory requirements for the internal
   ARKode memory block, it cannot be changed after the first call to
   :c:func:`ARKode()`, unless :c:func:`ARKStepReInit()` is called.



.. c:function:: int ERKStepSetOrder(void* arkode_mem, int ord)

   Specifies the order of accuracy for the integration method in the
   ERK timestepper module. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *ord* -- requested order of accuracy.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The allowed values are :math:`2 \le` *ord* :math:`\le
   8`.  Any illegal input will result in the default value of 4. 
   
   Since *ord* affects the memory requirements for the internal
   ARKode memory block, it cannot be changed after the first call to
   :c:func:`ARKode()`, unless :c:func:`ERKStepReInit()` is called.



.. c:function:: int ARKStepSetImEx(void* arkode_mem)

   Specifies that both the implicit and explicit portions
   of problem are enabled, and to use an additive Runge Kutta method.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This is automatically deduced when neither of the function
   pointers *fe* or *fi* passed to :c:func:`ARKStepCreate()` are
   ``NULL``, but may be set directly by the user if desired.



.. c:function:: int ARKStepSetExplicit(void* arkode_mem)

   Specifies that the implicit portion of problem is disabled,
   and to use an explicit RK method in the ARK timestepper module.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This is automatically deduced when the function pointer `fi`
   passed to :c:func:`ARKStepCreate()` is ``NULL``, but may be set
   directly by the user if desired.

   If the problem is posed in explicit form, i.e. :math:`\dot{y} =
   f(t,y)`, then we recommend that the ERKStep timestepper module be
   used instead.


.. c:function:: int ARKStepSetImplicit(void* arkode_mem)

   Specifies that the explicit portion of problem is disabled,
   and to use a diagonally implicit RK method.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This is automatically deduced when the function pointer `fe`
   passed to :c:func:`ARKStepCreate()` is ``NULL``, but may be set
   directly by the user if desired.



.. c:function:: int ARKStepSetARKTables(void* arkode_mem, int s, int q, int p, realtype* ci, realtype* ce, realtype* Ai, realtype* Ae, realtype* bi, realtype* be, realtype* b2i, realtype* b2e)

   Specifies a customized Butcher table pair for the additive RK method.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *s* -- number of stages in the RK method.
      * *q* -- global order of accuracy for the RK method.
      * *p* -- global order of accuracy for the embedded RK method.
      * *ci* -- array (of length *s*) of stage times for the implicit RK method.
      * *ce* -- array (of length *s*) of stage times for the explicit RK method.
      * *Ai* -- array of coefficients defining the implicit RK stages.  This should
        be stored as a 1D array of size *s*s*, in row-major order.
      * *Ae* -- array of coefficients defining the explicit RK stages.  This should
        be stored as a 1D array of size *s*s*, in row-major order.
      * *bi* -- array of implicit coefficients (of length *s*) defining the time step solution.
      * *be* -- array of explicit coefficients (of length *s*) defining the time step solution.
      * *b2i* -- array of implicit coefficients (of length *s*) defining the embedded solution.
      * *b2e* -- array of explicit coefficients (of length *s*) defining the embedded solution.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``   
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This automatically calls :c:func:`ARKStepSetImEx()`.
   
   No error checking is performed to ensure that either *p* or *q*
   correctly describe the coefficients that were input.
   
   Error checking is performed on both *Ai* and *Ae* to ensure
   that they specify DIRK and ERK methods, respectively.  
   
   If either the inputs *b2i* or *b2e* are set to ``NULL``, ARKode
   will run in fixed-step mode (see :c:func:`ARKodeSetFixedStep()`);
   if called in this manner the user *must* call either
   :c:func:`ARKodeSetFixedStep()` or :c:func:`ARKodeSetInitStep()` to
   set the desired time step size.



.. c:function:: int ARKStepSetERKTable(void* arkode_mem, int s, int q, int p, realtype* c, realtype* A, realtype* b, realtype* bembed)

   Specifies a customized Butcher table for the explicit portion of the system.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *s* -- number of stages in the RK method.
      * *q* -- global order of accuracy for the RK method.
      * *p* -- global order of accuracy for the embedded RK method.
      * *c* -- array (of length *s*) of stage times for the RK method.
      * *A* -- array of coefficients defining the RK stages.  This should
        be stored as a 1D array of size *s*s*, in row-major order.
      * *b* -- array of coefficients (of length *s*) defining the time step solution.
      * *bembed* -- array of coefficients (of length *s*) defining the embedded solution.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This automatically calls :c:func:`ARKStepSetExplicit()`.
   
   No error checking is performed to ensure that either *p* or *q*
   correctly describe the coefficients that were input.
   
   Error checking is performed to ensure that *A* is strictly
   lower-triangular (i.e. that it specifies an ERK method).
   
   An input *bembed* of ``NULL`` will signal that ARKode will run in
   fixed-step mode (see :c:func:`ARKodeSetFixedStep()`); if called in
   this manner the user *must* call either
   :c:func:`ARKodeSetFixedStep()` or :c:func:`ARKodeSetInitStep()` to
   set the desired time step size.

   If the problem is posed in explicit form, i.e. :math:`\dot{y} =
   f(t,y)`, then we recommend that the ERKStep timestepper module be
   used instead.

   

.. c:function:: int ERKStepSetERKTable(void* arkode_mem, int s, int q, int p, realtype* c, realtype* A, realtype* b, realtype* bembed)

   Specifies a customized Butcher table for the ERK timestepping
   module to use.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *s* -- number of stages in the RK method.
      * *q* -- global order of accuracy for the RK method.
      * *p* -- global order of accuracy for the embedded RK method.
      * *c* -- array (of length *s*) of stage times for the RK method.
      * *A* -- array of coefficients defining the RK stages.  This should
        be stored as a 1D array of size *s*s*, in row-major order.
      * *b* -- array of coefficients (of length *s*) defining the time step solution.
      * *bembed* -- array of coefficients (of length *s*) defining the embedded solution.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** No error checking is performed to ensure that either *p*
   or *q* correctly describe the coefficients that were input.
   
   Error checking is performed to ensure that *A* is strictly
   lower-triangular (i.e. that it specifies an ERK method).
   
   An input *bembed* of ``NULL`` will signal that ARKode will run in
   fixed-step mode (see :c:func:`ARKodeSetFixedStep()`); if called in
   this manner the user *must* call either
   :c:func:`ARKodeSetFixedStep()` or :c:func:`ARKodeSetInitStep()` to
   set the desired time step size.

   

.. c:function:: int ARKStepSetIRKTable(void* arkode_mem, int s, int q, int p, realtype* c, realtype* A, realtype* b, realtype* bembed)

   Specifies a customized Butcher table for the implicit portion of the system.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *s* -- number of stages in the RK method.
      * *q* -- global order of accuracy for the RK method.
      * *p* -- global order of accuracy for the embedded RK method.
      * *c* -- array (of length *s*) of stage times for the RK method.
      * *A* -- array of coefficients defining the RK stages.  This should
        be stored as a 1D array of size *s*s*, in row-major order.
      * *b* -- array of coefficients (of length *s*) defining the time step solution.
      * *bembed* -- array of coefficients (of length *s*) defining the embedded solution.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This automatically calls :c:func:`ARKStepSetImplicit()`.
   
   No error checking is performed to ensure that either *p* or *q*
   correctly describe the coefficients that were input.
   
   Error checking is performed to ensure that *A* is 
   lower-triangular with a nonzero value on at least one of the
   diagonal entries (i.e. that it specifies a DIRK method).
   
   An input *bembed* of ``NULL`` will signal that ARKode will run in
   fixed-step mode (see :c:func:`ARKodeSetFixedStep()`); if called in
   this manner the user *must* call either
   :c:func:`ARKodeSetFixedStep()` or :c:func:`ARKodeSetInitStep()` to
   set the desired time step size.



.. c:function:: int ARKStepSetARKTableNum(void* arkode_mem, int itable, int etable)

   Indicates to use specific built-in Butcher tables for the ImEx system.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *itable* -- index of the DIRK Butcher table.
      * *etable* -- index of the ERK Butcher table.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Both *itable* and *etable* should match an existing
   implicit/explicit pair, listed in the section :ref:`Butcher.additive`.   
   Error-checking is performed to ensure that the tables exist.
   Subsequent error-checking is automatically performed to ensure that
   the tables' stage times and solution coefficients match.  

   This automatically calls :c:func:`ARKStepSetImEx()`. 



.. c:function:: int ARKStepSetERKTableNum(void* arkode_mem, int etable)

   Indicates to use a specific built-in Butcher table for explicit
   integration of the problem by the ARK timestepping module.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *etable* -- index of the Butcher table.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** *etable* should match an existing explicit method from
   the section :ref:`Butcher.explicit`.  Error-checking is performed
   to ensure that the table exists, and is not implicit.  
   
   This automatically calls :c:func:`ARKStepSetExplicit()`. 

   If the problem is posed in explicit form, i.e. :math:`\dot{y} =
   f(t,y)`, then we recommend that the ERKStep timestepper module be
   used instead.


   
.. c:function:: int ERKStepSetERKTableNum(void* arkode_mem, int etable)

   Indicates to use a specific built-in Butcher table for use by the
   ERK timestepping module.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *etable* -- index of the Butcher table.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** *etable* should match an existing explicit method from
   the section :ref:`Butcher.explicit`.  Error-checking is performed
   to ensure that the table exists, and is not implicit.  
   


.. c:function:: int ARKStepSetIRKTableNum(void* arkode_mem, int itable)

   Indicates to use a specific built-in Butcher table for implicit
   integration of the problem.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *itable* -- index of the Butcher table.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** *itable* should match an existing implicit method from
   the section :ref:`Butcher.implicit`.  Error-checking is performed
   to ensure that the table exists, and is not explicit.  
   
   This automatically calls :c:func:`ARKStepSetImplicit()`. 







.. _CInterface.ARKodeAdaptivityInputTable:

Optional inputs for time step adaptivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of ARKode's time step adaptivity
algorithm, including how each of the parameters below is used within
the code, is provided in the section :ref:`Mathematics.Adaptivity`.


.. cssclass:: table-bordered

==============================================  ======================================  ========
Optional input                                  Function name                           Default
==============================================  ======================================  ========
Set a custom time step adaptivity function      :c:func:`ARKStepSetAdaptivityFn()`      internal
                                                :c:func:`ERKStepSetAdaptivityFn()`      internal
Choose an existing time step adaptivity method  :c:func:`ARKStepSetAdaptivityMethod()`  0
                                                :c:func:`ERKStepSetAdaptivityMethod()`  0
Explicit stability safety factor                :c:func:`ARKStepSetCFLFraction()`       0.5
                                                :c:func:`ERKStepSetCFLFraction()`       0.5
Time step error bias factor                     :c:func:`ARKStepSetErrorBias()`         1.5
                                                :c:func:`ERKStepSetErrorBias()`         1.5
Bounds determining no change in step size       :c:func:`ARKStepSetFixedStepBounds()`   1.0  1.5
                                                :c:func:`ERKStepSetFixedStepBounds()`   1.0  1.5
Maximum step growth factor on convergence fail  :c:func:`ARKStepSetMaxCFailGrowth()`    0.25
Maximum step growth factor on error test fail   :c:func:`ARKStepSetMaxEFailGrowth()`    0.3
                                                :c:func:`ERKStepSetMaxEFailGrowth()`    0.3
Maximum first step growth factor                :c:func:`ARKStepSetMaxFirstGrowth()`    10000.0
                                                :c:func:`ERKStepSetMaxFirstGrowth()`    10000.0
Maximum general step growth factor              :c:func:`ARKStepSetMaxGrowth()`         20.0
                                                :c:func:`ERKStepSetMaxGrowth()`         20.0
Time step safety factor                         :c:func:`ARKStepSetSafetyFactor()`      0.96
                                                :c:func:`ERKStepSetSafetyFactor()`      0.96
Error fails before MaxEFailGrowth takes effect  :c:func:`ARKStepSetSmallNumEFails()`    2
                                                :c:func:`ERKStepSetSmallNumEFails()`    2
Explicit stability function                     :c:func:`ARKStepSetStabilityFn()`       internal
                                                :c:func:`ERKStepSetStabilityFn()`       internal
==============================================  ======================================  ========



.. c:function:: int ARKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)

   Sets a user-supplied time-step adaptivity function.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hfun* -- name of user-supplied adaptivity function.
      * *h_data* -- pointer to user data passed to *hfun* every time
        it is called.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This function should focus on accuracy-based time step
   estimation; for stability based time steps the function
   :c:func:`ARKStepSetStabilityFn()` should be used instead.


      
.. c:function:: int ERKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)

   Sets a user-supplied time-step adaptivity function.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hfun* -- name of user-supplied adaptivity function.
      * *h_data* -- pointer to user data passed to *hfun* every time
        it is called.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This function should focus on accuracy-based time step
   estimation; for stability based time steps the function
   :c:func:`ERKStepSetStabilityFn()` should be used instead.


      
.. c:function:: int ARKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault, int pq, realtype* adapt_params)

   Specifies the method (and associated parameters) used for time step adaptivity.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *imethod* -- accuracy-based adaptivity method choice 
        (0 :math:`\le` `imethod` :math:`\le` 5): 
        0 is PID, 1 is PI, 2 is I, 3 is explicit Gustafsson, 4 is
        implicit Gustafsson, and 5 is the ImEx Gustafsson.
      * *idefault* -- flag denoting whether to use default adaptivity
	parameters (1), or that they will be supplied in the
	*adapt_params* argument (0).
      * *pq* -- flag denoting whether to use the embedding order of
	accuracy *p* (0) or the method order of accuracy *q* (1)
	within the adaptivity algorithm.  *p* is the ARKode default.
      * *adapt_params[0]* -- :math:`k_1` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[1]* -- :math:`k_2` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[2]* -- :math:`k_3` parameter within accuracy-based adaptivity algorithms.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** If custom parameters are supplied, they will be checked
   for validity against published stability intervals.  If other
   parameter values are desired, it is recommended to instead provide
   a custom function through a call to :c:func:`ARKStepSetAdaptivityFn()`.


      
.. c:function:: int ERKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault, int pq, realtype* adapt_params)

   Specifies the method (and associated parameters) used for time step adaptivity.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *imethod* -- accuracy-based adaptivity method choice 
        (0 :math:`\le` `imethod` :math:`\le` 5): 
        0 is PID, 1 is PI, 2 is I, 3 is explicit Gustafsson, 4 is
        implicit Gustafsson, and 5 is the ImEx Gustafsson.
      * *idefault* -- flag denoting whether to use default adaptivity
	parameters (1), or that they will be supplied in the
	*adapt_params* argument (0).
      * *pq* -- flag denoting whether to use the embedding order of
	accuracy *p* (0) or the method order of accuracy *q* (1)
	within the adaptivity algorithm.  *p* is the ARKode default.
      * *adapt_params[0]* -- :math:`k_1` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[1]* -- :math:`k_2` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[2]* -- :math:`k_3` parameter within accuracy-based adaptivity algorithms.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** If custom parameters are supplied, they will be checked
   for validity against published stability intervals.  If other
   parameter values are desired, it is recommended to instead provide
   a custom function through a call to :c:func:`ERKStepSetAdaptivityFn()`.


      
.. c:function:: int ARKStepSetCFLFraction(void* arkode_mem, realtype cfl_frac)

   Specifies the fraction of the estimated explicitly stable step to use.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *cfl_frac* -- maximum allowed fraction of explicitly stable step (default is 0.5).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default
   value.  
   

      
.. c:function:: int ERKStepSetCFLFraction(void* arkode_mem, realtype cfl_frac)

   Specifies the fraction of the estimated explicitly stable step to use.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *cfl_frac* -- maximum allowed fraction of explicitly stable step (default is 0.5).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default
   value.  
   

      
.. c:function:: int ARKStepSetErrorBias(void* arkode_mem, realtype bias)

   Specifies the bias to be applied to the error estimates within
   accuracy-based adaptivity strategies.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *bias* -- bias applied to error in accuracy-based time
        step estimation (default is 1.5).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value below 1.0 will imply a reset to the default value.  


      
.. c:function:: int ERKStepSetErrorBias(void* arkode_mem, realtype bias)

   Specifies the bias to be applied to the error estimates within
   accuracy-based adaptivity strategies.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *bias* -- bias applied to error in accuracy-based time
        step estimation (default is 1.5).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value below 1.0 will imply a reset to the default value.  


      
.. c:function:: int ARKStepSetFixedStepBounds(void* arkode_mem, realtype lb, realtype ub)

   Specifies the step growth interval in which the step size will remain unchanged.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lb* -- lower bound on window to leave step size fixed (default is 1.0).
      * *ub* -- upper bound on window to leave step size fixed (default is 1.5).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any interval *not* containing 1.0 will imply a reset to the default values.
   

      
.. c:function:: int ERKStepSetFixedStepBounds(void* arkode_mem, realtype lb, realtype ub)

   Specifies the step growth interval in which the step size will remain unchanged.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lb* -- lower bound on window to leave step size fixed (default is 1.0).
      * *ub* -- upper bound on window to leave step size fixed (default is 1.5).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any interval *not* containing 1.0 will imply a reset to the default values.
   

      
.. c:function:: int ARKStepSetMaxCFailGrowth(void* arkode_mem, realtype etacf)

   Specifies the maximum step size growth factor upon an algebraic
   solver convergence failure on a stage solve within a step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *etacf* -- time step reduction factor on a nonlinear solver
        convergence failure (default is 0.25).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value outside the interval :math:`(0,1]` will imply a reset to the default value.



.. c:function:: int ARKStepSetMaxEFailGrowth(void* arkode_mem, realtype etamxf)

   Specifies the maximum step size growth factor upon multiple successive
   accuracy-based error failures in the solver.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *etamxf* -- time step reduction factor on multiple error fails (default is 0.3).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value outside the interval :math:`(0,1]` will imply a reset to the default value.
   

      
.. c:function:: int ERKStepSetMaxEFailGrowth(void* arkode_mem, realtype etamxf)

   Specifies the maximum step size growth factor upon multiple successive
   accuracy-based error failures in the solver.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *etamxf* -- time step reduction factor on multiple error fails (default is 0.3).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value outside the interval :math:`(0,1]` will imply a reset to the default value.
   

      
.. c:function:: int ARKStepSetMaxFirstGrowth(void* arkode_mem, realtype etamx1)

   Specifies the maximum allowed step size change following the very
   first integration step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *etamx1* -- maximum allowed growth factor after the first time
        step (default is 10000.0).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value :math:`\le 1.0` will imply a reset to the default value.



.. c:function:: int ERKStepSetMaxFirstGrowth(void* arkode_mem, realtype etamx1)

   Specifies the maximum allowed step size change following the very
   first integration step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *etamx1* -- maximum allowed growth factor after the first time
        step (default is 10000.0).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value :math:`\le 1.0` will imply a reset to the default value.



.. c:function:: int ARKStepSetMaxGrowth(void* arkode_mem, realtype mx_growth)

   Specifies the maximum growth of the step size between consecutive
   steps in the integration process.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *growth* -- maximum allowed growth factor between consecutive time steps (default is 20.0).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value :math:`\le 1.0` will imply a reset to the default
   value.  



.. c:function:: int ERKStepSetMaxGrowth(void* arkode_mem, realtype mx_growth)

   Specifies the maximum growth of the step size between consecutive
   steps in the integration process.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *growth* -- maximum allowed growth factor between consecutive time steps (default is 20.0).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any value :math:`\le 1.0` will imply a reset to the default
   value.  



.. c:function:: int ARKStepSetSafetyFactor(void* arkode_mem, realtype safety)

   Specifies the safety factor to be applied to the accuracy-based
   estimated step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *safety* -- safety factor applied to accuracy-based time step (default is 0.96).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default
   value.  


      
.. c:function:: int ERKStepSetSafetyFactor(void* arkode_mem, realtype safety)

   Specifies the safety factor to be applied to the accuracy-based
   estimated step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *safety* -- safety factor applied to accuracy-based time step (default is 0.96).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default
   value.  


      
.. c:function:: int ARKStepSetSmallNumEFails(void* arkode_mem, int small_nef)

   Specifies the threshold for "multiple" successive error failures
   before the *etamxf* parameter from
   :c:func:`ARKStepSetMaxEFailGrowth()` is applied.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *small_nef* -- bound to determine 'multiple' for *etamxf* (default is 2).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ERKStepSetSmallNumEFails(void* arkode_mem, int small_nef)

   Specifies the threshold for "multiple" successive error failures
   before the *etamxf* parameter from
   :c:func:`ERKStepSetMaxEFailGrowth()` is applied.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *small_nef* -- bound to determine 'multiple' for *etamxf* (default is 2).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ARKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)

   Sets the problem-dependent function to estimate a stable
   time step size for the explicit portion of the ODE system.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *EStab* -- name of user-supplied stability function.
      * *estab_data* -- pointer to user data passed to *EStab* every time
        it is called.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This function should return an estimate of the absolute
   value of the maximum stable time step for the explicit portion of
   the ODE system.  It is not required, since accuracy-based
   adaptivity may be sufficient for retaining stability, but this can
   be quite useful for problems where the explicit right-hand side
   function :math:`f_E(t,y)` may contain stiff terms.



.. c:function:: int ERKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)

   Sets the problem-dependent function to estimate a stable
   time step size for the explicit portion of the ODE system.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *EStab* -- name of user-supplied stability function.
      * *estab_data* -- pointer to user data passed to *EStab* every time
        it is called.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This function should return an estimate of the absolute
   value of the maximum stable time step for the the ODE system.  It
   is not required, since accuracy-based adaptivity may be sufficient
   for retaining stability, but this can be quite useful for problems
   where the right-hand side function :math:`f(t,y)` may contain stiff
   terms.







.. _CInterface.CInterface.ARKodeSolverInputTable:

Optional inputs for implicit stage solves (ARKStep module only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation for the nonlinear solver strategies used
by the ARKStep timestepping module, including how each of the
parameters below is used within the code, is provided in the section
:ref:`Mathematics.Nonlinear`. 


.. cssclass:: table-bordered

=============================================  =========================================  ============
Optional input                                 Function name                              Default
=============================================  =========================================  ============
Specify use of the fixed-point stage solver    :c:func:`ARKStepSetFixedPoint()`           ``SUNFALSE``
Specify use of the Newton stage solver         :c:func:`ARKStepSetNewton()`               ``SUNTRUE``
Specify linearly implicit :math:`f_I`          :c:func:`ARKStepSetLinear()`               ``SUNFALSE``
Specify nonlinearly implicit :math:`f_I`       :c:func:`ARKStepSetNonlinear()`            ``SUNTRUE``
Implicit predictor method                      :c:func:`ARKStepSetPredictorMethod()`      0
Maximum number of nonlinear iterations         :c:func:`ARKStepSetMaxNonlinIters()`       3
Coefficient in the nonlinear convergence test  :c:func:`ARKStepSetNonlinConvCoef()`       0.1
Nonlinear convergence rate constant            :c:func:`ARKStepSetNonlinCRDown()`         0.3
Nonlinear residual divergence ratio            :c:func:`ARKStepSetNonlinRDiv()`           2.3
Max change in step signaling new :math:`J`     :c:func:`ARKStepSetDeltaGammaMax()`        0.2
Max steps between calls to new :math:`J`       :c:func:`ARKStepSetMaxStepsBetweenLSet()`  20
Maximum number of convergence failures         :c:func:`ARKStepSetMaxConvFails()`         10
=============================================  =========================================  ============




.. c:function:: int ARKStepSetFixedPoint(void* arkode_mem, long int fp_m)

   Specifies that the implicit portion of the problem should be solved
   using the accelerated fixed-point solver instead of the modified
   Newton iteration, and provides the maximum dimension of the
   acceleration subspace.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *fp_m* -- number of vectors to store within the Anderson
        acceleration subspace.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Since the accelerated fixed-point solver has a slower
   rate of convergence than the Newton iteration (but each iteration
   is typically much more efficient), it is recommended that the
   maximum nonlinear correction iterations be increased through a call
   to :c:func:`ARKStepSetMaxNonlinIters()`. 

   Since *fp_m* affects the memory requirements for the internal
   ARKode memory block, it cannot be changed after the first call to
   :c:func:`ARKode()`, unless :c:func:`ARKStepReInit()` is called.

      

.. c:function:: int ARKStepSetNewton(void* arkode_mem)

   Specifies that the implicit portion of the problem should be solved
   using the modified Newton solver.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This is the default behavior of ARKStep, so the function
   is primarily useful to undo a previous call to :c:func:`ARKStepSetFixedPoint()`.



.. c:function:: int ARKStepSetLinear(void* arkode_mem, int timedepend)

   Specifies that the implicit portion of the problem is linear.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *timedepend* -- flag denoting whether the Jacobian of
	:math:`f_I(t,y)` is time-dependent (1) or not (0).
	Alternately, when using an iterative linear solver this flag
	denotes time dependence of the preconditioner.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Tightens the linear solver tolerances and takes only a
   single Newton iteration.  Calls :c:func:`ARKStepSetDeltaGammaMax()`
   to enforce Jacobian recomputation when the step size ratio changes
   by more than 100 times the unit roundoff (since nonlinear
   convergence is not tested).  Only applicable when used in
   combination with the modified Newton iteration (not the fixed-point
   solver).



.. c:function:: int ARKStepSetNonlinear(void* arkode_mem)

   Specifies that the implicit portion of the problem is nonlinear.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** This is the default behavior of ARKode, so the function
   is primarily useful to undo a previous call to
   :c:func:`ARKStepSetLinear()`.  Calls
   :c:func:`ARKStepSetDeltaGammaMax()` to reset the step size ratio
   threshold to the default value. 



.. c:function:: int ARKStepSetPredictorMethod(void* arkode_mem, int method)

   Specifies the method to use for predicting implicit solutions.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *method* -- method choice (0 :math:`\le` *method* :math:`\le` 4): 

        * 0 is the trivial predictor, 

        * 1 is the maximum order (dense output) predictor, 

	* 2 is the variable order predictor, that decreases the
	  polynomial degree for more distant RK stages,

        * 3 is the cutoff order predictor, that uses the maximum order
	  for early RK stages, and a first-order predictor for distant
	  RK stages,

        * 4 is the bootstrap predictor, that uses a second-order
	  predictor based on only information within the current step.
   
        * 5 is the minimum correction predictor, that uses all
          preceding stage information within the current step for
          prediction.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value is 0.  If *method* is set to an
   undefined value, this default predictor will be used.



.. c:function:: int ARKStepSetMaxNonlinIters(void* arkode_mem, int maxcor)

   Specifies the maximum number of nonlinear solver
   iterations permitted per RK stage within each time step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *maxcor* -- maximum allowed solver iterations per stage :math:`(>0)`.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value is 3; set *maxcor* :math:`\le 0`
   to specify this default.



.. c:function:: int ARKStepSetNonlinConvCoef(void* arkode_mem, realtype nlscoef)

   Specifies the safety factor used within the nonlinear
   solver convergence test.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nlscoef* -- coefficient in nonlinear solver convergence test :math:`(>0.0)`.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value is 0.1; set *nlscoef* :math:`\le 0`
   to specify this default.



.. c:function:: int ARKStepSetNonlinCRDown(void* arkode_mem, realtype crdown)

   Specifies the constant used in estimating the nonlinear solver convergence rate.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *crdown* -- nonlinear convergence rate estimation constant (default is 0.3).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ARKStepSetNonlinRDiv(void* arkode_mem, realtype rdiv)

   Specifies the nonlinear correction threshold beyond which the
   iteration will be declared divergent.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *rdiv* -- tolerance on nonlinear correction size ratio to
	declare divergence (default is 2.3).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ARKStepSetDeltaGammaMax(void* arkode_mem, realtype dgmax)

   Specifies a scaled step size ratio tolerance, beyond which the
   linear solver setup routine will be signaled.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *dgmax* -- tolerance on step size ratio change before calling
        linear solver setup routine (default is 0.2).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:**  Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ARKStepSetMaxStepsBetweenLSet(void* arkode_mem, int msbp)

   Specifies the frequency of calls to the linear solver setup
   routine.  Positive values specify the number of time steps between
   setup calls; negative values force recomputation at each Newton
   step; zero values reset to the default.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *msbp* -- maximum number of time steps between linear solver
	setup calls, or flag to force recomputation at each Newton
	iteration (default is 20).
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
   


.. c:function:: int ARKStepSetMaxConvFails(void* arkode_mem, int maxncf)

   Specifies the maximum number of nonlinear solver convergence
   failures permitted during one step, before ARKode will return with
   an error.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *maxncf* -- maximum allowed nonlinear solver convergence failures
        per step :math:`(>0)`.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default value is 10; set *maxncf* :math:`\le 0`
   to specify this default.  

   Upon each convergence failure, ARKStep will first call the Jacobian
   setup routine and try again (if a Newton method is used).  If a
   convergence failure still occurs, the time step size is reduced by
   the factor *etacf* (set within :c:func:`ARKStepSetMaxCFailGrowth()`).





.. _CInterface.ARKDlsInputs:


Direct linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of ARKode's direct linear solver methods
is provided in the section :ref:`Mathematics.Linear`.


Table: Optional inputs for ARKDLS
"""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==========================  ================================  =============
Optional input              Function name                     Default
==========================  ================================  =============
Jacobian function           :c:func:`ARKDlsSetJacFn()`        ``DQ``
Mass matrix function        :c:func:`ARKDlsSetMassFn()`       none
==========================  ================================  =============

The ARKDLS solver interface needs a function to compute an
approximation to the Jacobian matrix :math:`J(t,y)`. This function
must be of type :c:func:`ARKDlsJacFn()`.  The user can supply a custom
Jacobian function, or if using a dense or banded :math:`J` can use the
default internal difference quotient approximation that comes with the
ARKDLS interface.  To specify a user-supplied Jacobian function *jac*,
ARKDLS provides the function :c:func:`ARKDlsSetJacFn()`. The ARKDLS
interface passes the user data pointer to the Jacobian function. This
allows the user to create an arbitrary structure with relevant problem
data and access it during the execution of the user-supplied Jacobian
function, without using global data in the program. The user
data pointer may be specified through :c:func:`ARKodeSetUserData()`.

..
   Similarly, if the ODE system involves a non-identity mass matrix,
   :math:`M\ne I`, the ARKDLS interface needs a function to compute an
   approximation to the mass matrix :math:`M(t)`. There is no default
Similarly, if the ODE system involves a non-identity mass matrix,
:math:`M\ne I`, the ARKDLS interface needs a function to compute an
approximation to the mass matrix :math:`M`. There is no default
difference quotient approximation, so this routine must be supplied by
the user. This function must be of type :c:func:`ARKDlsMassFn()`, and
should be set using the function :c:func:`ARKDlsSetMassFn()`.  We note
that the ARKDLS solver passes the user data pointer to the mass matrix
function. This allows the user to create an arbitrary structure with
relevant problem data and access it during the execution of the
user-supplied mass matrix function, without using global data in the
program. The pointer user data may be specified through
:c:func:`ARKodeSetUserData()`.



.. c:function:: int ARKDlsSetJacFn(void* arkode_mem, ARKDlsJacFn jac)

   Specifies the Jacobian approximation routine to
   be used for the ARKDLS interface.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *jac* -- name of user-supplied Jacobian approximation function.
   
   **Return value:** 
      * *ARKDLS_SUCCESS*  if successful
      * *ARKDLS_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** This routine must be called after the ARKDLS linear
   solver interface has been initialized through a call to 
   :c:func:`ARKDlsSetLinearSolver()`.

   By default, ARKDLS uses an internal difference quotient
   function for dense and band matrices.  If ``NULL`` is passed in for
   *jac*, this default is used. An error will occur if no *jac*
   is supplied when using a sparse matrix.
  
   The function type :c:func:`ARKDlsJacFn()` is described in the section
   :ref:`CInterface.UserSupplied`.



.. c:function:: int ARKDlsSetMassFn(void* arkode_mem, ARKDlsMassFn mass)

   Specifies the mass matrix approximation routine to be used for the
   ARKDLS interface. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *mass* -- name of user-supplied mass matrix approximation function.
   
   **Return value:** 
      * *ARKDLS_SUCCESS*  if successful
      * *ARKDLS_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARKDLS_MASSMEM_NULL* if the mass matrix solver memory was ``NULL``
   
   **Notes:** This routine must be called after the ARKDLS mass matrix
   solver interface has been initialized through a call to 
   :c:func:`ARKDlsSetMassLinearSolver()`. 
   
   Since there is no default difference quotient function for mass
   matrices, *mass* must be non-``NULL``.
  
   The function type :c:func:`ARKDlsMassFn()` is described in the section
   :ref:`CInterface.UserSupplied`.




.. _CInterface.ARKSpilsInputs:

Iterative linear solvers optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As described in the section :ref:`Mathematics.Linear`, when using
the ARKSPILS iterative linear solver interface, a user may supply a
preconditioning operator to aid in solution of the system.  This
operator consists of two user-supplied functions, *psetup* and
*psolve*, that are supplied to ARKode using either the function
:c:func:`ARKSpilsSetPreconditioner()` (for preconditioning the
Newton system), or the function
:c:func:`ARKSpilsSetMassPreconditioner()` (for preconditioning the
mass matrix system).  The *psetup* function should handle evaluation
and preprocessing of any Jacobian or mass-matrix data needed by the
user's preconditioner solve function, *psolve*.  The user data pointer
received through :c:func:`ARKodeSetUserData()` (or a pointer to
``NULL`` if user data was not specified) is passed to the *psetup* and
*psolve* functions.  This allows the user to create an arbitrary
structure with relevant problem data and access it during the
execution of the user-supplied preconditioner functions without using
global data in the program.  If preconditioning is supplied for both
the Newton and mass matrix linear systems, it is expected that the
user will supply different *psetup* and *psolve* function for each.

Additionally, when solving the Newton linear systems, the ARKSPILS
interface requires a *jtimes* function to compute an approximation to the
product between the Jacobian matrix :math:`J(t,y)` and a vector
:math:`v`. The user can supply a custom Jacobian-times-vector
approximation function, or use the default internal difference
quotient function that comes with the ARKSPILS interface.  A
user-defined Jacobian-vector function must be of type
:c:type:`ARKSpilsJacTimesVecFn` and can be specified through a call
to :c:func:`ARKSpilsSetJacTimes()` (see the section
:ref:`CInterface.UserSupplied` for specification details).  As with the
user-supplied preconditioner functions, the evaluation and
processing of any Jacobian-related data needed by the user's
Jacobian-times-vector function is done in the optional user-supplied
function of type :c:type:`ARKSpilsJacTimesSetupFn` (see the section
:ref:`CInterface.UserSupplied` for specification details).  As with
the preconditioner functions, a pointer to the user-defined 
data structure, *user_data*, specified through
:c:func:`ARKodeSetUserData()` (or a ``NULL`` pointer otherwise) is
passed to the Jacobian-times-vector setup and product functions each
time they are called.

..
   Similarly, if a problem involves a non-identity mass matrix,
   :math:`M\ne I`, then the ARKSPILS solvers require a *mtimes* function
   to compute an approximation to the product between the mass matrix
   :math:`M(t)` and a vector :math:`v`.  This function must be 
Similarly, if a problem involves a non-identity mass matrix,
:math:`M\ne I`, then the ARKSPILS solvers require a *mtimes* function
to compute an approximation to the product between the mass matrix
:math:`M` and a vector :math:`v`.  This function must be 
user-supplied, since there is no default value.  *mtimes* must be 
of type :c:func:`ARKSpilsMassTimesVecFn()`, and can be specified
through a call to the  :c:func:`ARKSpilsSetMassTimes()` routine.
As with the user-supplied preconditioner functions, the evaluation and
processing of any Jacobian-related data needed by the user's
mass-matrix-times-vector function is done in the optional user-supplied
function of type :c:type:`ARKSpilsMassTimesSetupFn` (see the section
:ref:`CInterface.UserSupplied` for specification details).

Finally, as described in the section :ref:`Mathematics.Error.Linear`, the
ARKSPILS interface requires that iterative linear solvers stop when
the norm of the preconditioned residual satisfies

.. math::
   \|r\| \le \frac{\epsilon_L \epsilon}{10}

where the default :math:`\epsilon_L = 0.05`, which may be modified by
the user through the :c:func:`ARKSpilsSetEpsLin()` function.



.. _CInterface.ARKSpilsInputTable:

Table: Optional inputs for ARKSPILS
"""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==================================================  =========================================  ==================
Optional input                                      Function name                              Default
==================================================  =========================================  ==================
:math:`Jv` functions (*jtimes* and *jtsetup*)       :c:func:`ARKSpilsSetJacTimes()`            DQ,  none
Newton linear and nonlinear tolerance ratio         :c:func:`ARKSpilsSetEpsLin()`              0.05
Newton preconditioning functions                    :c:func:`ARKSpilsSetPreconditioner()`      ``NULL``, ``NULL``
:math:`Mv` functions (*mtimes* and *mtsetup*)       :c:func:`ARKSpilsSetMassTimes()`           none, none
Mass matrix linear and nonlinear tolerance ratio    :c:func:`ARKSpilsSetMassEpsLin()`          0.05
Mass matrix preconditioning functions               :c:func:`ARKSpilsSetMassPreconditioner()`  ``NULL``, ``NULL``
==================================================  =========================================  ==================




.. c:function:: int ARKSpilsSetJacTimes(void* arkode_mem, ARKSpilsJacTimesSetupFn jtsetup, ARKSpilsJacTimesVecFn jtimes)

   Specifies the Jacobian-times-vector setup and product functions. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *jtsetup* -- user-defined Jacobian-vector setup function.
        Pass ``NULL`` if no setup is necessary.
      * *jtimes* -- user-defined Jacobian-vector product function.
   
   **Return value:** 
      * *ARKSPILS_SUCCESS* if successful.
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``.
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKSPILS_ILL_INPUT* if an input has an illegal value.
      * *ARKSPILS_SUNLS_FAIL* if an error occurred when setting up
        the Jacobian-vector product in the ``SUNLinearSolver``
        object used by the ARKSPILS interface. 

   **Notes:** The default is to use an internal finite difference
   quotient for *jtimes* and to leave out *jtsetup*.  If ``NULL`` is
   passed to *jtimes*, these defaults are used.  A user may
   specify non-``NULL`` *jtimes* and ``NULL`` *jtsetup*.

   This function must be called *after* the ARKSPILS system solver
   interface has been initialized through a call to
   :c:func:`ARKSpilsSetLinearSolver()`.
   
   The function types :c:type:`ARKSpilsJacTimesSetupFn` and
   :c:type:`ARKSpilsJacTimesVecFn` is described in the section
   :ref:`CInterface.UserSupplied`. 



.. c:function:: int ARKSpilsSetEpsLin(void* arkode_mem, realtype eplifac)

   Specifies the factor by which the tolerance on the nonlinear
   iteration is multiplied to get a tolerance on the linear
   iteration. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *eplifac* -- linear convergence safety factor :math:`(\ge 0.0)`.
   
   **Return value:** 
      * *ARKSPILS_SUCCESS* if successful.
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``.
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKSPILS_ILL_INPUT* if an input has an illegal value.
   
   **Notes:** Passing a value *eplifac* of 0.0 indicates to use the
   default value of 0.05. 

   This function must be called *after* the ARKSPILS system solver
   interface has been initialized through a call to
   :c:func:`ARKSpilsSetLinearSolver()`.
   


.. c:function:: int ARKSpilsSetPreconditioner(void* arkode_mem, ARKSpilsPrecSetupFn psetup, ARKSpilsPrecSolveFn psolve)

   Specifies the user-supplied preconditioner setup and solve functions.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *psetup* -- user defined preconditioner setup function.  Pass
        ``NULL`` if no setup is needed.
      * *psolve* -- user-defined preconditioner solve function.
   
   **Return value:** 
      * *ARKSPILS_SUCCESS* if successful.
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``.
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKSPILS_ILL_INPUT* if an input has an illegal value.
      * *ARKSPILS_SUNLS_FAIL* if an error occurred when setting up
        preconditioning in the ``SUNLinearSolver`` object used
        by the ARKSPILS interface. 
   
   **Notes:** The default is ``NULL`` for both arguments (i.e. no
   preconditioning).

   This function must be called *after* the ARKSPILS system solver
   interface has been initialized through a call to
   :c:func:`ARKSpilsSetLinearSolver()`.
   
   Both of the function types :c:func:`ARKSpilsPrecSetupFn()` and
   :c:func:`ARKSpilsPrecSolveFn()` are described in the section
   :ref:`CInterface.UserSupplied`. 



.. c:function:: int ARKSpilsSetMassTimes(void* arkode_mem, ARKSpilsMassTimesSetupFn mtsetup, ARKSpilsMassTimesVecFn mtimes, void* mtimes_data)

   Specifies the mass matrix-times-vector setup and product functions. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *mtsetup* -- user-defined mass matrix-vector setup function.
        Pass ``NULL`` if no setup is necessary.
      * *mtimes* -- user-defined mass matrix-vector product function.
      * *mtimes_data* -- a pointer to user data, that will be supplied
        to both the *mtsetup* and *mtimes* functions. 
   
   **Return value:** 
      * *ARKSPILS_SUCCESS* if successful.
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``.
      * *ARKSPILS_MASSMEM_NULL* if the mass matrix solver memory was ``NULL``.
      * *ARKSPILS_ILL_INPUT* if an input has an illegal value.
      * *ARKSPILS_SUNLS_FAIL* if an error occurred when setting up
        the mass-matrix-vector product in the ``SUNLinearSolver``
        object used by the ARKSPILS interface. 

   **Notes:** There is no default finite difference quotient for
   *mtimes*, so if using the ARKSPILS mass matrix solver interface and
   this routine is not called with non-``NULL`` *mtimes*, and error
   will occur.  A user may specify ``NULL`` for *mtsetup*.

   This function must be called *after* the ARKSPILS mass
   matrix solver interface has been initialized through a call to 
   :c:func:`ARKSpilsSetMassLinearSolver()`.
   
   The function types :c:type:`ARKSpilsMassTimesSetupFn` and
   :c:type:`ARKSpilsMassTimesVecFn` are described in the section
   :ref:`CInterface.UserSupplied`. 



.. c:function:: int ARKSpilsSetMassEpsLin(void* arkode_mem, realtype eplifac)

   Specifies the factor by which the tolerance on the nonlinear
   iteration is multiplied to get a tolerance on the mass matrix
   linear iteration. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *eplifac* -- linear convergence safety factor :math:`(\ge 0.0)`.
   
   **Return value:** 
      * *ARKSPILS_SUCCESS* if successful.
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``.
      * *ARKSPILS_MASSMEM_NULL* if the mass matrix solver memory was ``NULL``.
      * *ARKSPILS_ILL_INPUT* if an input has an illegal value.
   
   **Notes:**  This function must be called *after* the ARKSPILS mass
   matrix solver interface has been initialized through a call to 
   :c:func:`ARKSpilsSetMassLinearSolver()`.

   Passing a value *eplifac* of 0.0 indicates to use the default value
   of 0.05.

   

.. c:function:: int ARKSpilsSetMassPreconditioner(void* arkode_mem, ARKSpilsMassPrecSetupFn psetup, ARKSpilsMassPrecSolveFn psolve)

   Specifies the mass matrix preconditioner setup and solve functions.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *psetup* -- user defined preconditioner setup function.  Pass
        ``NULL`` if no setup is to be done.
      * *psolve* -- user-defined preconditioner solve function.
   
   **Return value:** 
      * *ARKSPILS_SUCCESS* if successful.
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``.
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``.
      * *ARKSPILS_ILL_INPUT* if an input has an illegal value.
      * *ARKSPILS_SUNLS_FAIL* if an error occurred when setting up
        preconditioning in the ``SUNLinearSolver`` object used
        by the ARKSPILS interface. 
   
   **Notes:** This function must be called *after* the ARKSPILS mass
   matrix solver interface has been initialized through a call to 
   :c:func:`ARKSpilsSetMassLinearSolver()`.

   The default is ``NULL`` for both arguments (i.e. no
   preconditioning).
    
   Both of the function types :c:func:`ARKSpilsMassPrecSetupFn()` and
   :c:func:`ARKSpilsMassPrecSolveFn()` are described in the section
   :ref:`CInterface.UserSupplied`. 






Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions can be called to set optional inputs to
control the rootfinding algorithm, the mathematics of which are
described in the section :ref:`Mathematics.Rootfinding`.


.. cssclass:: table-bordered

======================================  =======================================  ==================
Optional input                          Function name                            Default
======================================  =======================================  ==================
Direction of zero-crossings to monitor  :c:func:`ARKodeSetRootDirection()`       both
Disable inactive root warnings          :c:func:`ARKodeSetNoInactiveRootWarn()`  enabled
======================================  =======================================  ==================



.. c:function:: int ARKodeSetRootDirection(void* arkode_mem, int* rootdir)

   Specifies the direction of zero-crossings to be located and returned.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *rootdir* -- state array of length *nrtfn*, the number of root
        functions :math:`g_i`, as specified in the call to the function
        :c:func:`ARKodeRootInit()`.  If ``rootdir[i] == 0`` then
	crossing in either direction for :math:`g_i` should be
	reported.  A value of +1 or -1 indicates that the solver
	should report only zero-crossings where :math:`g_i` is
	increasing or decreasing, respectively.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value
   
   **Notes:** The default behavior is to monitor for both zero-crossing
      directions.



.. c:function:: int ARKodeSetNoInactiveRootWarn(void* arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.
  
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
   
   **Notes:** ARKode will not report the initial conditions as a
   possible zero-crossing (assuming that one or more components
   :math:`g_i` are zero at the initial time).  However, if it appears
   that some :math:`g_i` is identically zero at the initial time
   (i.e., :math:`g_i` is zero at the initial time *and* after the
   first step), ARKode will issue a warning which can be disabled with
   this optional input function. 





.. _CInterface.InterpolatedOutput:

Interpolated output function
--------------------------------

An optional function :c:func:`ARKodeGetDky()` is available to obtain
additional values of solution-related quantities.  This function
should only be called after a successful return from
:c:func:`ARKode()`, as it provides interpolated values either of
:math:`y` or of its derivatives (up to the 3rd derivative)
interpolated to any value of :math:`t` in the last internal step taken
by :c:func:`ARKode()`.  Internally, this *dense output* algorithm is
identical to the algorithm used for the maximum order implicit
predictors, described in the section
:ref:`Mathematics.Predictors.Max`, except that derivatives of the
polynomial model may be evaluated upon request. 



.. c:function:: int ARKodeGetDky(void* arkode_mem, realtype t, int k, N_Vector dky)

   Computes the *k*-th derivative of the function
   :math:`y` at the time *t*,
   i.e. :math:`\frac{d^{(k)}}{dt^{(k)}}y(t)`, for values of the
   independent variable satisfying :math:`t_n-h_n \le t \le t_n`, with
   :math:`t_n` as current internal time reached, and :math:`h_n` is
   the last internal step size successfully used by the solver.  The
   user may request *k* in the range {0,1,2,3}.  This routine uses an
   interpolating polynomial of degree *max(dord, k)*, where *dord* is the
   argument provided to :c:func:`ARKodeSetDenseOrder()`.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *t* -- the value of the independent variable at which the
        derivative is to be evaluated.
      * *k* -- the derivative order requested.
      * *dky* -- output vector (must be allocated by the user).
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_BAD_K* if *k* is not in the range {0,1,2,3}.
      * *ARK_BAD_T* if *t* is not in the interval :math:`[t_n-h_n, t_n]`
      * *ARK_BAD_DKY* if the *dky* vector was ``NULL``
      * *ARK_MEM_NULL* if the ARKode memory is ``NULL``
   
   **Notes:** It is only legal to call this function after a successful
   return from :c:func:`ARKode()`.  

   A user may access the values :math:`t_n` and :math:`h_n` via the
   functions :c:func:`ARKodeGetCurrentTime()` and
   :c:func:`ARKodeGetLastStep()`, respectively.




.. _CInterface.OptionalOutputs:

Optional output functions
------------------------------

ARKode provides an extensive set of functions that can be used to
obtain solver performance information.  We organize these into groups:

1. SUNDIALS version information accessor routines are in the subsection
   :ref:`CInterface.SUNVersionInfo`, 
2. General ARKode output routines are in the subsection
   :ref:`CInterface.ARKodeMainOutputs`, 
3. ARKStep-specific general output routines are in the subsection
   :ref:`CInterface.ARKStepMainOutputs`, 
4. ERKStep-specific general output routines are in the subsection
   :ref:`CInterface.ERKStepMainOutputs`, 
5. ARKStep implicit solver output routines are in the subsection
   :ref:`CInterface.ARKodeImplicitSolverOutputs`, 
6. Output routines regarding root-finding results are in the subsection
   :ref:`CInterface.ARKodeRootOutputs`, 
7. Dense linear solver output routines are in the subsection
   :ref:`CInterface.ARKDlsOutputs` and 
8. Iterative linear solver output routines are in the subsection
   :ref:`CInterface.ARKSpilsOutputs`.
9. General usability routines (e.g. to print the current ARKode
   parameters, or output the current Butcher table(s)) are in the
   subsection :ref:`CInterface.ARKodeExtraOutputs`.

Following each table, we elaborate on each function.

Some of the optional outputs, especially the various counters, can be
very useful in determining the efficiency of various methods inside
ARKode.  For example:

* The counters *nsteps*, *nfe_evals*, *nfi_evals* and *nf_evals*
  provide a rough measure of the overall cost of a given run, and can
  be compared between runs with different solver options to suggest
  which set of options is the most efficient. 

* The ratio *nniters/nsteps* measures the performance of the
  nonlinear iteration in solving the nonlinear systems at each stage,
  providing a measure of the degree of nonlinearity in the problem.
  Typical values of this for a Newton solver on a general problem
  range from 1.1 to 1.8.

* When using a Newton nonlinear solver, the ratio *njevals/nniters*
  (in the case of a direct linear solver), and the ratio
  *npevals/nniters* (in the case of an iterative linear solver)
  can measure the overall degree of nonlinearity in the problem,
  since these are updated infrequently, unless the Newton method
  convergence slows.

* When using a Newton nonlinear solver, the ratio *njevals/nniters*
  (when using a direct linear solver), and the ratio
  *nliters/nniters* (when using an iterative linear solver) can
  indicate the quality of the approximate Jacobian or preconditioner being
  used.  For example, if this ratio is larger for a user-supplied
  Jacobian or Jacobian-vector product routine than for the
  difference-quotient routine, it can indicate that the user-supplied
  Jacobian is inaccurate. 

* The ratio *expsteps/accsteps* can measure the quality of the ImEx
  splitting used, since a higher-quality splitting will be dominated
  by accuracy-limited steps.

* The ratio *nsteps/step_attempts* can measure the quality of the
  time step adaptivity algorithm, since a poor algorithm will result
  in more failed steps, and hence a lower ratio.

It is therefore recommended that users retrieve and output these
statistics following each run, and take some time to investigate
alternate solver options that will be more optimal for their
particular problem of interest.
 


.. _CInterface.SUNVersionInfo:

SUNDIALS version information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions provide a way to get SUNDIALS version
information at runtime.


.. c:function:: int SUNDIALSGetVersion(char *version, int len)

   This routine fills a string with SUNDIALS version information.

   **Arguments:**
      * *version* -- character array to hold the SUNDIALS version information.
      * *len* -- allocated length of the *version* character array.
   
   **Return value:**  
      * 0 if successful
      * -1 if the input string is too short to store the SUNDIALS version
   
   **Notes:** An array of 25 characters should be sufficient to hold
   the version information. 



.. c:function:: int SUNDIALSGetVersionNumber(int *major, int *minor, int *patch, char *label, int len)

   This routine The function sets integers for the SUNDIALS major,
   minor, and patch release numbers and fills a string with the
   release label if applicable.

   **Arguments:**
      * *major* -- SUNDIALS release major version number.
      * *minor* -- SUNDIALS release minor version number.
      * *patch* -- SUNDIALS release patch version number.
      * *label* -- string to hold the SUNDIALS release label.
      * *len* -- allocated length of the *label* character array.
   
   **Return value:**  
      * 0 if successful
      * -1 if the input string is too short to store the SUNDIALS label
   
   **Notes:** An array of 10 characters should be sufficient to hold
   the label information. If a label is not used in the release
   version, no information is copied to *label*.


   


.. _CInterface.ARKodeMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ============================================
Optional output                                      Function name
===================================================  ============================================
Size of ARKode real and integer workspaces           :c:func:`ARKodeGetWorkSpace()`
Cumulative number of internal steps                  :c:func:`ARKodeGetNumSteps()`
Actual initial time step size used                   :c:func:`ARKodeGetActualInitStep()`
Step size used for the last successful step          :c:func:`ARKodeGetLastStep()`
Step size to be attempted on the next step           :c:func:`ARKodeGetCurrentStep()`
Current internal time reached by the solver          :c:func:`ARKodeGetCurrentTime()`
Suggested factor for tolerance scaling               :c:func:`ARKodeGetTolScaleFactor()`
Error weight vector for state variables              :c:func:`ARKodeGetErrWeights()`
Residual weight vector                               :c:func:`ARKodeGetResWeights()`
Single accessor to many statistics at once           :c:func:`ARKodeGetStepStats()`
Name of constant associated with a return flag       :c:func:`ARKodeGetReturnFlagName()`
===================================================  ============================================ 




.. c:function:: int ARKodeGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)

   Returns the ARKode real and integer workspace sizes.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lenrw* -- the number of ``realtype`` values in the ARKode workspace.
      * *leniw* -- the number of integer values in the ARKode workspace.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumSteps(void* arkode_mem, long int* nsteps)

   Returns the cumulative number of internal steps taken by
   the solver (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nsteps* -- number of steps taken in the solver.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetActualInitStep(void* arkode_mem, realtype* hinused)

   Returns the value of the integration step size used on the first step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hinused* -- actual value of initial step size.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** Even if the value of the initial integration step was
   specified by the user through a call to
   :c:func:`ARKodeSetInitStep()`, this value may have been changed by
   ARKode to ensure that the step size fell within the prescribed
   bounds :math:`(h_{min} \le h_0 \le h_{max})`, or to satisfy the
   local error test condition, or to ensure convergence of the
   nonlinear solver.



.. c:function:: int ARKodeGetLastStep(void* arkode_mem, realtype* hlast)

   Returns the integration step size taken on the last successful
   internal step. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hlast* -- step size taken on the last internal step.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetCurrentStep(void* arkode_mem, realtype* hcur)

   Returns the integration step size to be attempted on the next internal step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *hcur* -- step size to be attempted on the next internal step.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetCurrentTime(void* arkode_mem, realtype* tcur)

   Returns the current internal time reached by the solver.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *tcur* -- current internal time reached.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetTolScaleFactor(void* arkode_mem, realtype* tolsfac)

   Returns a suggested factor by which the user's
   tolerances should be scaled when too much accuracy has been
   requested for some internal step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *tolsfac* -- suggested scaling factor for user-supplied tolerances.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetErrWeights(void* arkode_mem, N_Vector eweight)

   Returns the current error weight vector.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *eweight* -- solution error weights at the current time.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The user must allocate space for *eweight*, that will be
   filled in by this function.



.. c:function:: int ARKodeGetResWeights(void* arkode_mem, N_Vector rweight)

   Returns the current residual weight vector.  
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *rweight* -- residual error weights at the current time.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The user must allocate space for *rweight*, that will be
   filled in by this function.



.. c:function:: int ARKodeGetStepStats(void* arkode_mem, long int* nsteps, realtype* hinused, realtype* hlast, realtype* hcur, realtype* tcur)

   Returns many of the most useful optional outputs in a single call.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nsteps* -- number of steps taken in the solver.
      * *hinused* -- actual value of initial step size.
      * *hlast* -- step size taken on the last internal step.
      * *hcur* -- step size to be attempted on the next internal step.
      * *tcur* -- current internal time reached.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: char *ARKodeGetReturnFlagName(long int flag)

   Returns the name of the ARKode constant corresponding to *flag*.
   
   **Arguments:**
      * *flag* -- a return flag from an ARKode function.
   
   **Return value:**  
   The return value is a string containing the name of
   the corresponding constant. 




.. _CInterface.ARKStepMainOutputs:

ARKStep timestepper module optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ============================================
Optional output                                      Function name
===================================================  ============================================
No. of explicit stability-limited steps              :c:func:`ARKStepGetNumExpSteps()`
No. of accuracy-limited steps                        :c:func:`ARKStepGetNumAccSteps()`
No. of attempted steps                               :c:func:`ARKStepGetNumStepAttempts()`
No. of calls to *fe* and *fi* functions              :c:func:`ARKStepGetNumRhsEvals()`
No. of local error test failures that have occurred  :c:func:`ARKStepGetNumErrTestFails()`
Current ERK and DIRK Butcher tables                  :c:func:`ARKStepGetCurrentButcherTables()`
Estimated local truncation error vector              :c:func:`ARKStepGetEstLocalErrors()`
Single accessor to many statistics at once           :c:func:`ARKStepGetTimestepperStats()`
===================================================  ============================================ 




.. c:function:: int ARKStepGetNumExpSteps(void* arkode_mem, long int* expsteps)

   Returns the cumulative number of stability-limited steps
   taken by the solver (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKStepGetNumAccSteps(void* arkode_mem, long int* accsteps)

   Returns the cumulative number of accuracy-limited steps
   taken by the solver (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *step_attempts* -- number of steps attempted by solver.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKStepGetNumRhsEvals(void* arkode_mem, long int* nfe_evals, long int* nfi_evals)

   Returns the number of calls to the user's right-hand
   side functions, :math:`f_E` and :math:`f_I` (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nfe_evals* -- number of calls to the user's :math:`f_E(t,y)` function.
      * *nfi_evals* -- number of calls to the user's :math:`f_I(t,y)` function.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The *nfi_evals* value does not account for calls made to
   :math:`f_I` by a linear solver or preconditioner module.



.. c:function:: int ARKStepGetNumErrTestFails(void* arkode_mem, long int* netfails)

   Returns the number of local error test failures that
   have occured (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *netfails* -- number of error test failures.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKStepGetCurrentButcherTables(void* arkode_mem, ARKodeButcherTable *Bi, ARKodeButcherTable *Be)

   Returns the explicit and implicit Butcher tables
   currently in use by the solver.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *Bi* -- pointer to implicit Butcher table structure.
      * *Be* -- pointer to explicit Butcher table structure.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:**  The *ARKodeButcherTable* data structure is defined in
   the header file ``arkode/arkode_butcher.h``.  It is defined as a
   pointer to the following C structure:

   .. code-block:: c

      typedef struct ARKodeButcherTableMem {

        int q;           /* method order of accuracy       */
        int p;           /* embedding order of accuracy    */
        int stages;      /* number of stages               */
        realtype **A;    /* Butcher table coefficients     */
        realtype *c;     /* canopy node coefficients       */
        realtype *b;     /* root node coefficients         */
        realtype *d;     /* embedding coefficients         */

      } *ARKodeButcherTable;


.. c:function:: int ARKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele)

   Returns the vector of estimated local truncation errors
   for the current step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *ele* -- vector of estimated local truncation errors.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:**  The user must allocate space for *ele*, that will be
   filled in by this function.
   
   The values returned in *ele* are valid only after a successful call
   to :c:func:`ARKode()` (i.e. it returned a non-negative value).
   
   The *ele* vector, together with the *eweight* vector from
   :c:func:`ARKodeGetErrWeights()`, can be used to determine how the
   various components of the system contributed to the estimated local
   error test.  Specifically, that error test uses the RMS norm of a
   vector whose components are the products of the components of these
   two vectors.  Thus, for example, if there were recent error test
   failures, the components causing the failures are those with largest
   values for the products, denoted loosely as ``eweight[i]*ele[i]``.



.. c:function:: int ARKStepGetTimestepperStats(void* arkode_mem, long int* expsteps, long int* accsteps, long int* step_attempts, long int* nfe_evals, long int* nfi_evals, long int* nlinsetups, long int* netfails)

   Returns many of the most useful timestepper statistics in a single call.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.
      * *step_attempts* -- number of steps attempted by the solver.
      * *nfe_evals* -- number of calls to the user's :math:`f_E(t,y)` function.
      * *nfi_evals* -- number of calls to the user's :math:`f_I(t,y)` function.
      * *nlinsetups* -- number of linear solver setup calls made.
      * *netfails* -- number of error test failures.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``

        

        
.. _CInterface.ERKStepMainOutputs:

ERKStep timestepper module optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ============================================
Optional output                                      Function name
===================================================  ============================================
No. of explicit stability-limited steps              :c:func:`ERKStepGetNumExpSteps()`
No. of accuracy-limited steps                        :c:func:`ERKStepGetNumAccSteps()`
No. of attempted steps                               :c:func:`ERKStepGetNumStepAttempts()`
No. of calls to *fe* and *fi* functions              :c:func:`ERKStepGetNumRhsEvals()`
No. of local error test failures that have occurred  :c:func:`ERKStepGetNumErrTestFails()`
Current ERK and DIRK Butcher tables                  :c:func:`ERKStepGetCurrentButcherTable()`
Estimated local truncation error vector              :c:func:`ERKStepGetEstLocalErrors()`
Single accessor to many statistics at once           :c:func:`ERKStepGetTimestepperStats()`
===================================================  ============================================ 




.. c:function:: int ERKStepGetNumExpSteps(void* arkode_mem, long int* expsteps)

   Returns the cumulative number of stability-limited steps
   taken by the solver (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ERKStepGetNumAccSteps(void* arkode_mem, long int* accsteps)

   Returns the cumulative number of accuracy-limited steps
   taken by the solver (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ERKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *step_attempts* -- number of steps attempted by solver.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ERKStepGetNumRhsEvals(void* arkode_mem, long int* nf_evals)

   Returns the number of calls to the user's right-hand
   side functions, :math:`f_E` and :math:`f_I` (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nf_evals* -- number of calls to the user's :math:`f(t,y)` function.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   


.. c:function:: int ERKStepGetNumErrTestFails(void* arkode_mem, long int* netfails)

   Returns the number of local error test failures that
   have occured (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *netfails* -- number of error test failures.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ERKStepGetCurrentButcherTable(void* arkode_mem, ARKodeButcherTable *B)

   Returns the Butcher table currently in use by the solver.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *B* -- pointer to Butcher table structure.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:**  The *ARKodeButcherTable* data structure is defined in
   the header file ``arkode/arkode_butcher.h``.  It is defined as a
   pointer to the following C structure:

   .. code-block:: c

      typedef struct ARKodeButcherTableMem {

        int q;           /* method order of accuracy       */
        int p;           /* embedding order of accuracy    */
        int stages;      /* number of stages               */
        realtype **A;    /* Butcher table coefficients     */
        realtype *c;     /* canopy node coefficients       */
        realtype *b;     /* root node coefficients         */
        realtype *d;     /* embedding coefficients         */

      } *ARKodeButcherTable;


.. c:function:: int ERKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele)

   Returns the vector of estimated local truncation errors
   for the current step.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *ele* -- vector of estimated local truncation errors.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:**  The user must allocate space for *ele*, that will be
   filled in by this function.
   
   The values returned in *ele* are valid only after a successful call
   to :c:func:`ARKode()` (i.e. it returned a non-negative value).
   
   The *ele* vector, together with the *eweight* vector from
   :c:func:`ARKodeGetErrWeights()`, can be used to determine how the
   various components of the system contributed to the estimated local
   error test.  Specifically, that error test uses the RMS norm of a
   vector whose components are the products of the components of these
   two vectors.  Thus, for example, if there were recent error test
   failures, the components causing the failures are those with largest
   values for the products, denoted loosely as ``eweight[i]*ele[i]``.



.. c:function:: int ERKStepGetTimestepperStats(void* arkode_mem, long int* expsteps, long int* accsteps, long int* step_attempts, long int* nf_evals, long int* netfails)

   Returns many of the most useful timestepper statistics in a single call.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.
      * *step_attempts* -- number of steps attempted by the solver.
      * *nf_evals* -- number of calls to the user's :math:`f(t,y)` function.
      * *netfails* -- number of error test failures.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``





                
.. _CInterface.ARKodeImplicitSolverOutputs:

Implicit solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ============================================
Optional output                                      Function name
===================================================  ============================================
No. of calls to linear solver setup function         :c:func:`ARKStepGetNumLinSolvSetups()`
No. of nonlinear solver iterations                   :c:func:`ARKStepGetNumNonlinSolvIters()`
No. of nonlinear solver convergence failures         :c:func:`ARKStepGetNumNonlinSolvConvFails()`
Single accessor to all nonlinear solver statistics   :c:func:`ARKStepGetNonlinSolvStats()`
===================================================  ============================================ 




.. c:function:: int ARKStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups)

   Returns the number of calls made to the linear solver's
   setup routine (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nlinsetups* -- number of linear solver setup calls made.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters)

   Returns the number of nonlinear solver iterations
   performed (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nniters* -- number of nonlinear iterations performed.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nncfails)

   Returns the number of nonlinear solver convergence
   failures that have occurred (so far).
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nncfails* -- number of nonlinear convergence failures.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``



.. c:function:: int ARKStepGetNonlinSolvStats(void* arkode_mem, long int* nniters, long int* nncfails)

   Returns all of the nonlinear solver statistics in a single call.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nniters* -- number of nonlinear iterations performed.
      * *nncfails* -- number of nonlinear convergence failures.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``




.. _CInterface.ARKodeRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===================================================  ==========================================
Optional output                                      Function name
===================================================  ==========================================
Array showing roots found                            :c:func:`ARKodeGetRootInfo()`
No. of calls to user root function                   :c:func:`ARKodeGetNumGEvals()`
===================================================  ========================================== 



.. c:function:: int ARKodeGetRootInfo(void* arkode_mem, int* rootsfound)

   Returns an array showing which functions were found to
   have a root.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *rootsfound* -- array of length *nrtfn* with the indices of the
        user functions :math:`g_i` found to have a root.  For 
	:math:`i = 0 \ldots` *nrtfn*-1, ``rootsfound[i]`` is nonzero
        if :math:`g_i` has a root, and 0 if not. 
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The user must allocate space for *rootsfound* prior to
   calling this function. 
   
   For the components of :math:`g_i` for which a root was found, the
   sign of ``rootsfound[i]`` indicates the direction of
   zero-crossing.  A value of +1 indicates that :math:`g_i` is
   increasing, while a value of -1 indicates a decreasing :math:`g_i`.



.. c:function:: int ARKodeGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *ngevals* -- number of calls made to :math:`g` so far.
   
   **Return value:**  
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKode memory was ``NULL``




.. _CInterface.ARKDlsOutputs:

Direct linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the ARKDLS
modules: workspace requirements, number of calls to the Jacobian
routine, number of calls to the mass matrix routine, number of calls
to the implicit right-hand side routine for finite-difference Jacobian
approximation, and last return value from an ARKDLS function.  Note
that, where the name of an output would otherwise conflict with the
name of an optional output from the main solver, a suffix LS (for
Linear Solver) or MLS (for Mass Linear Solver) has been added here
(e.g. *lenrwLS*).  


.. cssclass:: table-bordered

===================================================  ===================================
Optional output                                      Function name
===================================================  ===================================
Size of real and integer workspaces                  :c:func:`ARKDlsGetWorkSpace()`
Size of mass real and integer workspaces             :c:func:`ARKDlsGetMassWorkSpace()`
No. of Jacobian evaluations                          :c:func:`ARKDlsGetNumJacEvals()`
No. of mass matrix setups                            :c:func:`ARKDlsGetNumMassSetups()`
No. of mass matrix solves                            :c:func:`ARKDlsGetNumMassSolves()`
No. of mass matrix multiplies                        :c:func:`ARKDlsGetNumMassMult()`
No. of `fi` calls for finite diff. Jacobian evals    :c:func:`ARKDlsGetNumRhsEvals()`
Last return flag from a linear solver function       :c:func:`ARKDlsGetLastFlag()`
Last return flag from a mass matrix solver function  :c:func:`ARKDlsGetLastMassFlag()`
Name of constant associated with a return flag       :c:func:`ARKDlsGetReturnFlagName()`
===================================================  =================================== 



    
.. c:function:: int ARKDlsGetWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)

   Returns the real and integer workspace used by the
   ARKDLS linear solver interface.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lenrwLS* -- the number of ``realtype`` values in the ARKDLS workspace.
      * *leniwLS* -- the number of integer values in the ARKDLS workspace.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** The workspace requirements reported by this routine
   correspond only to memory allocated within this interface and to
   memory allocated by the ``SUNLinearSolver`` object attached
   to it.  The template Jacobian matrix allocated by the user outside
   of ARKDLS is not included in this report.



.. c:function:: int ARKDlsGetMassWorkSpace(void* arkode_mem, long int* lenrwMLS, long int* leniwMLS)

   Returns the real and integer workspace used by the
   ARKDLS mass matrix linear solver interface.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lenrwMLS* -- the number of ``realtype`` values in the ARKDLS workspace.
      * *leniwMLS* -- the number of integer values in the ARKDLS workspace.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** The workspace requirements reported by this routine
   correspond only to memory allocated within this interface and to
   memory allocated by the ``SUNLinearSolver`` object attached
   to it.  The template mass matrix allocated by the user outside
   of ARKDLS is not included in this report.


.. c:function:: int ARKDlsGetNumJacEvals(void* arkode_mem, long int* njevals)

   Returns the number of calls made to the ARKDLS Jacobian
   approximation routine. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *njevals* -- number of calls to the Jacobian function.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKDlsGetNumMassSetups(void* arkode_mem, long int* nmsetups)

   Returns the number of calls made to the ARKDLS mass matrix solver 'setup' routine. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmsetups* -- number of calls to the mass matrix solver setup routine.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKDlsGetNumMassSolves(void* arkode_mem, long int* nmsolves)

   Returns the number of calls made to the ARKDLS mass matrix solver 'solve' routine. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmsolvess* -- number of calls to the mass matrix solver solve routine.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKDlsGetNumMassMult(void* arkode_mem, long int* nmmults)

   Returns the number of calls made to the ARKDLS mass matrix 'matvec' routine.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmmults* -- number of calls to the mass matrix solver matrix-times-vector routine.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKDlsGetNumRhsEvals(void* arkode_mem, long int* nfevalsLS)

   Returns the number of calls made to the user-supplied
   :math:`f_I` routine due to the finite difference Jacobian
   approximation. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nfevalsLS* -- the number of calls made to the user-supplied
        :math:`f_I` function.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** The value of *nfevalsLS* is incremented only if one of
   the default internal difference quotient functions (dense or
   banded) is used.



.. c:function:: int ARKDlsGetLastFlag(void* arkode_mem, long int* lsflag)

   Returns the last return value from an ARKDLS routine.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lsflag* -- the value of the last return flag from an ARKDLS function.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:**   If the ``SUNLINSOL_DENSE`` or ``SUNLINSOL_BAND`` setup
   function failed (ARKode returned ``ARK_LSETUP_FAIL``), then the
   value of *lsflag* is equal to the column index (numbered from one)
   at which a zero diagonal element was encountered during the LU
   factorization of the (dense or banded) Jacobian matrix.  For all
   other failures, *lsflag* is negative.



.. c:function:: int ARKDlsGetLastMassFlag(void* arkode_mem, long int* mlsflag)

   Returns the last return value from an ARKDLS mass matrix solve routine.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *mlsflag* -- the value of the last return flag from an ARKDLS
	mass matrix solver function.
   
   **Return value:**  
      * *ARKDLS_SUCCESS* if successful
      * *ARKDLS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKDLS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** If the ``SUNLINSOL_DENSE`` or ``SUNLINSOL_BAND`` setup
   function failed (ARKode returned ``ARK_LSETUP_FAIL``), then the
   value of *lsflag* is equal to the column index (numbered from one)
   at which a zero diagonal element was encountered during the LU
   factorization of the (dense or banded) mass matrix.  For all
   other failures, *lsflag* is negative.



.. c:function:: char *ARKDlsGetReturnFlagName(long int lsflag)

   Returns the name of the ARKDLS constant
   corresponding to *lsflag*.
   
   **Arguments:**
      * *lsflag* -- a return flag from an ARKDLS function.
   
   **Return value:**  The return value is a string containing the name of
   the corresponding constant. If 1 :math:`\le` `lsflag` :math:`\le
   n` (LU factorization failed), this routine returns "NONE". 





.. _CInterface.ARKSpilsOutputs:

Iterative linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the ARKSPILS
modules: workspace requirements, number of linear iterations, number
of linear convergence failures, number of calls to the preconditioner
setup and solve routines, number of calls to the Jacobian-vector
setup and product routines, number of calls to the mass-matrix-vector
setup and product routines, number of calls to the implicit right-hand
side routine for finite-difference Jacobian-vector product
approximation, and last return value from a linear solver function.
Note that, where the name of an output would otherwise conflict with
the name of an optional output from the main solver, a suffix LS (for
Linear Solver) or MLS (for Mass Linear Solver) has been added here
(e.g. *lenrwLS*).


.. cssclass:: table-bordered

===========================================================  ======================================== 
Optional output                                              Function name
===========================================================  ========================================
Size of real and integer workspaces                          :c:func:`ARKSpilsGetWorkSpace()`
No. of preconditioner evaluations                            :c:func:`ARKSpilsGetNumPrecEvals()`
No. of preconditioner solves                                 :c:func:`ARKSpilsGetNumPrecSolves()`
No. of linear iterations                                     :c:func:`ARKSpilsGetNumLinIters()`
No. of linear convergence failures                           :c:func:`ARKSpilsGetNumConvFails()`
No. of Jacobian-vector setup evaluations                     :c:func:`ARKSpilsGetNumJTSetupEvals()`
No. of Jacobian-vector product evaluations                   :c:func:`ARKSpilsGetNumJtimesEvals()`
No. of *fi* calls for finite diff. Jacobian-vector evals.    :c:func:`ARKSpilsGetNumRhsEvals()`
Last return from a linear solver function                    :c:func:`ARKSpilsGetLastFlag()`
Size of real and integer mass matrix solver workspaces       :c:func:`ARKSpilsGetMassWorkSpace()`
No. of mass matrix preconditioner evaluations                :c:func:`ARKSpilsGetNumMassPrecEvals()`
No. of mass matrix preconditioner solves                     :c:func:`ARKSpilsGetNumMassPrecSolves()`
No. of mass matrix linear iterations                         :c:func:`ARKSpilsGetNumMassIters()`
No. of mass matrix solver convergence failures               :c:func:`ARKSpilsGetNumMassConvFails()`
No. of mass-matrix-vector setup evaluations                  :c:func:`ARKSpilsGetNumMTSetupEvals()`
No. of mass-matrix-vector product evaluations                :c:func:`ARKSpilsGetNumMtimesEvals()`
Last return from a mass matrix solver function               :c:func:`ARKSpilsGetLastMassFlag()`
Name of constant associated with a return flag               :c:func:`ARKSpilsGetReturnFlagName()`
===========================================================  ========================================




.. c:function:: int ARKSpilsGetWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)

   Returns the global sizes of the ARKSPILS real and integer workspaces.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lenrwLS* -- the number of ``realtype`` values in the ARKSPILS workspace.
      * *leniwLS* -- the number of integer values in the ARKSPILS workspace.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** The workspace requirements reported by this routine
   correspond only to memory allocated within this interface and to
   memory allocated by the ``SUNLinearSolver`` object attached
   to it.
   
   In a parallel setting, the above values are global (i.e. summed over all
   processors).



.. c:function:: int ARKSpilsGetNumPrecEvals(void* arkode_mem, long int* npevals)

   Returns the total number of preconditioner evaluations,
   i.e. the number of calls made to *psetup* with *jok* = ``SUNFALSE``.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *npevals* -- the current number of calls to *psetup*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumPrecSolves(void* arkode_mem, long int* npsolves)

   Returns the number of calls made to the preconditioner
   solve function, *psolve*.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *npsolves* -- the number of calls to *psolve*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumLinIters(void* arkode_mem, long int* nliters)

   Returns the cumulative number of linear iterations.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nliters* -- the current number of linear iterations.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumConvFails(void* arkode_mem, long int* nlcfails)

   Returns the cumulative number of linear convergence failures.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nlcfails* -- the current number of linear convergence failures.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumJTSetupEvals(void* arkode_mem, long int* njtsetup)

   Returns the cumulative number of calls made to the
   Jacobian-vector setup function, *jtsetup*.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *njtsetup* -- the current number of calls to *jtsetup*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumJtimesEvals(void* arkode_mem, long int* njvevals)

   Returns the cumulative number of calls made to the
   Jacobian-vector product function, *jtimes*.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *njvevals* -- the current number of calls to *jtimes*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumRhsEvals(void* arkode_mem, long int* nfevalsLS)

   Returns the number of calls to the user-supplied implicit
   right-hand side function :math:`f_I` for finite difference
   Jacobian-vector product approximation. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nfevalsLS* -- the number of calls to the user implicit
        right-hand side function.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** The value *nfevalsLS* is incremented only if the default
   internal difference quotient function is used.



.. c:function:: int ARKSpilsGetLastFlag(void* arkode_mem, long int* lsflag)

   Returns the last return value from an ARKSPILS routine.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lsflag* -- the value of the last return flag from an
        ARKSPILS function.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** If the ARKSPILS setup function failed (:c:func:`ARKode()`
   returned *ARK_LSETUP_FAIL*), then *lsflag* will be
   *SUNLS_PSET_FAIL_UNREC*, *SUNLS_ASET_FAIL_UNREC* or
   *SUNLS_PACKAGE_FAIL_UNREC*.
   
   If the ARKSPILS solve function failed (:c:func:`ARKode()`
   returned *ARK_LSOLVE_FAIL*), then *lsflag* contains the error
   return flag from the ``SUNLinearSolver`` object, which will
   be one of:
   *SUNLS_MEM_NULL*, indicating that the ``SUNLinearSolver``
   memory is ``NULL``;
   *SUNLS_ATIMES_FAIL_UNREC*, indicating an unrecoverable failure in
   the :math:`Jv` function;
   *SUNLS_PSOLVE_FAIL_UNREC*, indicating that the preconditioner solve
   function failed unrecoverably;
   *SUNLS_GS_FAIL*, indicating a failure in the Gram-Schmidt procedure
   (SPGMR and SPFGMR only);  
   *SUNLS_QRSOL_FAIL*, indicating that the matrix :math:`R` was found
   to be singular during the QR solve phase (SPGMR and SPFGMR only); or
   *SUNLS_PACKAGE_FAIL_UNREC*, indicating an unrecoverable failure in
   an external iterative linear solver package. 


   
.. c:function:: char *ARKSpilsGetReturnFlagName(long int lsflag)

   Returns the name of the ARKSPILS constant
   corresponding to *lsflag*.
   
   **Arguments:**
      * *lsflag* -- a return flag from an ARKSPILS function.
   
   **Return value:**  
   The return value is a string containing the name of
   the corresponding constant.
   



.. c:function:: int ARKSpilsGetMassWorkSpace(void* arkode_mem, long int* lenrwMLS, long int* leniwMLS)

   Returns the global sizes of the ARKSPILS real and integer workspaces.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *lenrwMLS* -- the number of ``realtype`` values in the ARKSPILS workspace.
      * *leniwMLS* -- the number of integer values in the ARKSPILS workspace.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** The workspace requirements reported by this routine
   correspond only to memory allocated within this interface and to
   memory allocated by the ``SUNLinearSolver`` object attached
   to it.
   
   In a parallel setting, the above values are global (i.e. summed over all
   processors).



.. c:function:: int ARKSpilsGetNumMassPrecEvals(void* arkode_mem, long int* nmpevals)

   Returns the total number of mass matrix preconditioner evaluations,
   i.e. the number of calls made to *psetup*.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmpevals* -- the current number of calls to *psetup*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumMassPrecSolves(void* arkode_mem, long int* nmpsolves)

   Returns the number of calls made to the mass matrix preconditioner
   solve function, *psolve*.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmpsolves* -- the number of calls to *psolve*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumMassIters(void* arkode_mem, long int* nmiters)

   Returns the cumulative number of mass matrix solver iterations.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmiters* -- the current number of mass matrix solver linear iterations.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumMassConvFails(void* arkode_mem, long int* nmcfails)

   Returns the cumulative number of mass matrix solver convergence failures.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmcfails* -- the current number of mass matrix solver convergence failures.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumMTSetupEvals(void* arkode_mem, long int* nmtsetup)

   Returns the cumulative number of calls made to the
   mass-matrix-vector setup function, *mtsetup*.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmtsetup* -- the current number of calls to *mtsetup*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumMtimesEvals(void* arkode_mem, long int* nmvevals)

   Returns the cumulative number of calls made to the
   mass-matrix-vector product function, *mtimes*.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *nmvevals* -- the current number of calls to *mtimes*.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetLastMassFlag(void* arkode_mem, long int* msflag)

   Returns the last return value from an ARKSPILS mass matrix solver routine.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *msflag* -- the value of the last return flag from an
        ARKSPILS mass matrix solver function.
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
      * *ARKSPILS_LMEM_NULL* if the linear solver memory was ``NULL``
   
   **Notes:** The values of *msflag* for each of the various solvers
   will match those described above for the function
   :c:func:`ARKSpilsGetLastFlag()`.   




.. _CInterface.ARKodeExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional routines may be called by a user to inquire
about existing solver parameters, to retrieve stored Butcher tables,
write the current Butcher table(s), or even to test a provided Butcher
table to determine its analytical order of accuracy.  While none of
these would typically be called during the course of solving an
initial value problem, these may be useful for users wishing to better
understand ARKode and/or specific Runge-Kutta methods.


.. cssclass:: table-bordered

===========================================================  ======================================== 
Optional routine                                             Function name
===========================================================  ========================================
Output all ARKode solver parameters                          :c:func:`ARKodeWriteParameters()`
Output all ARKStep timestepper parameters                    :c:func:`ARKStepWriteParameters()`
Output all ERKStep timestepper parameters                    :c:func:`ERKStepWriteParameters()`
Retrieve a given Butcher table by its unique name            :c:func:`ARKodeLoadButcherTable()`
Output the current Butcher table(s)                          :c:func:`ARKStepWriteButcher()`
                                                             :c:func:`ERKStepWriteButcher()`
===========================================================  ========================================




.. c:function:: int ARKodeWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all ARKode solver parameters to the provided file pointer.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *fp* -- pointer to use for printing the solver parameters
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The *fp* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.
   
   When run in parallel, only one process should set a non-NULL value
   for this pointer, since parameters for all processes would be
   identical.


.. c:function:: int ARKStepWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all ARKStep timestepper parameters to the provided file pointer.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *fp* -- pointer to use for printing the solver parameters
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The *fp* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.
   
   When run in parallel, only one process should set a non-NULL value
   for this pointer, since parameters for all processes would be
   identical.


.. c:function:: int ERKStepWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all ERKStep timestepper parameters to the provided file pointer.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *fp* -- pointer to use for printing the solver parameters
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The *fp* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.
   
   When run in parallel, only one process should set a non-NULL value
   for this pointer, since parameters for all processes would be
   identical.


.. c:function:: ARKodeButcherTable ARKodeLoadButcherTable(int imethod)

   Retrieves a specified Butcher table.  The *ARKodeButcherTable* data
   structure is defined in the header file
   ``arkode/arkode_butcher.h``, and is described above in the notes
   for the function :c:func:`ARKodeGetCurrentButcherTables()`.
   
   **Arguments:**
      * *imethod* -- integer input specifying the given Butcher table --
        valid values match those for the functions
        :c:func:`ARKodeSetERKTableNum()` and :c:func:`ARKodeSetIRKTableNum()`
   
   **Return value:**  
      * *ARKodeButcherTable* structure if successful
      * *NULL* pointer if *imethod* was invalid
   

.. c:function:: int ARKStepWriteButcher(void* arkode_mem, FILE *fp)

   Outputs the current Butcher table(s) to the provided file pointer.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *fp* -- pointer to use for printing the Butcher table(s)
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The *fp* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.

   If ARKStep is currently configured to run in purely explicit or
   purely implicit mode, this will output a single Butcher table;  if
   configured to run an ImEx method then both tables will be output.
   
   When run in parallel, only one process should set a non-NULL value
   for this pointer, since tables for all processes would be
   identical.
   

.. c:function:: int ERKStepWriteButcher(void* arkode_mem, FILE *fp)

   Outputs the current Butcher table to the provided file pointer.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *fp* -- pointer to use for printing the Butcher table
   
   **Return value:**  
      * *ARKSPILS_SUCCESS* if successful
      * *ARKSPILS_MEM_NULL* if the ARKode memory was ``NULL``
   
   **Notes:** The *fp* argument can be ``stdout`` or ``stderr``, or it
   may point to a specific file created using ``fopen``.

   When run in parallel, only one process should set a non-NULL value
   for this pointer, since tables for all processes would be
   identical.






.. _CInterface.Reinitialization:

ARKode reinitialization functions
-------------------------------------

To reinitialize the timestepper module and main ARKode solver for the
solution of a new problem, where a prior call to
:c:func:`ARKStepCreate()` has been made, the user must call the
function :c:func:`ARKStepReInit()`.  Similarly, if using the ERK
timestepper module, then the user must call the
function :c:func:`ERKStepReInit()` to reinitialize ERKStep and ARKode
for a new problem.  In both cases, the new problem must have the same
size as the previous one.  These routines perform the same input
checking and initializations that are done in
:c:func:`ARKStepCreate()` and :c:func:`ERKStepCreate()`, but they
perform no memory allocation as they assume that the existing internal
memory is sufficient for the new problem.  A call to either of these
re-initialization routines deletes the solution history that was
stored internally during the previous integration.  Following a
successful call to :c:func:`ARKStepReInit()` or
:c:func:`ERKStepReInit()`, call :c:func:`ARKode()` again for the
solution of the new problem. 

The use of :c:func:`ARKStepReInit()` or :c:func:`ERKStepReInit()`
requires that the number of Runge Kutta stages, denoted by *s*, be no
larger for the new problem than for the previous problem.  This
condition is automatically fulfilled if the method order *q* and the
problem type (explicit, implicit, ImEx) are left unchanged.

When using the ARKStep timestepping module, if there are changes to
the linear solver specifications, the user should make the appropriate
calls to either the linear solver objects themselves, or to the ARKDLS
or ARKSPILS interface routines, as described  in the section 
:ref:`CInterface.LinearSolvers`. Otherwise, all solver inputs set
previously remain in effect. 

One important use of the :c:func:`ARKStepReInit()` and
:c:func:`ERKStepReInit()` functions is in the treating of jump
discontinuities in the RHS function.  Except in cases of fairly small
jumps, it is usually more efficient to stop at each point of
discontinuity and restart the integrator with a readjusted ODE model,
using a call to one of these two routines.  To stop when the location
of the discontinuity is known, simply make that location a value of
``tout``.  To stop when the location of the discontinuity is
determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS function *not* incorporate the
discontinuity, but rather have a smooth extention over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
function (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int ARKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, realtype t0, N_Vector y0)

   Provides required problem specifications and reinitializes the
   ARKStep timestepper module.
   
   ..
      **Arguments:**
         * *arkode_mem* -- pointer to the ARKode memory block.
         * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
           defining the explicit portion of the right-hand side function in 
           :math:`M(t)\, \dot{y} = f_E(t,y) + f_I(t,y)`.
         * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
           defining the implicit portion of the right-hand side function in 
           :math:`M(t)\, \dot{y} = f_E(t,y) + f_I(t,y)`.
         * *t0* -- the initial value of :math:`t`.
         * *y0* -- the initial condition vector :math:`y(t_0)`.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *fe* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit portion of the right-hand side function in 
        :math:`M\, \dot{y} = f_E(t,y) + f_I(t,y)`.
      * *fi* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit portion of the right-hand side function in 
        :math:`M\, \dot{y} = f_E(t,y) + f_I(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument has an illegal value.
   
   **Notes:** If an error occurred, :c:func:`ARKStepReInit()` also
   sends an error message to the error handler function.


.. c:function:: int ERKStepReInit(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0)

   Provides required problem specifications and reinitializes the
   ERKStep timestepper module. 
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *f* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function in :math:`\dot{y} = f(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument has an illegal value.
   
   **Notes:** If an error occurred, :c:func:`ERKStepReInit()` also
   sends an error message to the error handler function.




.. _CInterface.Resizing:

ARKode system resize function
-------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatially-adaptive
PDE simulations under a method-of-lines approach), the ARKode
integrator may be "resized" between integration steps, through calls
to the :c:func:`ARKodeResize()` function. This function modifies
ARKode's internal memory structures to use the new problem size,
without destruction of the temporal adaptivity heuristics.  It is
assumed that the dynamical time scales before and after the vector
resize will be comparable, so that all time-stepping heuristics prior
to calling :c:func:`ARKodeResize()` remain valid after the call.  If
instead the dynamics should be recomputed from scratch, the ARKode
memory structure should be deleted with a call to
:c:func:`ARKodeFree()`, and recreated with calls to
:c:func:`ARKodeCreate()` and either :c:func:`ARKStepInit()` or
:c:func:`ERKStepInit()`.

To aid in the vector resize operation, the user can supply a vector
resize function that will take as input a vector with the previous
size, and transform it in-place to return a corresponding vector of
the new size.  If this function (of type :c:func:`ARKVecResizeFn()`)
is not supplied (i.e. is set to ``NULL``), then all existing vectors
internal to ARKode will be destroyed and re-cloned from the new input
vector. 

In the case that the dynamical time scale should be modified slightly
from the previous time scale, an input *hscale* is allowed, that will
rescale the upcoming time step by the specified factor.  If a value
*hscale* :math:`\le 0` is specified, the default of 1.0 will be used.



.. c:function:: int ARKodeResize(void* arkode_mem, N_Vector ynew, realtype hscale, realtype t0, ARKVecResizeFn resize, void* resize_data)

   Re-initializes ARKode with a different state vector but with
   comparable dynamical time scale.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *ynew* -- the newly-sized solution vector, holding the current
	dependent variable values :math:`y(t_0)`. 
      * *hscale* -- the desired scaling factor for the dynamical time
	scale (i.e. the next step will be of size *h\*hscale*).
      * *t0* -- the current value of the independent variable
	:math:`t_0` (this must be consistent with *ynew*.
      * *resize* -- the user-supplied vector resize function (of type
	:c:func:`ARKVecResizeFn()`.
      * *resize_data* -- the user-supplied data structure to be passed
	to *resize* when modifying internal ARKode vectors.
   
   **Return value:** 
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ARKode memory was ``NULL``
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if an argument has an illegal value.
   
   **Notes:** If an error occurred, :c:func:`ARKodeResize()` also sends an error
   message to the error handler function.



Resizing the linear solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using any of the SUNDIALS-provided linear solver modules, the
linear solver memory structures must also be resized.  At present,
none of these include a solver-specific 'resize' function, so the linear
solver memory must be destroyed and re-allocated **following** each
call to :c:func:`ARKodeResize()`.  Moreover, the existing ARKDLS or
ARKSPILS interface should then be deleted and recreated by attaching
the updated ``SUNLinearSolver`` (and possibly
``SUNMatrix``) object(s) through calls to
:c:func:`ARKDlsSetLinearSolver()`, 
:c:func:`ARKSpilsSetLinearSolver()`,
:c:func:`ARKDlsSetMassLinearSolver()` and
:c:func:`ARKSpilsSetMassLinearSolver()`. 

If any user-supplied routines are provided to aid the linear solver
(e.g. Jacobian construction, Jacobian-vector product,
mass-matrix-vector product, preconditioning), then the corresponding
"set" routines must be called again **following** the solver
re-specification.


Resizing the absolute tolerance array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If using array-valued absolute tolerances, the absolute tolerance
vector will be invalid after the call to :c:func:`ARKodeResize()`, so
the new absolute tolerance vector should be re-set **following** each
call to :c:func:`ARKodeResize()` through a new call to
:c:func:`ARKodeSVtolerances()` (and similarly to
:c:func:`ARKodeResVtolerance()` if that was used for the original
problem).

If scalar-valued tolerances or a tolerance function was specified
through either :c:func:`ARKodeSStolerances()` or
:c:func:`ARKodeWFtolerances()`, then these will remain valid. and no
further action is necessary. 


.. note:: For an example of :c:func:`ARKodeResize()` usage, see the
	  supplied serial C example problem, ``ark_heat1D_adapt.c``.
