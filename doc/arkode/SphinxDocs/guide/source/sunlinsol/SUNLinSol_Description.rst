..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _SUNLinSol.Description:

Description of the SUNLinearSolver Module
===========================================

For problems that involve the solution of linear systems of equations,
the SUNDIALS solvers operate using generic linear solver modules
(of type ``SUNLinearSolver``), through a set of operations defined by 
the particular ``SUNLinearSolver`` implementation.  These work in coordination
with the SUNDIALS generic ``N_Vector`` and ``SUNMatrix`` modules to
provide a set of compatible data structures and solvers for the
solution of linear systems using direct or iterative
methods. Moreover, users can provide their own specific
``SUNLinearSolver`` implementation to each SUNDIALS solver,
particularly in cases where they provide their own ``N_Vector`` and/or
``SUNMatrix`` modules, and the customized linear solver leverages
these additional data structures to create highly efficient and/or
scalable solvers for their particular problem.  Additionally, SUNDIALS
provides native implementations ``SUNLinearSolver`` modules, as well
as ``SUNLinearSolver`` modules that interface between SUNDIALS and
external linear solver libraries.

The various SUNDIALS solvers have been designed to specifically
leverage the use of either *direct linear solvers*
or matrix-free, *scaled, preconditioned, iterative linear solvers*,
through their "Dls" and "Spils" interfaces, respectively.
Additionally, SUNDIALS solvers can make use of user-supplied custom
linear solvers, whether these are problem-specific or come from
external solver libraries.

For iterative (and possibly custom) linear solvers, the SUNDIALS
solvers leverage scaling and preconditioning, as applicable, to
balance error between solution components and to accelerate
convergence of the linear solver.  To this end, instead of solving the 
linear system :math:`Ax = b` directly, we apply the underlying
iterative algorithm to the transformed system  

.. math::
   :label: eq:transformed_linear_system
   
   \tilde{A} \tilde{x} = \tilde{b}

where

.. math::
  :label: eq:transformed_linear_system_components
   
  \tilde{A} &= S_1 P_1^{-1} A P_2^{-1} S_2^{-1},\\
  \tilde{b} &= S_1 P_1^{-1} b,\\
  \tilde{x} &= S_2 P_2 x,

and where

* :math:`P_1` is the left preconditioner,
  
* :math:`P_2` is the right preconditioner,
    
* :math:`S_1` is a diagonal matrix of scale factors for
  :math:`P_1^{-1} b`,
      
* :math:`S_2` is a diagonal matrix of scale factors for :math:`P_2 x`.

The SUNDIALS solvers request that iterative linear solvers stop
based on the 2-norm of the scaled preconditioned residual meeting a
prescribed tolerance

.. math::
   
   \left\| \tilde{b} - \tilde{A} \tilde{x} \right\|_2  <  \text{tol}.

We note that not all of the iterative linear solvers implemented in
SUNDIALS support the full range of the above options.  Similarly,
some of the SUNDIALS integrators only utilize a subset of these
options.  Exceptions to the operators shown above are described in
the documentation for each ``SUNLinearSolver`` implementation, or for each
SUNDIALS solver "Spils" interface.

The generic ``SUNLinearSolver`` type has been modeled after the
object-oriented style of the generic ``N_Vector`` type.
Specifically, a generic ``SUNLinearSolver`` is a pointer to a structure
that has an implementation-dependent *content* field containing
the description and actual data of the linear solver, and an *ops*
field pointing to a structure with generic linear solver operations.
The type ``SUNLinearSolver`` is defined as

.. code-block:: c

   typedef struct _generic_SUNLinearSolver *SUNLinearSolver;

   struct _generic_SUNLinearSolver {
     void *content;
     struct _generic_SUNLinearSolver_Ops *ops;
   };

The ``_generic_SUNLinearSolver_Ops`` structure is essentially a
list of pointers to the various actual linear solver operations, and
is defined as 

.. code-block:: c

   struct _generic_SUNLinearSolver_Ops {
     SUNLinearSolver_Type (*gettype)(SUNLinearSolver);
     int                  (*setatimes)(SUNLinearSolver, void*, ATimesFn);
     int                  (*setpreconditioner)(SUNLinearSolver, void*, 
                                               PSetupFn, PSolveFn);
     int                  (*setscalingvectors)(SUNLinearSolver,
                                               N_Vector, N_Vector);
     int                  (*initialize)(SUNLinearSolver);
     int                  (*setup)(SUNLinearSolver, SUNMatrix);
     int                  (*solve)(SUNLinearSolver, SUNMatrix, N_Vector, 
                                   N_Vector, realtype);
     int                  (*numiters)(SUNLinearSolver);
     realtype             (*resnorm)(SUNLinearSolver);
     long int             (*lastflag)(SUNLinearSolver);
     int                  (*space)(SUNLinearSolver, long int*, long int*);
     N_Vector             (*resid)(SUNLinearSolver);
     int                  (*free)(SUNLinearSolver);
   };


The generic ``SUNLinearSolver`` module defines and implements the
linear solver operations acting on ``SUNLinearSolver`` objects.
These routines are in fact only wrappers for the linear solver operations
defined by a particular ``SUNLinearSolver`` implementation, which are
accessed through the {\em ops} field of the ``SUNLinearSolver``
structure. To illustrate this point we show below the implementation
of a typical linear solver operation from the generic ``SUNLinearSolver``
module, namely ``SUNLinSolInitialize``, which initializes a
``SUNLinearSolver`` object for use after it has been created and configured,
and returns a flag denoting a successful/failed operation:

.. code-block:: c

   int SUNLinSolInitialize(SUNLinearSolver S)
   {
     return ((int) S->ops->initialize(S));
   }

The subsection :ref:`SUNLinSol.Ops` contains a complete list of all
linear solver operations defined by the generic ``SUNLinearSolver``
module.  In order to support both direct and iterative linear solver
types, the generic ``SUNLinearSolver`` module defines linear solver
routines (or arguments) that may be specific to individual use cases.
As such, for each routine we specify its intended use.  If a custom 
``SUNLinearSolver`` module is provided, the function pointers for
non-required routines may be set to ``NULL`` to indicate that they
are not provided.

A particular implementation of the ``SUNLinearSolver`` module must:

* Specify the *content* field of the ``SUNLinearSolver`` object.
  
* Define and implement a minimal subset of the linear solver
  operations. See the documentation for each SUNDIALS linear solver
  interface to determine which ``SUNLinearSolver`` operations they
  require. 

  Note that the names of these routines should be unique to that
  implementation in order to permit using more than one
  ``SUNLinearSolver`` module (each with different ``SUNLinearSolver``
  internal data representations) in the same code.

* Define and implement user-callable constructor and destructor
  routines to create and free a ``SUNLinearSolver`` with
  the new *content* field and with *ops* pointing to the
  new linear solver operations.

* Optionally, define and implement additional user-callable routines
  acting on the newly defined ``SUNLinearSolver`` (e.g., routines to
  set various configuration options for tuning the linear solver to a
  particular problem).

* Optionally, provide functions as needed for that particular
  implementation to access different parts in the *content* field of
  the newly defined ``SUNLinearSolver`` object (e.g., routines to
  return various statistics from the solver). 

Each ``SUNLinearSolver`` implementation included in SUNDIALS has a "type"
identifier specified in enumeration and shown in Table :ref:`SUNLinSol.linsolIDs`.
It is recommended that a user-supplied ``SUNLinearSolver`` implementation set
this identifier based on the SUNDIALS solver interface they intend
to use: "Dls" interfaces require the ``SUNLINEARSOLVER_DIRECT``
``SUNLinearSolver`` objects, "Spils" interfaces require the
``SUNLINEARSOLVER_ITERATIVE`` objects. 


.. _SUNLinSol.linsolIDs:

Identifiers associated with linear solver kernels supplied with SUNDIALS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=========================  =================  ========
Linear Solver ID           Solver type        ID Value
=========================  =================  ========
SUNLINEARSOLVER_DIRECT     Direct solvers     0
SUNLINEARSOLVER_ITERATIVE  Iterative solvers  1
SUNLINEARSOLVER_CUSTOM     Custom solvers     2
=========================  =================  ========



