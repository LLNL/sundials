..
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2020, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

:tocdepth: 3

.. _SUNNonlinSol.ARKode:

====================================
ARKode SUNNonlinearSolver interface
====================================

As discussed in :ref:`Mathematics` integration steps often require the
(approximate) solution of a nonlinear system. This system can be formulated as
the rootfinding problem

.. math::
   G(z_i) \equiv M z_i - M y_n - \gamma f(t_i, z_i) - a_i = 0

or, when :math:`M=I`, as the fixed-point problem

.. math::
   z_i = y_n + \gamma f(t_i, z_i) + a_i

where :math:`z_i` is the i-th stage at time :math:`t_i` and :math:`a_i` is known
data that depends on the integration method.

Rather than solving the above nonlinear systems for the stage value :math:`z_i`
ARKode modules solve for the correction :math:`z_{cor}` to the predicted stage
value :math:`z_{pred}` so that :math:`z_i = z_{pred} + z_{cor}`. The nonlinear
systems rewritten in terms of :math:`z_{cor}` are

.. math::
   G(z_{cor}) \equiv M z_{cor} + M z_{pred} - M y_n - \gamma f(t_i, z_{cor} + z_{pred}) - a_i = 0

for the rootfinding problem and

.. math::
   z_{cor} = y_n - z_{pred} - G\left( z_{cor} + z{pred} \right) + a_i

for the fixed-point problem.

The nonlinear system functions provided by ARKode modules to the nonlinear
solver module internally update the current value of the stage based on the
input correction vector. The updated vector is used when calling the ODE
right-hand side function and when setting up linear solves (e.g., updating the
Jacobian or preconditioner).

ARKode modules also provide several advanced function that will not be needed by
most users, but might be useful for users who choose to provide their own
implementation of the SUNNonlinearSolver API. For example, such a user might
need access to the current value of :math:`\gamma` to comput Jacboian data.

ARKStep advanced output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int ARKStepGetCurrentState(void* arkode_mem, N_Vector* y)

   Returns the current state vector

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *y* -- N_Vector pointer that will get set to the current state vector

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory was ``NULL``

.. c:function:: int ARKStepGetCurrentGamma(void* arkode_mem, realtype* gamma)

   Returns the current value of the scalar :math:`\gamma`

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *gamma* -- the current value of the scalar :math:`\gamma` appearing in the
      Newton equation :math:`A = I - \gamma J` or :math:`A = M - \gamma J`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ARKStep memory was ``NULL``


MRIStep advanced output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int MRIStepGetCurrentState(void* arkode_mem, N_Vector* y)

   Returns the current state vector

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *y* -- N_Vector pointer that will get set to the current state vector

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the MRIStep memory was ``NULL``
