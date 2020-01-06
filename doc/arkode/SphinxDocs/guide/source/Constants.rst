..
   Programmer(s): Daniel R. Reynolds @ SMU
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

.. _Constants:

===========================
Appendix: ARKode Constants
===========================

Below we list all input and output constants used by the main solver,
timestepper, and linear solver modules, together with their numerical
values and a short description of their meaning.


ARKode input constants
==========================

Shared ARKode input constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`ARK_NORMAL` (1):
   Solver returns at a specified output time.

:index:`ARK_ONE_STEP`  (2):
   Solver returns after each successful step.


Explicit Butcher table specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`HEUN_EULER_2_1_2 <Heun-Euler-2-1-2 ERK method>`  (0):
   Use the Heun-Euler-2-1-2 ERK method

:index:`BOGACKI_SHAMPINE_4_2_3 <Bogacki-Shampine-4-2-3 ERK method>`  (1):
   Use the Bogacki-Shampine-4-2-3 ERK method

:index:`ARK324L2SA_ERK_4_2_3 <ARK-4-2-3 ERK method>`  (2):
   Use the ARK-4-2-3 ERK method

:index:`ZONNEVELD_5_3_4 <Zonneveld-5-3-4 ERK method>`  (3):
   Use the Zonneveld-5-3-4 ERK method

:index:`ARK436L2SA_ERK_6_3_4 <ARK-6-3-4 ERK method>`  (4):
   Use the ARK-6-3-4 ERK method

:index:`SAYFY_ABURUB_6_3_4 <Sayfy-Aburub-6-3-4 ERK method>`  (5):
   Use the Sayfy-Aburub-6-3-4 ERK method

:index:`CASH_KARP_6_4_5 <Cash-Karp-6-4-5 ERK method>`  (6):
   Use the Cash-Karp-6-4-5 ERK method

:index:`FEHLBERG_6_4_5 <Fehlberg-6-4-5 ERK method>`  (7):
   Use the Fehlberg-6-4-5 ERK method

:index:`DORMAND_PRINCE_7_4_5 <Dormand-Prince-7-4-5 ERK method>`  (8):
   Use the Dormand-Prince-7-4-5 ERK method

:index:`ARK548L2SA_ERK_8_4_5 <ARK-8-4-5 ERK method>`  (9):
   Use the ARK-8-4-5 ERK method

:index:`VERNER_8_5_6 <Verner-8-5-6 ERK method>`  (10):
   Use the Verner-8-5-6 ERK method

:index:`FEHLBERG_13_7_8 <Fehlberg-13-7-8 ERK method>`  (11):
   Use the Fehlberg-13-7-8 ERK method

:index:`KNOTH_WOLKE_3_3 <Knoth-Wolke-3-3 ERK method>`  (12):
   Use the Knoth-Wolke-3-3 ERK method


:index:`DEFAULT_ERK_2`  (HEUN_EULER_2_1_2):
   Use the default second-order ERK method

:index:`DEFAULT_ERK_3`  (BOGACKI_SHAMPINE_4_2_3):
   Use the default third-order ERK method

:index:`DEFAULT_ERK_4`  (ZONNEVELD_5_3_4):
   Use the default fourth-order ERK method

:index:`DEFAULT_ERK_5`  (CASH_KARP_6_4_5):
   Use the default fifth-order ERK method

:index:`DEFAULT_ERK_6`  (VERNER_8_5_6):
   Use the default sixth-order ERK method

:index:`DEFAULT_ERK_8`  (FEHLBERG_13_7_8):
   Use the default eighth-order ERK method




Implicit Butcher table specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`SDIRK_2_1_2 <SDIRK-2-1-2 method>`  (12):
   Use the SDIRK-2-1-2 SDIRK method

:index:`BILLINGTON_3_3_2 <Billington-3-3-2 SDIRK method>`  (13):
   Use the Billington-3-3-2 SDIRK method

:index:`TRBDF2_3_3_2 <TRBDF2-3-3-2 ESDIRK method>`  (14):
   Use the TRBDF2-3-3-2 ESDIRK method

:index:`KVAERNO_4_2_3 <Kvaerno-4-2-3 ESDIRK method>`  (15):
   Use the Kvaerno-4-2-3 ESDIRK method

:index:`ARK324L2SA_DIRK_4_2_3 <ARK-4-2-3 ESDIRK method>`  (16):
   Use the ARK-4-2-3 ESDIRK method

:index:`CASH_5_2_4 <Cash-5-2-4 SDIRK method>`  (17):
   Use the Cash-5-2-4 SDIRK method

:index:`CASH_5_3_4 <Cash-5-3-4 SDIRK method>`  (18):
   Use the Cash-5-3-4 SDIRK method

:index:`SDIRK_5_3_4 <SDIRK-5-3-4 method>`  (19):
   Use the SDIRK-5-3-4 SDIRK method

:index:`KVAERNO_5_3_4 <Kvaerno-5-3-4 ESDIRK method>`  (20):
   Use the Kvaerno-5-3-4 ESDIRK method

:index:`ARK436L2SA_DIRK_6_3_4 <ARK-6-3-4 ESDIRK method>`  (21):
   Use the ARK-6-3-4 ESDIRK method

:index:`KVAERNO_7_4_5 <Kvaerno-7-4-5 ESDIRK method>`  (22):
   Use the Kvaerno-7-4-5 ESDIRK method

:index:`ARK548L2SA_DIRK_8_4_5 <ARK-8-4-5 ESDIRK method>`  (23):
   Use the ARK-8-4-5 ESDIRK method


:index:`DEFAULT_DIRK_2`  (SDIRK_2_1_2):
   Use the default second-order DIRK method

:index:`DEFAULT_DIRK_3`  (ARK324L2SA_DIRK_4_2_3):
   Use the default third-order DIRK method

:index:`DEFAULT_DIRK_4`  (SDIRK_5_3_4):
   Use the default fourth-order DIRK method

:index:`DEFAULT_DIRK_5`  (ARK548L2SA_DIRK_8_4_5):
   Use the default fifth-order DIRK method



ImEx Butcher table specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`ARK324L2SA_ERK_4_2_3 and ARK324L2SA_DIRK_4_2_3 <ARK-4-2-3 ARK method>` (2 and 16):
   Use the ARK-4-2-3 ARK method

:index:`ARK436L2SA_ERK_6_3_4 and ARK436L2SA_DIRK_6_3_4 <ARK-6-3-4 ARK method>` (4 and 21):
   Use the ARK-6-3-4 ARK method

:index:`ARK548L2SA_ERK_8_4_5 and ARK548L2SA_DIRK_8_4_5 <ARK-8-4-5 ARK method>` (9 and 23):
   Use the ARK-8-4-5 ARK method


:index:`DEFAULT_ARK_ETABLE_3` and :index:`DEFAULT_ARK_ITABLE_3` (ARK324L2SA_[ERK,DIRK]_4_2_3):
   Use the default third-order ARK method

:index:`DEFAULT_ARK_ETABLE_4` and :index:`DEFAULT_ARK_ITABLE_4` (ARK436L2SA_[ERK,DIRK]_6_3_4):
   Use the default fourth-order ARK method

:index:`DEFAULT_ARK_ETABLE_5` and :index:`DEFAULT_ARK_ITABLE_5` (ARK548L2SA_[ERK,DIRK]_8_4_5):
   Use the default fifth-order ARK method




ARKode output constants
==========================

Shared ARKode output constants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`ARK_SUCCESS`  (0):
   Successful function return.

:index:`ARK_TSTOP_RETURN`  (1):
   ARKode succeeded by reaching the specified
   stopping point.

:index:`ARK_ROOT_RETURN`  (2):
   ARKode succeeded and found one more more roots.

:index:`ARK_WARNING`  (99):
   ARKode succeeded but an unusual situation occurred.

:index:`ARK_TOO_MUCH_WORK`  (-1):
   The solver took ``mxstep`` internal steps
   but could not reach ``tout``.

:index:`ARK_TOO_MUCH_ACC`  (-2):
   The solver could not satisfy the accuracy
   demanded by the user for some internal step.

:index:`ARK_ERR_FAILURE`  (-3):
   Error test failures occurred too many times
   during one internal time step, or the minimum step size was
   reached.

:index:`ARK_CONV_FAILURE`  (-4):
   Convergence test failures occurred too many
   times during one internal time step, or the minimum step size was
   reached.

:index:`ARK_LINIT_FAIL`  (-5):
   The linear solver's initialization function failed.

:index:`ARK_LSETUP_FAIL`  (-6):
   The linear solver's setup function failed in
   an unrecoverable manner.

:index:`ARK_LSOLVE_FAIL`  (-7):
   The linear solver's solve function failed in
   an unrecoverable manner.

:index:`ARK_RHSFUNC_FAIL`  (-8):
   The right-hand side function failed in an
   unrecoverable manner.

:index:`ARK_FIRST_RHSFUNC_ERR`  (-9):
   The right-hand side function failed
   at the first call.

:index:`ARK_REPTD_RHSFUNC_ERR`  (-10):
   The right-hand side function had
   repeated recoverable errors.

:index:`ARK_UNREC_RHSFUNC_ERR`  (-11):
   The right-hand side function had a
   recoverable error, but no recovery is possible.

:index:`ARK_RTFUNC_FAIL`  (-12):
   The rootfinding function failed in an
   unrecoverable manner.

:index:`ARK_LFREE_FAIL`  (-13):
   The linear solver's memory deallocation function failed.

:index:`ARK_MASSINIT_FAIL`  (-14):
   The mass matrix linear solver's initialization function failed.

:index:`ARK_MASSSETUP_FAIL`  (-15):
   The mass matrix linear solver's setup function failed in
   an unrecoverable manner.

:index:`ARK_MASSSOLVE_FAIL`  (-16):
   The mass matrix linear solver's solve function failed in
   an unrecoverable manner.

:index:`ARK_MASSFREE_FAIL`  (-17):
   The mass matrix linear solver's memory deallocation function failed.

:index:`ARK_MASSMULT_FAIL`  (-18):
   The mass matrix-vector product function failed.

:index:`ARK_CONSTR_FAIL`  (-19):
   The inequality constraint test failed repeatedly or failed with the minimum
   step size.

:index:`ARK_MEM_FAIL`  (-20):
   A memory allocation failed.

:index:`ARK_MEM_NULL`  (-21):
   The ``arkode_mem`` argument was ``NULL``.

:index:`ARK_ILL_INPUT`  (-22):
   One of the function inputs is illegal.

:index:`ARK_NO_MALLOC`  (-23):
   The ARKode memory block was not allocated by
   a call to :c:func:`ARKodeMalloc()`.

:index:`ARK_BAD_K`  (-24):
   The derivative order :math:`k` is larger than allowed.

:index:`ARK_BAD_T`  (-25):
   The time :math:`t` is outside the last step taken.

:index:`ARK_BAD_DKY`  (-26):
   The output derivative vector is ``NULL``.

:index:`ARK_TOO_CLOSE`  (-27):
   The output and initial times are too close to
   each other.

..
   :index:`ARK_POSTPROCESS_FAIL`  (-28):
      An error occurred when calling the user-provided ``ARKPostProcessStepFn`` routine.

:index:`ARK_VECTOROP_ERR`  (-29):
   An error occurred when calling an NVECTOR routine.

:index:`ARK_NLS_INIT_FAIL`  (-30):
   An error occurred when initializing a SUNNonlinearSolver module.

:index:`ARK_NLS_SETUP_FAIL`  (-31):
   A non-recoverable error occurred when setting up a
   SUNNonlinearSolver module.

:index:`ARK_NLS_SETUP_RECVR`  (-32):
   A recoverable error occurred when setting up a SUNNonlinearSolver module.

:index:`ARK_NLS_OP_ERR`  (-33):
   An error occurred when calling a set/get routine in a SUNNonlinearSolver
   module.

:index:`ARK_INNERSTEP_ATTACH_ERR`  (-34):
   An error occurred when attaching the inner stepper module.

:index:`ARK_INNERSTEP_FAIL`  (-35):
   An error occurred in the inner stepper module.

:index:`ARK_PREINNERFN_FAIL`  (-36):
   An error occurred in the MRIStep pre inner integrator function.

:index:`ARK_POSTINNERFN_FAIL`  (-37):
   An error occurred in the MRIStep post inner integrator function.

:index:`ARK_UNRECOGNIZED_ERROR` (-99):
   An unknown error was encountered.



ARKLS linear solver modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`ARKLS_SUCCESS`  (0):
   Successful function return.

:index:`ARKLS_MEM_NULL`  (-1):
   The ``arkode_mem`` argument was ``NULL``.

:index:`ARKLS_LMEM_NULL`  (-2):
   The ARKLS linear solver interface has not been initialized.

:index:`ARKLS_ILL_INPUT`  (-3):
   The ARKLS solver interface is not compatible with
   the current NVECTOR module, or an input value was illegal.

:index:`ARKLS_MEM_FAIL`  (-4):
   A memory allocation request failed.

:index:`ARKLS_PMEM_NULL`  (-5):
   The preconditioner module has not been initialized.

:index:`ARKLS_MASSMEM_NULL`  (-6):
   The ARKLS mass-matrix linear solver interface has not been initialized.

:index:`ARKLS_JACFUNC_UNRECVR`  (-7):
   The Jacobian function failed in an unrecoverable manner.

:index:`ARKLS_JACFUNC_RECVR`  (-8):
   The Jacobian function had a recoverable error.

:index:`ARKLS_MASSFUNC_UNRECVR`  (-9):
   The mass matrix function failed in an unrecoverable manner.

:index:`ARKLS_MASSFUNC_RECVR`  (-10):
   The mass matrix function had a recoverable error.

:index:`ARKLS_SUNMAT_FAIL`  (-11):
   An error occurred with the current SUNMATRIX module.

:index:`ARKLS_SUNLS_FAIL`  (-12):
   An error occurred with the current SUNLINSOL module.
