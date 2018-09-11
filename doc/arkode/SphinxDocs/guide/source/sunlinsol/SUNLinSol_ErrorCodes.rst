..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2017, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _SUNLinSol.ErrorCodes:

Error Codes returned from SUNLinearSolver implementations
============================================================

The functions within the SUNDIALS-provided ``SUNLinearSolver``
implementations return a common set of error codes, listed here:

* ``SUNLS_SUCCESS`` (0) -- successful call or converged solve
* ``SUNLS_MEM_NULL`` (-1) -- the memory argument to the function is ``NULL``
* ``SUNLS_ILL_INPUT`` (-2) -- an illegal input has been provided to the function 
* ``SUNLS_MEM_FAIL`` (-3) -- failed memory access or allocation
* ``SUNLS_ATIMES_FAIL_UNREC`` (-4) -- an unrecoverable failure occurred in the ``ATimes`` routine
* ``SUNLS_PSET_FAIL_UNREC`` (-5) -- an unrecoverable failure occurred in the ``Pset`` routine
* ``SUNLS_PSOLVE_FAIL_UNREC`` (-6) -- an unrecoverable failure occurred in the ``Psolve`` routine
* ``SUNLS_PACKAGE_FAIL_UNREC`` (-7) -- an unrecoverable failure occurred in an external linear solver package
* ``SUNLS_GS_FAIL`` (-8) -- a failure occurred during Gram-Schmidt orthogonalization (SPGMR/SPFGMR)
* ``SUNLS_QRSOL_FAIL`` (-9) -- a singular $R$ matrix was encountered in a QR factorization (SPGMR/SPFGMR)
* ``SUNLS_RES_REDUCED`` (1) -- an iterative solver reduced the residual, but did not converge to the desired tolerance
* ``SUNLS_CONV_FAIL`` (2) -- an iterative solver did not converge (and the residual was not reduced)
* ``SUNLS_ATIMES_FAIL_REC`` (3) -- a recoverable failure occurred in the ``ATimes`` routine
* ``SUNLS_PSET_FAIL_REC`` (4) -- a recoverable failure occurred in the ``Pset`` routine
* ``SUNLS_PSOLVE_FAIL_REC`` (5) -- a recoverable failure occurred in the ``Psolve`` routine
* ``SUNLS_PACKAGE_FAIL_REC`` (6) -- a recoverable failure occurred in an external linear solver package
* ``SUNLS_QRFACT_FAIL`` (7) -- a singular matrix was encountered during a QR factorization (SPGMR/SPFGMR)
* ``SUNLS_LUFACT_FAIL`` (8) -- a singular matrix was encountered during a LU factorization
