..
   Copyright (c) 2002-2026, Lawrence Livermore National Security and
   the SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause

.. _CVODE.Examples.Accelerators:

External and accelerator vector examples
========================================

hypre: cvAdvDiff_non_ph
-----------------------

This example solves the same nonstiff problem as ``cvAdvDiff_non_p`` using a
hypre ParVector wrapped as an ``N_Vector``. hypre owns the vector storage and
partition; the wrapper is destroyed before the underlying hypre objects.

CUDA: cvAdvDiff_kry_cuda
------------------------

This example solves the 1D advection-diffusion problem with a CUDA vector and
unpreconditioned GMRES. The right-hand side is implemented by a CUDA kernel;
host access is used only for initialization and output.

.. literalinclude:: ../../../../examples/cvode/cuda/cvAdvDiff_kry_cuda.out
   :language: none

RAJA: cvAdvDiff_kry_raja
------------------------

This C++ example solves the same problem with the RAJA vector and
unpreconditioned GMRES. RAJA execution policies provide the device-independent
right-hand-side kernel.

.. literalinclude:: ../../../../examples/cvode/raja/cvAdvDiff_kry_raja.out
   :language: none
