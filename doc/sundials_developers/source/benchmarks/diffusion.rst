..
   Author(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _Benchmarks.Diffusion:


Diffusion Benchmark
-------------------

This benchmark problem implements a 2D diffusion equation using either MPI,
MPI + CUDA, or MPI + HIP parallelism. Note a GPU-aware MPI implementation is
required.


Problem description
^^^^^^^^^^^^^^^^^^^

This code simulates the anisotropic 2D heat equation,

.. math::

    u_t = \nabla \cdot (D \nabla u) + b(t,\mathbf{x}),

where :math:`D` is a diagonal matrix with entries :math:`k_x` and :math:`k_y`.
The system is evolved for :math:`t \in [0, t_f]` on the rectangular domain
:math:`(x,y) \equiv \mathbf{x} \in [\mathbf{0}, \mathbf{x_{\text{max}}}]^2`,
with the initial condition

.. math::

   u(0,\mathbf{x}) = \sin^2(\pi x) \sin^2(\pi y),

and stationary boundary conditions

.. math::

   u_t(t,0,y) = u_t(t,x_{\text{max}},y) =
   u_t(t,x,0) = u_t(t,x,y_{\text{max}}) = 0.

The source term is given by

.. math::

   b(t,\mathbf{x}) = & -2 \pi \sin^2(\pi x) \sin^2(\pi y) \sin(\pi t) \cos(\pi t) \\
   & - k_x 2 \pi^2 (\cos^2(\pi x) - \sin^2(\pi x)) \sin^2(\pi y) \cos^2(\pi t) \\
   & - k_y 2 \pi^2 (\cos^2(\pi y) - \sin^2(\pi y)) \sin^2(\pi x) \cos^2(\pi t).

Under this setup, the problem has the analytical solution

.. math::

   u(t,\mathbf{x}) = \sin^2(\pi x) \sin^2(\pi y) \cos^2(\pi t).

Spatial derivatives are computed using second-order centered differences on a
uniform spatial grid. The problem can be evolved in time with ARKODE, CVODE, or
IDA. With ARKODE, an adaptive step diagonally implicit Runge-Kutta (DIRK) method
is applied. When using CVODE or IDA, adaptive order and step BDF methods are
used.

In all cases, the nonlinear system(s) in each time step are solved using an
inexact Newton method paired with a matrix-free PCG or GMRES linear solver and a
Jacobi preconditioner.


Options
^^^^^^^

Several command line options are available to change the problem parameters
as well as the integrator and solver options. A summary of the options are
listed in :numref:`Benchmarks.Table.2D_diffusion_options`.

.. tabularcolumns:: |\Y{0.32}|\Y{0.53}|\Y{0.15}|

.. _Benchmarks.Table.2D_diffusion_options:

.. Table:: 2D Diffusion Benchmark Command Line Options

   +-------------------------------+--------------------------------------+---------------+
   | Option                        | Description                          | Default       |
   +===============================+======================================+===============+
   | ``--help``                    | Print the command line options       | --            |
   |                               | and description                      |               |
   +-------------------------------+--------------------------------------+---------------+
   | Problem Configuration Options                                                        |
   +-------------------------------+--------------------------------------+---------------+
   | ``--npx <int>``               | Number of MPI tasks in the           | 0             |
   |                               | x-direction (0 forces MPI to decide) |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--npy <int>``               | Number of MPI tasks in the           | 0             |
   |                               | y-direction (0 forces MPI to decide) |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--nx <int>``                | Number of mesh points in the         | 32            |
   |                               | x-direction                          |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--ny <int>``                | Number of mesh points in the         | 32            |
   |                               | y-direction                          |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--xu <realtype>``           | The domain upper bound in the        | 1.0           |
   |                               | x-direction (:math:`x_\text{max}`)   |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--yu <realtype>``           | The domain upper bound in the        | 1.0           |
   |                               | y-direction :math:`y_\text{max}`     |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--kx <realtype>``           | Diffusion coefficient in the         | 1.0           |
   |                               | x-direction :math:`k_x`              |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--ky <realtype>``           | Diffusion coefficient in the         | 1.0           |
   |                               | y-direction :math:`k_y`              |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--tf <realtype>``           | The final time :math:`t_f`           | 1.0           |
   +-------------------------------+--------------------------------------+---------------+
   | ``--noforcing``               | Disable the forcing term             | Enabled       |
   +-------------------------------+--------------------------------------+---------------+
   | Output Options                                                                       |
   +-------------------------------+--------------------------------------+---------------+
   | ``--output <int>``            | Output level: ``0`` no output,       | 1             |
   |                               | ``1`` output progress and stats,     |               |
   |                               | ``2`` write solution to disk         |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--nout <int>``              | Number of output times               | 20            |
   +-------------------------------+--------------------------------------+---------------+
   | Common Integrator and Solver Options                                                 |
   +-------------------------------+--------------------------------------+---------------+
   | ``--rtol <realtype>``         | Relative tolerance                   | 1e-5          |
   +-------------------------------+--------------------------------------+---------------+
   | ``--atol <realtype>``         | Absolute tolerance                   | 1e-10         |
   +-------------------------------+--------------------------------------+---------------+
   | ``--maxsteps <int>``          | Max number of steps between outputs  | 0             |
   |                               | (0 uses the integrator default)      |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--onstep <int>``            | Number of steps to run using         | 0             |
   |                               | ``ONE_STEP`` mode for debugging      |               |
   |                               | (0 uses ``NORMAL`` mode)             |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--gmres``                   | Use GMRES rather than PCG            | PCG           |
   +-------------------------------+--------------------------------------+---------------+
   | ``--lsinfo``                  | Output linear solver diagnostics     | Off           |
   +-------------------------------+--------------------------------------+---------------+
   | ``--liniters <int>``          | Number of linear iterations          | 20            |
   +-------------------------------+--------------------------------------+---------------+
   | ``--epslin <realtype>``       | Linear solve tolerance factor        | 0             |
   |                               | (0 uses the integrator default)      |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--msbp <int>``              | The linear solver setup frequency    | 0             |
   |                               | (CVODE and ARKODE only, 0 uses the   |               |
   |                               | integrator default)                  |               |
   +-------------------------------+--------------------------------------+---------------+
   | Additional ARKODE Options                                                            |
   +-------------------------------+--------------------------------------+---------------+
   | ``--order <int>``             | Methods order                        | 3             |
   +-------------------------------+--------------------------------------+---------------+
   | ``--controller <int>``        | Error controller option              | 0             |
   +-------------------------------+--------------------------------------+---------------+
   | ``--nonlinear``               | Treat the problem as nonlinearly     | Linear        |
   |                               | implicit                             |               |
   +-------------------------------+--------------------------------------+---------------+
   | ``--diagnostics``             | Output integrator diagnostics        | Off           |
   +-------------------------------+--------------------------------------+---------------+


Building and Running
^^^^^^^^^^^^^^^^^^^^

To build the benchmark executables SUNDIALS should be configured with ARKODE,
CVODE, or IDA enabled and with MPI support on. Additionally, either CUDA or HIP
support must be on to build executables utilizing NVIDIA or AMD GPUs. See the
installation guide for more details on configuring, building, and installing
SUNDIALS.

Based on the configuration, executables for each integrator and backend option
are built and installed in the ``<install prefix>/bin/benchmarks/diffusion_2D``
directory. The executables follow the naming convention
``<package>_diffusion_2D_<parallelism>`` where ``<package>`` is ``arkode``,
``cvode``, or ``ida`` and ``<parallelism>`` is ``mpi`` for MPI only parallelism,
``mpicuda`` for MPI + CUDA, and ``mpihip`` for MPI + HIP.

On Summit, with the default environment

* Compiler: xl/16.1.1-5
* MPI: spectrum-mpi/10.3.1.2-20200121
* CUDA: cuda/10.1.243

an example ``jsrun`` command using CUDA-aware MPI is

.. code-block:: none

   jsrun --smpiargs="-gpu" -n 2 -a 1 -c 1 -g 1 ./cvode_diffusion_2D_mpicuda


On Lassen, with the environment

* Compiler: gcc/8.3.1
* MPI: mvapich2/2021.05.28-cuda-11.1.1
* CUDA: cuda/11.1.1

an example ``jsrun`` command using CUDA-aware MPI

.. code-block:: none

   jsrun -n 2 -a 1 -c 1 -g 1 ./cvode_diffusion_2D_mpicuda
