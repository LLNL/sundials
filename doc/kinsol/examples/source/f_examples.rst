.. _KINSOL.Examples.Fortran:

Fortran example problems
========================

The Fortran example problem programs supplied with the KINSOL package are all
written using the Fortran 2003 standard, use double precision arithmetic, and
32- or 64-bit integers for indexing. See the :ref:`SUNDIALS.Fortran` section in
the KINSOL User Documentation for details.

A serial example: fkinDiagon_kry
--------------------------------

The ``fkinDiagon_kry`` program solves the problem

.. math::

   F(u) = 0, \quad \text{where } f_i(u) = u_i^2 - i^2, \quad 1 \le i \le N.

using the serial vector module.

The main program begins by calling ``fnvinits`` to initialize computations with
the nvecs module. Next, the array ``uu`` is set to contain the initial guess
:math:`u_i = 2 i`, the array ``scale`` is set with all components equal to
``1.0`` (meaning that no scaling is done), and the array ``constr`` is set with
all components equal to ``0.0`` to indicate that no inequality constraints
should be imposed on the solution vector.

The KINSOL solver is initialized and memory for it is allocated by calling
``fkinmalloc``, which also specifies the ``iout`` and ``rout`` arrays which are
used to store integer and real outputs, respectively (see Table
``t:fkinsol_out``). Also, various integer, real, and vector parameters are
specified by calling the ``fkinsetiin``, ``fkinsetrin``, and ``fkinsetvin``
subroutines, respectively. In particular, the maximum number of iterations
between calls to the preconditioner setup routine (``msbpre = 5``), the
tolerance for stopping based on the function norm (``fnormtol = 1.0e-5``), and
the tolerance for stopping based on the step length (``scsteptol = 1.0e-4``) are
specified.

Next, the ``sunlinsolspgmr`` linear solver module is attached to KINSOL by
calling ``fsunspgmrinit``, which also specifies the maximum Krylov subspace
dimension (``maxl = 10``). This is then attached to KINSOL by calling
``fkinlsinit``. The maximum number of restarts allowed for SPGMR is then updated
to ``maxlrst = 2`` by calling ``fsunspgmrsetmaxrs``. The ``sunlinsolspgmr``
module is then directed to use the supplied preconditioner by calling the
``fkinlssetprec`` routine with a first argument equal to ``1``. The solution of
the nonlinear system is obtained after a successful return from ``fkinsol``,
which is then printed to unit 6 (stdout). Finally, memory allocated for the
KINSOL solver is released by calling ``fkinfree``.

The user-supplied routine ``fkfun`` contains a straightforward transcription of
the nonlinear system function :math:`f`, while the routine ``fkpset`` sets the
array ``pp`` (in the common block ``pcom``) to contain an approximation to the
reciprocals of the Jacobian diagonal elements. The components of ``pp`` are then
used in ``fkpsol`` to solve the preconditioner linear system :math:`P x = v`
through simple multiplications.

The following is sample output from ``fkinDiagon_kry``, using :math:`N = 128`.

.. literalinclude:: ../../../../examples/kinsol/F2003_serial/kinDiagon_kry_f2003.out
   :language: none

A parallel example: fkinDiagon_kry_p
------------------------------------

The program ``fkinDiagon_kry_p`` is a straightforward modification of
``fkinDiagon_kry`` to use the MPI vector module.

After initialization of MPI, the nvecp module is initialized by calling
``fnvinitp`` with the default MPI communicator ``mpi_comm_world`` and the local
and global vector sizes as its first three arguments. The rank of the local
process, ``mype``, is used in both the initial guess and the system function,
inasmuch as the global and local indices to the vector ``u`` are related by the
equation ``iglobal = ilocal + mype*nlocal``. In other respects, the problem
setup (KINSOL initialization, ``sunlinsolspgmr`` specification) and solution
steps are the same as in ``fkinDiagon_kry``. Upon successful return from
``fkinsol``, the solution segment local to the process with id equal to 0 is
printed to unit 6. Finally, the KINSOL memory is released and the MPI
environment is terminated.

For this simple example, no inter-process communication is required to evaluate
the nonlinear system function :math:`f` or the preconditioner. As a consequence,
the user-supplied routines ``fkfun``, ``fkpset``, and ``fkpsol`` are basically
identical to those in ``fkinDiagon_kry``.

Sample output from ``fkinDiagon_kry_p``, for :math:`N = 128`, follows.

.. literalinclude:: ../../../../examples/kinsol/F2003_parallel/kin_diagon_kry_f2003_8.out
   :language: none
