      program fkinDiagon_kry_p
c     --------------------------------------------------------------------
c     Programmer(s): Allan G. Taylor, Alan C. Hindmarsh and
c                    Radu Serban @ LLNL
c     --------------------------------------------------------------------
C     SUNDIALS Copyright Start
C     Copyright (c) 2002-2022, Lawrence Livermore National Security
C     and Southern Methodist University.
C     All rights reserved.
C
C     See the top-level LICENSE and NOTICE files for details.
C
C     SPDX-License-Identifier: BSD-3-Clause
C     SUNDIALS Copyright End
c     --------------------------------------------------------------------
c     Simple diagonal test with Fortran interface, using user-supplied
c     preconditioner setup and solve routines (supplied in Fortran).
c
c     This example does a basic test of the solver by solving the
c     system:
c               f(u) = 0  for
c               f(u) = u(i)^2 - i^2
c
c     No scaling is done.
c     An approximate diagonal preconditioner is used.
c
c      Execution command: mpirun -np 4 fkinDiagon_kry_p
c     --------------------------------------------------------------------
c
      implicit none

      include "mpif.h"

      integer localsize
      parameter(localsize=32)
      integer baseadd, i, ii
      integer size, rank, mype, npes
c The following declaration specification should match C type long int.
      integer*8 neq, nlocal, iout(16), msbpre
      integer ier, globalstrat, prectype, maxl, maxlrst
      double precision pp, fnormtol, scsteptol
      double precision rout(2), uu(localsize), scale(localsize)
      double precision constr(localsize)

      common /pcom/ pp(localsize), nlocal, mype, npes, baseadd

      nlocal = localsize
      neq = 4 * nlocal
      globalstrat = 0
      fnormtol = 1.0d-5
      scsteptol = 1.0d-4
      prectype = 2
      maxl = 10
      maxlrst = 2
      msbpre  = 5

c     The user MUST call mpi_init, Fortran binding, for the fkinsol package
c     to work. The communicator, MPI_COMM_WORLD, is the only one common 
c     between the Fortran and C bindings. So in the following, the communicator
c     MPI_COMM_WORLD is used in calls to mpi_comm_size and mpi_comm_rank
c     to determine the total number of processors and the rank (0 ... size-1) 
c     number of this process.
      
      call mpi_init(ier)
      if (ier .ne. 0) then
         write(6,1210) ier
 1210    format('MPI_ERROR: MPI_INIT returned IER = ', i4)
         stop
      endif

      call fnvinitp(mpi_comm_world, 3, nlocal, neq, ier)
      if (ier .ne. 0) then
         write(6,1220) ier
 1220    format('SUNDIALS_ERROR: FNVINITP returned IER = ', i4)
         call mpi_finalize(ier)
         stop
      endif
      
      call mpi_comm_size(mpi_comm_world, size, ier)
      if (ier .ne. 0) then
         write(6,1222) ier
 1222    format('MPI_ERROR: MPI_COMM_SIZE returned IER = ', i4)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
      
      if (size .ne. 4) then
         write(6,1223)
 1223    format('MPI_ERROR: must use 4 processes')
         call mpi_finalize(ier)
         stop
      endif
      npes = size

      call mpi_comm_rank(mpi_comm_world, rank, ier)
      if (ier .ne. 0) then
         write(6,1224) ier
 1224    format('MPI_ERROR: MPI_COMM_RANK returned IER = ', i4)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif

      mype = rank
      baseadd = mype * nlocal

      do 20 ii = 1, nlocal
         i = ii + baseadd
         uu(ii) = 2.0d0 * i
         scale(ii) = 1.0d0
         constr(ii) = 0.0d0
 20   continue
      
      call fkincreate(ier)
      if (ier .ne. 0) then
         write(6,1230) ier
 1230    format('SUNDIALS_ERROR: FKINCREATE returned IER = ', i4)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif

      call fkinsetiin('MAX_SETUPS', msbpre, ier)
      if (ier .ne. 0) then
         write(6,1231) ier
 1231    format('SUNDIALS_ERROR: FKINSETIIN returned IER = ', i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif

      call fkinsetrin('FNORM_TOL', fnormtol, ier)
      if (ier .ne. 0) then
         write(6,1232) ier
 1232    format('SUNDIALS_ERROR: FKINSETRIN returned IER = ', i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif

      call fkinsetrin('SSTEP_TOL', scsteptol, ier)
      if (ier .ne. 0) then
         write(6,1232) ier
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif

      call fkinsetvin('CONSTR_VEC', constr, ier)
      if (ier .ne. 0) then
         write(6,1233) ier
 1233    format('SUNDIALS_ERROR: FKINSETVIN returned IER = ', i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c Initialize KINSOL
c
      call fkininit(iout, rout, ier)
      if (ier .ne. 0) then
         write(6,1234) ier
 1234    format('SUNDIALS_ERROR: FKININIT returned IER = ', i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c Initialize SPGMR linear solver module with right preconditioning
c and maximum Krylov dimension maxl
c
      call fsunspgmrinit(3, prectype, maxl, ier)
      if (ier .ne. 0) then
         write(6,1235) ier
 1235    format('SUNDIALS_ERROR: FSUNSPGMRLINSOLINIT returned IER = ',
     1          i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c Attach SPGMR linear solver module to KINSOL
c
      call fkinlsinit(ier)
      if (ier .ne. 0) then
         write(6,1236) ier
 1236    format('SUNDIALS_ERROR: FKINLSINIT returned IER = ', i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c Set the maximum number of SPGMR restarts to maxlrst
c
      call fsunspgmrsetmaxrs(3, maxlrst, ier)
      if (ier .ne. 0) then
         write(6,1237) ier
 1237    format('SUNDIALS_ERROR: FSUNSPGRM_SETMATRS returned IER = ',
     1          i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c Set preconditioner routines
c
      call fkinlssetprec(1, ier)
      if (ier .ne. 0) then
         write(6,1238) ier
 1238    format('SUNDIALS_ERROR: FKINLSSETPREC returned IER = ',
     1          i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif

      if (mype .eq. 0) write(6,1240)
 1240 format('Example program fkinDiagon_kry_p:'//
     1       ' This FKINSOL example',
     2       ' solves a 128 eqn diagonal algebraic system.'/
     3       ' Its purpose is to demonstrate the use of the Fortran',
     4       ' interface'/' in a parallel environment.')

      call fkinsol(uu, globalstrat, scale, scale, ier)
      if (ier .lt. 0) then
         write(6,1242) ier, iout(9)
 1242    format('SUNDIALS_ERROR: FKINSOL returned IER = ', i4, /,
     1          '                Linear Solver returned IER = ', i4)
         call fkinfree
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif

      if (mype .eq. 0) write(6,1245) ier
 1245 format(/' FKINSOL return code is ', i4)

      if (mype .eq. 0) write(6,1246)
 1246 format(/' The resultant values of uu (process 0) are:'/)

      do 30 i = 1, nlocal, 4
         if(mype .eq. 0) write(6,1256) i + baseadd, uu(i), uu(i+1),
     1                                 uu(i+2), uu(i+3)
 1256    format(i4, 4(1x, f10.6))
 30   continue

      if (mype .eq. 0) write(6,1267) iout(3), iout(15), iout(4),
     1                               iout(13), iout(14), iout(16)
 1267 format(/'Final statistics:'//
     1     '  nni = ', i3, ',  nli  = ', i3, /,
     2     '  nfe = ', i3, ',  npe  = ', i3, /,
     3     '  nps = ', i3, ',  ncfl = ', i3)

      call fkinfree

c     An explicit call to mpi_finalize (Fortran binding) is required by
c     the constructs used in fkinsol.
      call mpi_finalize(ier)

      stop
      end
      

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The function defining the system f(u) = 0 must be defined by a Fortran
c     function with the following name and form.
      
      subroutine fkfun(uu, fval, ier)

      implicit none

      integer ier, i
c The following declaration specification should match C type long int.
      integer*8 nlocal
      integer mype, npes, baseadd
      integer localsize
      parameter(localsize=32)
      double precision pp
      double precision fval(*), uu(*)

      common /pcom/ pp(localsize), nlocal, mype, npes, baseadd

      do 10 i = 1, nlocal
 10      fval(i) = uu(i) * uu(i) - (i + baseadd) * (i + baseadd)

      ier = 0

      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpreco is the preconditioner setup routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpset(udata, uscale, fdata, fscale, ier)

      implicit none

      integer ier, i
c The following declaration specification should match C type long int.
      integer*8 nlocal
      integer mype, npes, baseadd
      integer localsize
      parameter(localsize=32)
      double precision pp
      double precision udata(*), uscale(*), fdata(*), fscale(*)

      common /pcom/ pp(localsize), nlocal, mype, npes, baseadd

      do 10 i = 1, nlocal
 10      pp(i) = 0.5d0 / (udata(i) + 5.0d0)

      ier = 0

      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpsol is the preconditioner solve routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpsol(udata, uscale, fdata, fscale, vv, ier)

      implicit none

      integer ier, i
c The following declaration specification should match C type long int.
      integer*8 nlocal
      integer mype, npes, baseadd
      integer localsize
      parameter(localsize=32)
      double precision pp
      double precision udata(*), uscale(*), fdata(*), fscale(*), vv(*)

      common /pcom/ pp(localsize), nlocal, mype, npes, baseadd

      do 10 i = 1, nlocal
 10      vv(i) = vv(i) * pp(i)

      ier = 0

      return
      end
