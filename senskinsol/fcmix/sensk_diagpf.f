      program senskdiagpf
c   ***************************************************************************
c   *                                                                         *
c   * File: sensk_diagpf.f                                                    *
c   * Programmers: Allan G. Taylor, Alan C. Hindmarsh, and                    *
c   *              Keith E. Grant @ LLNL                                      *
c   *                                                                         *
c   * Version of 29 Nov 2000                                                  *
c   *   Modified as a FORTRAN interface example of parameter sensitivity      *
c   *   analysis. A parameter of nominal magnitude one is multiplied into the *
c   *   'i^2' term.                                                           *
c   *                                                                         *
c   * Version of 25 Mar 1999                                                  *
c   *   simple diagonal test with Fortran interface, using user-supplied      *
c   *   preconditioner setup and solve routines (supplied in Fortran, below)  *
c   *   This example does a basic test of the code suite by solving the       *
c   *   system                                                                *
c   *                        f(u) = 0.  for                                   *
c   *                        f(u) = u(i)^2 - p*(i^2) .                        *
c   *   No scaling is done.                                                   *
c   *                                                                         *
c   *   execute line:  mpirun -np 4 sensk_diagpf                              *
c   *                                                                         *
c   ***************************************************************************

      implicit none
c
c.... The following include is required
      include "mpif.h"

      integer localsize, nopt
      parameter(localsize=32)
      parameter(nopt=40)

      integer mype, npes, baseadd, nlocal, comm, rank

      integer neq, ier, size, globalstrat
      integer inopt, maxl, maxlrst, msbpre, mu, ml
      integer i, ii, ip

c **********************************************************************
c.... IOPT is an array of length 40 for integer optional inputs and
c.... outputs from FSENSKINOL. Declare it as INTEGER*4 or INTEGER*8
c.... according to the size of C type "long int" on the current machine
c.... This length can be determined with the following C program:
c
c        #include <stdio.h>
c        main () {printf("%d\n", sizeof(long int));}
c
c     integer*4 iopt(nopt)
      integer*8 iopt(nopt)
c **********************************************************************

      double precision fnormtol, scsteptol
      double precision ropt(nopt)
      double precision uu(localsize), scale(localsize)
      double precision constr(localsize)
      double precision ri
      double precision pp

      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal

c
c.... Define parameter and nominal magnitude variables for
c.... sensitiivty analysis
      integer np, idiffrhs
      double precision diffscl
      parameter (
     &   np = 1,
     &   idiffrhs = 0,
     &   diffscl = 1.0)
      double precision params(np), pbar(np), ww(localsize)

      comm = 0
      neq = 4*localsize
      nlocal = localsize
      globalstrat = 0
      fnormtol = 0.00001
      scsteptol = 0.0001
      inopt = 0
      maxl = 10
      maxlrst = 2
      msbpre  = 5
      mu = 0
      ml = 0

c
c.... Define the parameter and nominal magnitude values for
c.... sensitivity analysis
      do ip = 1,np
         params(ip) = 1.0
         pbar(ip)   = 1.0
      enddo

c     the user MUST call mpi_init , Fortran binding, for the fkinsol
c     package to work. The communicator, MPI_COMM_WORLD, is the only
c     one common between the Fortran and C bindings. So, in the
c     following, the communicator MPI_COMM_WORLD is used in calls to
c     mpi_comm_size and mpi_comm_rank to determine the total number of
c     processors and the rank (0 ... size-1) number of this process.

      call mpi_init(ier)

c
c.... Initialize KINSOL usage of MPI
      call fkinitmpi(nlocal, neq, ier)

c
c.... Verify that the number of processors allocated matches what
c.... this program expects
      call mpi_comm_size(MPI_COMM_WORLD,size,ier)
      if (size .ne. 4) then
         write(6,'(a/)')
     &      ' *** # of processors not 4. set to 4 and try again'
         call fkfreempi
         stop
      endif

      npes = size
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ier)
      mype = rank
      baseadd = mype * nlocal

c
c.... Allocate and initialize SensKINSOL memory
      call fpsenskinmalloc(neq, np, idiffrhs, diffscl, params,
     &   pbar, ier)
      if (ier .ne. 0) then
         write(6,'(a,i2/)') ' *** FPSENSKINMALLOC failed with ', ier
         stop
      endif

c
c.... Initialize the linear solver
      call fsenskinspgmr20(maxl, maxlrst, msbpre)

      if (mype .eq. 0) then
         write(6,'(//a/a/a/a//)')
     &      'SENS_DIAGPF demo case. This kinsol example code does ',
     &      'a 128 eqn diagonal algebraic system. Its purpose is to ',
     &      'demonstrate the use of the Fortran interface in a',
     &      'parallel environment.'
         write(6, '(a//)') ' globalstrategy = INEXACT_NEWTON'
      endif

      do ii = 1,localsize
         ri = ii + baseadd
         uu(ii) = 2*ri
         scale(ii) = 1.
         constr(ii) = 0.
      enddo

c
c.... Do the nonlinear system solution
      call fsenskinsol(neq, uu, 0, scale, scale, fnormtol,
     &   scsteptol, constr, inopt, iopt, ropt, ier)


c
c.... Output the nonlinear system results for processor zero
      if(mype .eq. 0) then

         write(6,'(a,i5)') ' FSENSKINSOL returned flag = ', ier
         write(6,'(a/)')
     &      ' The resultant values of uu (processor 0) are:'

         do i = 1,nlocal,4
            write(6,'(i4,4(1x,f10.6))') i+baseadd,
     &         uu(i), uu(i+1), uu(i+2), uu(i+3)
         enddo

         write(6,'(/6(1x,a,i4))')
     &      'nni=', iopt(4),  'nli=', iopt(11), 'nfe=',  iopt(5),
     &      'npe=', iopt(12), 'nps=', iopt(13), 'ncfl=', iopt(14)
      endif

c
c.... Initialize for parameter sensitivity analysis
      call fsenskinlininit (ier)

      do ip=1,np
c
c....    Calculate the sensitivity for the ip'th parameter
         call fsenskinlinsolve (ip, ww, ier)

c
c....    Write out the ip'th sensitivity vector
         if(mype .eq. 0) then
            write (6,'(//a/)')
     &      'The Sensitivity values s (processor 0) are:'
            do i = 1,nlocal,4
               write(6,'(i4,4(1x,f10.6))') i+baseadd,
     &            ww(i)/pbar(ip), ww(i+1)/pbar(ip), 
     &            ww(i+2)/pbar(ip), ww(i+3)/pbar(ip)
            enddo

         write(6,'(/6(1x,a,i4))')
     &      'nni=', iopt(4),  'nli=', iopt(11), 'nfe=',  iopt(5),
     &      'npe=', iopt(12), 'nps=', iopt(13), 'ncfl=', iopt(14)
         endif

      enddo


      call fsenskinfree
      call fkfreempi
c
c.... An explicit call to mpi_finalize (Fortran binding) is required by
c.... the constructs used in fkinsol.
      call mpi_finalize(ier)

      stop
      end

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The function defining the system f(u) = 0. must be defined by a Fortran
c     function of the following form. To use this for a sensitivity example,
c     a single scaling parameter of nominal value 1.0 is multiplied into the
c     'i^2' term.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      subroutine KFUN(nloc, uu, fval, np, params)

      implicit none

      integer nloc, np
      double precision fval(*), uu(*), params(np)

      integer i

c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      integer localsize
      parameter(localsize=32)
      double precision pp
      integer mype, npes, baseadd, nlocal
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *


      do i = 1,nlocal
         fval(i) = uu(i)*uu(i) - params(1)*(i+baseadd)*(i+baseadd)
      enddo

      return
      end


c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   the routine kpreco is the preconditioner setup routine. It is required that
c   that specific name be used in order that the c code can find and link to it
c   the argument list must also be as illustrated below:

      subroutine kpreco(neq, udata, uscale, fdata, fscale,
     1                  vtemp1, vtemp2, uround, nfe, ier)

      integer neq, nfe, ier

      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)
      double precision uround

      integer i

c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      integer localsize
      parameter(localsize=32)
      double precision pp
      integer mype, npes, baseadd, nlocal
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do i = 1,nlocal
         pp(i) = 0.5/(udata(i)+5.)
      enddo

      return
      end


c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   the routine kpsol is the preconditioner solve routine. It is required that
c   that specific name be used in order that the c code can find and link to it
c   the argument list must also be as illustrated below:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      subroutine kpsol(neq, udata, uscale, fdata, fscale,
     &   vtem, ftem, uround, nfe, ier)

      implicit none

      integer neq, nfe, ier

      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtem(*), ftem(*)
      double precision uround

      integer i

c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      integer localsize
      parameter(localsize=32)
      double precision pp
      integer mype, npes, baseadd, nlocal
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do i = 1,nlocal
         vtem(i) = vtem(i) * pp(i)
      enddo

      return
      end
