      program kindiagpf
c ------------------------------------------------------------------
c $Revision: 1.2 $
c $Date: 2004-05-04 16:35:48 $
c ------------------------------------------------------------------
c Programmers : Allan G. Taylor, Alan C. Hindmarsh, and
c               Radu Serban @ LLNL
c ------------------------------------------------------------------
c
c    File: kindiagpf.f                                               
c
c      Simple diagonal test with Fortran interface, using user-supplied      
c      preconditioner setup and solve routines (supplied in Fortran, below). 
c      This example does a basic test of the solver by solving the system    
c                           f(u) = 0   for                                   
c                           f(u) = u(i)^2 - i^2 .                            
c      No scaling is done.                                                   
c      An approximate diagonal preconditioner is used.                       
c      Execute line:  mpirun -np 4 kindiagpf                                 
c
      include "mpif.h"
c
      integer ier, size, globalstrat, rank, inopt
      integer mype, npes, maxl, maxlrst
      integer*8 localsize
      parameter(localsize=32)
      integer*8 neq, nlocal, msbpre, baseadd, i, ii
      integer*8 iopt(40)
      double precision fnormtol, scsteptol
      double precision ropt(40)
      double precision uu(localsize), scale(localsize)
      double precision constr(localsize)
      double precision pp
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c
      nlocal = localsize
      neq = 4*nlocal
      globalstrat = 0
      fnormtol = 0.00001
      scsteptol = 0.0001
      inopt = 0
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
 1210    format('MPI_ERROR: MPI_INIT returned IER =',i2)
         stop
      endif

      call fnvinitp(nlocal, neq, ier)
      if (ier .ne. 0) then
         write(6,1220) ier
 1220    format('SUNDIALS_ERROR: FNVINITP returned IER =',i2)
         call mpi_finalize(ier)
         stop
      endif
      
      call mpi_comm_size(MPI_COMM_WORLD,size,ier)
      if (ier .ne. 0) then
         write(6,1222) ier
 1222    format('MPI_ERROR: MPI_COMM_SIZE returned IER =',i2)
         call mpi_abort(MPI_COMM_WORLD, 1, ier)
         stop
      endif
      
      if (size .ne. 4) then
         write(6,1230)
 1230    format('MPI_ERROR: must use 4 processes')
         call mpi_finalize(ier)
         stop
      endif
      npes = size

      call mpi_comm_rank(MPI_COMM_WORLD, rank, ier)
      if (ier .ne. 0) then
         write(6,1224) ier
 1224    format('MPI_ERROR: MPI_COMM_RANK returned IER =',i2)
         call mpi_abort(MPI_COMM_WORLD, 1, ier)
         stop
      endif

      mype = rank
      baseadd = mype * nlocal 
      
      do 20 ii = 1,nlocal
         i = ii + baseadd
         uu(ii) = 2.0*i
         scale(ii) = 1.0
         constr(ii) = 0.0
 20   continue
      
      call fkinmalloc(msbpre, fnormtol, scsteptol, 
     &                constr, inopt, iopt, ropt, ier)
      
      if (ier .ne. 0) then
         write(6,1231)ier
 1231    format('SUNDIALS_ERROR: FKINMALLOC returned IER =',i2)
         call mpi_abort(MPI_COMM_WORLD, 1, ier)
         stop
      endif
      
      call fkinspgmr(maxl, maxlrst, ier)
      call fkinspgmrsetpsol(1, ier)
      call fkinspgmrsetpset(1, ier)
      
      if(mype .eq. 0)write(6,1240)
 1240 format('Example program kindiagpf'/' This fkinsol example code',
     1       ' solves a 128 eqn diagonal algebraic system.'/
     2       ' Its purpose is to demonstrate the use of the Fortran',
     3       ' interface'/' in a parallel environment.'/
     4       ' globalstrategy = INEXACT_NEWTON'/)

      call fkinsol(uu, globalstrat, scale, scale, ier)
      if (ier .lt. 0) then
         write(6,1242) ier
 1242    format('SUNDIALS_ERROR: FKINSOL returned IER =',i2)
         call mpi_abort(MPI_COMM_WORLD, 1, ier)
         stop
      endif

      if (mype .eq. 0) write(6,1245)ier
 1245 format(/' fkinsol return code is ',i5)

      if (mype .eq. 0) write(6,1246)
 1246 format(/' The resultant values of uu (processor 0) are:'//)
      
      do 30 i = 1,nlocal,4
         if(mype .eq. 0) write(6,1256)i+baseadd, uu(i), uu(i+1),
     1        uu(i+2), uu(i+3)
 1256    format(i4,4(1x,f10.6))
 30   continue

      if (mype .eq. 0) write(6,1267)iopt(4),iopt(11),iopt(5),iopt(12),
     1                 iopt(13),iopt(14)
 1267 format(//' nni=',i4,',  nli=',i4,',  nfe=',i4,',  npe=',i4,
     1       ',  nps=',i4,',  ncfl=',i4)

      call fkinfree
      call fnvfreep
      
c     An explicit call to mpi_finalize (Fortran binding) is required by 
c     the constructs used in fkinsol. 
      call mpi_finalize(ier)
      
      stop
      end
      

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The function defining the system f(u) = 0 must be defined by a Fortran
c     function with the following name and form. 
      
      subroutine fkfun(uu, fval)
      implicit none
      integer mype, npes
      integer*8 baseadd, nlocal, i
      integer*8 localsize
      parameter(localsize=32)
      double precision fval(*), uu(*)
      double precision pp
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c
      do 10 i = 1,nlocal
 10      fval(i) = uu(i)*uu(i) - (i+baseadd)*(i+baseadd)
      
      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpreco is the preconditioner setup routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpset(udata, uscale, fdata, fscale, 
     1                  vtemp1, vtemp2, ier)
      implicit none
      integer ier, mype, npes
      integer*8 baseadd, nlocal, i
      integer*8 localsize
      parameter(localsize=32)
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)
      double precision pp
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c
      do 10 i = 1,nlocal
 10      pp(i) = 0.5/(udata(i)+5.0)

      ier = 0
      
      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpsol is the preconditioner solve routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpsol(udata, uscale, fdata, fscale, 
     1                  vv, ftem, ier)
      implicit none
      integer ier, mype, npes
      integer*8 baseadd, nlocal, i
      integer*8 localsize
      parameter(localsize=32)
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vv(*), ftem(*)
      double precision pp
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c
      do 10 i = 1,nlocal
 10      vv(i) = vv(i) * pp(i)

      ier = 0

      return
      end
