      program diagpf
c   ***************************************************************************
c   *                                                                         *
c   * File: kindiagpf.f                                                       *
c   * Programmers: Allan G. Taylor and Alan C. Hindmarsh @ LLNL               *
c   * Version of 25 Mar 1999                                                  *
c   *   simple diagonal test with Fortran interface, using user-supplied      *
c   *   preconditioner setup and solve routines (supplied in Fortran, below)  *
c   *   This example does a basic test of the code suite by solving the system*
c   *                        f(u) = 0.  for                                   *
c   *                        f(u) = u(i)^2 - i^2 .                            *
c   *   No scaling is done.                                                   *
c   *   execute line:  mpirun -np 4 kindiagpf                                 *
c   *--------------------------------------------------------------------------
c     the following include is required
      include "mpif.h"

      parameter(localsize=32)
      integer neq, ier, size, globalstrat, rank
      integer mype, npes, baseadd, nlocal, comm
      integer iopt(40)
      
      double precision fnormtol, scsteptol
      double precision ropt(40)
      double precision uu(localsize), scale(localsize)
      double precision constr(localsize)
      double precision ri
      double precision pp
      
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal

      comm = 0
      neq = 4*localsize
      nlocal = 32
      globalstrat = 0
      fnormtol = 0.00001
      scsteptol = 0.0001
      inopt = 0
      maxl = 10
      maxlrst = 2
      msbpre  = 5
      mu = 0
      ml = 0

c    the user MUST call mpi_init , Fortran binding, for the fkinsol package
c    to work. The communicator, MPI_COMM_WORLD, is the only one common 
c    between the Fortran and C bindings. So, in the following, the communicator
c    MPI_COMM_WORLD is used in calls to mpi_comm_size and mpi_comm_rank
c    to determine the total number of processors and the rank (0 ... size-1) 
c    number of this process.

      call mpi_init(ierr)
c * * * * * * * * * * * * * * * * * * * * * * *
     
      call fkinitmpi(nlocal, neq, ier)

      call mpi_comm_size(MPI_COMM_WORLD,size,ier)

       if(size .ne. 4) then
         write(6,1234)
 1234    format(1x,' # of processors not 4. set to 4 and try again')
         call fkfreempi
         stop
       endif
       npes = size
       call mpi_comm_rank(MPI_COMM_WORLD, rank, ier)

       mype = rank
       baseadd = mype * nlocal 

       call fpkinmalloc(neq, ier)

       if(ier .ne. 0)then
          write(6,1230)ier
 1230     format('fpkinmalloc failed',i2)
          stop
       endif

       call fkinspgmr20(maxl, maxlrst, msbpre)

       if(mype .eq. 0)write(6,1240)
 1240  format(//,'kindiagpf demo case. This kinsol example code does a')
       if(mype .eq. 0) write(6,1241)
 1241  format('128 eqn diagonal algebraic system. Its purpose is to ')
       if(mype .eq. 0)write(6,1242)
 1242  format('demonstrate the use of the Fortran interface in a')
       if(mype .eq. 0)write(6,1243)
 1243  format('parallel environment.',//)
       if(mype .eq. 0)write(6,1244)
 1244  format('globalstrategy = INEXACT_NEWTON' ,//)

      do 20 ii = 1 , 32
         i = ii + baseadd
         ri = i
         uu(ii) = 2*ri
         scale(ii) = 1.
         constr(ii) = 0.
 20   continue

      call fkinsol(neq, uu, 0, scale, scale, fnormtol, 
     1              scsteptol, constr, inopt, iopt, ropt, ier)


      if(mype .eq. 0) write(6,1245)ier
 1245 format(1x,' return code is ',i5)


      if(mype .eq. 0)write(6,1246)
 1246 format(1x,'The resultant values of uu (processor 0) are:',///)

      do 30 i = 1, 32, 4
        if(mype .eq. 0) write(6,1256)i+baseadd, uu(i), uu(i+1),
     1                                          uu(i+2), uu(i+3)
 1256    format(i4,4(1x,f10.6))
 30   continue


      if(mype .eq. 0) write(6,1267) iopt(4), iopt(11), iopt(5), 
     1    iopt(12), iopt(13),  iopt(14)
 1267 format(1x,'nni=',i4,' nli=',i4,' nfe=',i4,' npe=',i4,
     1       ' nps=',i4,' ncfl=',i4)

      call fkinfree

      call fkfreempi

c     an explicit call to mpi_finalize (Fortran binding) is required by 
c     the constructs used in fkinsol. 

      call mpi_finalize(ier)

      stop

      end

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The function defining the system f(u) = 0. must be defined by a Fortran
c     function of the following form. 

      subroutine KFUN(nloc, uu, fval)
c
      double precision fval(*), uu(*)

c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      parameter(localsize=32)
      double precision pp
      integer mype, npes, baseadd, nlocal
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *


      do 10 i = 1 , nlocal
         fval(i) = uu(i)*uu(i) - (i+baseadd)*(i+baseadd)
 10   continue

      return
      end


c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   the routine kpreco is the preconditioner setup routine. It is required that
c   that specific name be used in order that the c code can find and link to it
c   the argument list must also be as illustrated below:

      subroutine kpreco(neq, udata, uscale, fdata, fscale, 
     1                  vtemp1, vtemp2, uround, nfe, ier)

      integer nfe, ier

      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*) 

      double precision uround
      
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      parameter(localsize=32)
      double precision pp
      integer mype, npes, baseadd, nlocal
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do 10 i = 1 , nlocal
         pp(i) = 0.5/(udata(i)+5.)
 10   continue

      return
      end


c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   the routine kpsol is the preconditioner solve routine. It is required that
c   that specific name be used in order that the c code can find and link to it
c   the argument list must also be as illustrated below:

      subroutine kpsol(neq, udata, uscale, fdata, fscale, 
     1                 vtem, ftem, uround, nfe, ier)

      integer nfe, ier

      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtem(*), ftem(*) 

      double precision uround
      
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      parameter(localsize=32)
      double precision pp
      integer mype, npes, baseadd, nlocal
      common /pcom/ pp(localsize), mype, npes, baseadd, nlocal
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do 10 i = 1 , nlocal
         vtem(i) = vtem(i) * pp(i)
 10   continue

      return
      end
