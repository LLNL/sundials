      program diagsf
c   ***************************************************************************
c   *                                                                         *
c   * File: kindiagsf.f                                                       *
c   * Programmers: Allan G. Taylor and Alan C. Hindmarsh @ LLNL               *
c   * Version of 25 Mar 1999                                                  *
c   *   simple diagonal test with Fortran interface, using user-supplied      *
c   *   preconditioner setup and solve routines (supplied in Fortran, below)  *
c   *   This example does a basic test of the code suite by solving the system*
c   *                        f(u) = 0.  for                                   *
c   *                        f(u) = u(i)^2 - i^2 .                            *
c   *   No scaling is done.                                                   *
c   *   execute line: kindiagsf                                               *
c   *--------------------------------------------------------------------------

      parameter(PROBSIZE=128)
      integer neq, ier, globalstrat
      integer iopt(40)
      
      double precision fnormtol, scsteptol
      double precision ropt(40)
      double precision uu(PROBSIZE), scale(PROBSIZE)
      double precision constr(PROBSIZE)
      double precision pp
      
      common /pcom/ pp(PROBSIZE)

      neq = PROBSIZE
      globalstrat = 0
      fnormtol = 0.00001
      scsteptol = 0.0001
      inopt = 0
      maxl = 10
      maxlrst = 2
      msbpre  = 5
      mu = 0
      ml = 0

c * * * * * * * * * * * * * * * * * * * * * *
      call fskinmalloc(neq, ier)

      if(ier .ne. 0)then
         
 1230    format('fskinmalloc failed',i2)
         stop
      endif

      call fkinspgmr20(maxl, maxlrst, msbpre)

      do 20 i = 1 , PROBSIZE
         uu(i) = 2.*i
         scale(i) = 1.
         constr(i) = 0.
 20   continue
      call fkinsol(neq, uu, 0, scale, scale, fnormtol, 
     1              scsteptol, constr, inopt, iopt, ropt, ier)

       write(6,1240)
 1240  format('kindiagsf example case. This kinsol example code does a')
       write(6,1241)
 1241  format('128 eqn diagonal algebraic system. Its purpose is to ')
       write(6,1242)
 1242  format('demonstrate the use of the Fortran interface in a')
       write(6,1243)
 1243  format('serial environment.')
       write(6,1244)
 1244  format('globalstrategy = INEXACT_NEWTON')

       write(6,1245)ier
 1245 format(1x,' return code is ',i5)


      write(6,1246)
 1246 format(1x,'The resultant values of uu are:',///)

      do 30 i = 1, PROBSIZE, 4
         write(6,1256)i, uu(i), uu(i+1), uu(i+2), uu(i+3)
 1256    format(i4,4(1x,f10.6))
 30   continue


       write(6,1267) iopt(4), iopt(11), iopt(5), 
     1    iopt(12), iopt(13),  iopt(14)
 1267 format(1x,'nni=',i4,' nli=',i4,' nfe=',i4,' npe=',i4,
     1       ' nps=',i4,' ncfl=',i4)

      call fkinfree

      stop

      end

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The function defining the system f(u) = 0. must be defined by a Fortran
c     function of the following form. 

      subroutine KFUN(nloc, uu, fval)
c
      double precision fval(*), uu(*)

      parameter(PROBSIZE=128)


      do 10 i = 1 , PROBSIZE
         fval(i) = uu(i)*uu(i) - i*i
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
      parameter(PROBSIZE=128)
      double precision pp
      common /pcom/ pp(PROBSIZE)
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do 10 i = 1 , PROBSIZE
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
      parameter(PROBSIZE=128)
      double precision pp
      common /pcom/ pp(PROBSIZE)
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do 10 i = 1 , PROBSIZE
         vtem(i) = vtem(i) * pp(i)
 10   continue

      return
      end
