      program kindiagsf
c   ***************************************************************************
c   * File        : kindiagsf.f                                               *
c   * Programmers : Allan G. Taylor and Alan C. Hindmarsh @ LLNL              *
c   * Version of  : 30 July 2002                                              *
c   *   Simple diagonal test with Fortran interface, using user-supplied      *
c   *   preconditioner setup and solve routines (supplied in Fortran, below). *
c   *   This example does a basic test of the solver by solving the system    *
c   *                        f(u) = 0   for                                   *
c   *                        f(u) = u(i)^2 - i^2 .                            *
c   *   No scaling is done.                                                   *
c   *   An approximate diagonal preconditioner is used.                       *
c   *   Execute line: kindiagsf                                               *
c   *-------------------------------------------------------------------------*
c   * Modified by Radu Serban to work with new serial NVECTOR, 12 March 2002. *
c   *-------------------------------------------------------------------------*

      integer PROBSIZE
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

c * * * * * * * * * * * * * * * * * * * * * *

      call fmenvinits(neq, ier)
      if (ier .ne. 0) then
         write(6,1220),ier
 1220    format('fmenvinits failed, ier =',i2)
         stop
      endif

      call fkinmalloc(neq, ier)
      if (ier .ne. 0) then
         write(6,1230),ier
 1230    format('fkinmalloc failed, ier =',i2)
         stop
      endif

      call fkinspgmr20(maxl, maxlrst, msbpre, ier)

      if (ier .ne. 0) then
         write(6,1232)ier
 1232    format('fkinspgmr20 failed, ier =',i2)
         stop
      endif

      do 20 i = 1,neq
         uu(i) = 2.0*i
         scale(i) = 1.0
         constr(i) = 0.0
  20  continue

      write(6,1240)
 1240 format('Example program kindiagsf'/' This fkinsol example code',
     1       ' solves a 128 eqn diagonal algebraic system.'/
     2       ' Its purpose is to demonstrate the use of the Fortran',
     3       ' interface'/' in a serial environment.'/
     4       ' globalstrategy = INEXACT_NEWTON'/)

      call fkinsol(neq, uu, 0, scale, scale, fnormtol, 
     1             scsteptol, constr, inopt, iopt, ropt, ier)

      write(6,1245)ier
 1245 format(/' fkinsol return code is ',i5)

      write(6,1246)
 1246 format(/' The resultant values of uu are:'/)

      do 30 i = 1,neq,4
         write(6,1256)i, uu(i), uu(i+1), uu(i+2), uu(i+3)
 1256    format(i4,4(1x,f10.6))
 30   continue

      write(6,1267)iopt(4),iopt(11),iopt(5),iopt(12),iopt(13),iopt(14)
 1267 format(//' nni=',i4,',  nli=',i4,',  nfe=',i4,',  npe=',i4,
     1       ',  nps=',i4,',  ncfl=',i4)

      call fkinfree
      call fmenvfrees

      stop
      end
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The function defining the system f(u) = 0 must be defined by a Fortran
c     function of the following form.
      
      subroutine KFUN(neq, uu, fval)

      double precision fval(*), uu(*)

      do 10 i = 1,neq
 10      fval(i) = uu(i)*uu(i) - i*i
      
      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpreco is the preconditioner setup routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine kpreco(neq, udata, uscale, fdata, fscale, 
     1                  vtemp1, vtemp2, uround, nfe, ier)
      
      integer nfe, ier
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*), uround
      
      double precision pp
      common /pcom/ pp(128)
      
      do 10 i = 1,neq
 10      pp(i) = 0.5/(udata(i)+5.)
      
      ier = 0
      
      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpsol is the preconditioner solve routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine kpsol(neq, udata, uscale, fdata, fscale, 
     1                 vv, ftem, uround, nfe, ier)
      
      integer nfe, ier
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vv(*), ftem(*), uround
      
      double precision pp
      common /pcom/ pp(128)
      
      do 10 i = 1,neq
 10      vv(i) = vv(i) * pp(i)
      
      ier = 0
      
      return
      end
