      program kindiagsf
c
c     File        : kindiagsf.f                                             
c     Programmers : Allan G. Taylor, Alan C. Hindmarsh, Radu Serban @ LLNL  
c     Version of  : 27 January 2004
c     Simple diagonal test with Fortran interface, using user-supplied      
c     preconditioner setup and solve routines (supplied in Fortran, below). 
c     This example does a basic test of the solver by solving the system    
c          f(u) = 0   for                                   
c          f(u) = u(i)^2 - i^2 .                            
c     No scaling is done.                                                   
c     An approximate diagonal preconditioner is used.                       
c     Execute line: kindiagsf                                               
c   

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
      common /psize/ neq

      neq = PROBSIZE
      globalstrat = 0
      fnormtol = 0.00001
      scsteptol = 0.0001
      inopt = 0
      maxl = 10
      maxlrst = 2
      msbpre  = 5

c * * * * * * * * * * * * * * * * * * * * * *

      call fnvinits(neq, ier)
      if (ier .ne. 0) then
         write(6,1220) ier
 1220    format('fnvspecinits failed, ier =',i2)
         stop
      endif

      do 20 i = 1,neq
         uu(i) = 2.0*i
         scale(i) = 1.0
         constr(i) = 0.0
  20  continue

      call fkinmalloc(msbpre, fnormtol, scsteptol, 
     &                constr, inopt, iopt, ropt, ier)
      if (ier .ne. 0) then
         write(6,1230) ier
 1230    format('fkinmalloc failed, ier =',i2)
         stop
      endif

      call fkinspgmr(maxl, maxlrst, ier)
      call fkinspgmrsetpsol(1, ier)
      call fkinspgmrsetpset(1, ier)

      write(6,1240)
 1240 format('Example program kindiagsf'/' This fkinsol example code',
     1       ' solves a 128 eqn diagonal algebraic system.'/
     2       ' Its purpose is to demonstrate the use of the Fortran',
     3       ' interface'/' in a serial environment.'/
     4       ' globalstrategy = INEXACT_NEWTON'/)

      call fkinsol(uu, 0, scale, scale, ier)

      write(6,1245) ier
 1245 format(/' fkinsol return code is ',i5)

      write(6,1246)
 1246 format(/' The resultant values of uu are:'/)

      do 30 i = 1,neq,4
         write(6,1256) i, uu(i), uu(i+1), uu(i+2), uu(i+3)
 1256    format(i4,4(1x,f10.6))
 30   continue

      write(6,1267) iopt(4),iopt(11),iopt(5),iopt(12),iopt(13),iopt(14)
 1267 format(//' nni=',i4,',  nli=',i4,',  nfe=',i4,',  npe=',i4,
     1       ',  nps=',i4,',  ncfl=',i4)

      call fkinfree
      call fnvfrees

      stop
      end
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The function defining the system f(u) = 0 must be defined by a Fortran
c     function of the following form.
      
      subroutine fkfun(uu, fval)
      
      double precision fval(*), uu(*)
      integer neq
      common /psize/ neq

      do 10 i = 1,neq
         fval(i) = uu(i)*uu(i) - i*i
 10   continue
      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpreco is the preconditioner setup routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpset(udata, uscale, fdata, fscale, 
     1                  vtemp1, vtemp2, ier)
      
      integer ier
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)
      
      double precision pp
      common /pcom/ pp(128)
      
      integer neq
      common /psize/ neq

      do 10 i = 1,neq
         pp(i) = 0.5/(udata(i)+5.)
 10   continue
      ier = 0
      
      return
      end
      
      
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     The routine kpsol is the preconditioner solve routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:
      
      subroutine fkpsol(udata, uscale, fdata, fscale, 
     1                  vv, ftem, ier)
      
      integer ier
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vv(*), ftem(*)
      
      double precision pp
      common /pcom/ pp(128)

      integer neq
      common /psize/ neq

      do 10 i = 1,neq
         vv(i) = vv(i) * pp(i)
 10   continue
      ier = 0
      
      return
      end
