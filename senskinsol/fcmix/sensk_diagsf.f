      program senskdiagsf
c   ************************************************************************
c   *                                                                      *
c   * File: sensk_diagsf.f                                                 *
c   * Programmers: Allan G. Taylor, Alan C. Hindmarsh, and                 *
c   *              Keith E. Grant @ LLNL                                   *
c   *                                                                      *
c   * Version of 06 Dec 2000                                               *
c   *   Modified as a FORTRAN interface example of parameter sensitivity   *
c   *   analysis. A parameter of nominal magnitude one is multiplied into  *
c   *   the 'i^2' term.                                                    *
c   *                                                                      *
c   * Version of 25 Mar 1999                                               *
c   *   simple diagonal test with Fortran interface, using user-supplied   *
c   *   preconditioner setup and solve routines (supplied in Fortran,      *
c   *   below) This example does a basic test of the code suite by         *
c   *   solving the system                                                 *
c   *                        f(u) = 0.  for                                *
c   *                        f(u) = u(i)^2 - p*(i^2) .                     *
c   *   No scaling is done.                                                *
c   *                                                                      *
c   *   execute line:  sensk_diagsf                                        *
c   *                                                                      *
c   ************************************************************************
c
      implicit none

      integer PROBSIZE, nopt
      parameter(PROBSIZE=128)
      parameter(nopt=40)

c
c.... Define to minimize nonsignificant code differences relative to
c.... the parallel version
      integer mype
      parameter (mype=0)

      integer neq, ier, globalstrat
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
      integer*4 iopt(nopt)
c     integer*8 iopt(nopt)
c **********************************************************************

      double precision fnormtol, scsteptol
      double precision ropt(nopt)
      double precision uu(PROBSIZE), scale(PROBSIZE)
      double precision constr(PROBSIZE)
      double precision ri
      double precision pp

      common /pcom/ pp(PROBSIZE)

c
c.... Define parameter and nominal magnitude variables for
c.... sensitiivty analysis
      integer np, idiffrhs
      double precision diffscl
      parameter (
     &   np = 1,
     &   idiffrhs = 0,
     &   diffscl = 1.0)
      double precision params(np), pbar(np), ww(PROBSIZE)

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

c
c.... Define the parameter and nominal magnitude values for
c.... sensitivity analysis
      do ip = 1,np
         params(ip) = 1.0
         pbar(ip)   = 1.0
      enddo

c
c.... Allocate and initialize SensKINSOL memory
      call fssenskinmalloc(neq, np, idiffrhs, diffscl, params,
     &   pbar, ier)
      if (ier .ne. 0) then
         write(6,'(a,i2/)') ' *** FSSENSKINMALLOC failed with ', ier
         stop
      endif

c
c.... Initialize the linear solver
      call fsenskinspgmr20(maxl, maxlrst, msbpre)

      if (mype .eq. 0) then
         write(6,'(//a/a/a/a//)')
     &      'SENS_DIAGSF demo case. This kinsol example code does',
     &      'a 128 eqn diagonal algebraic system. Its purpose is to',
     &      'demonstrate the use of the Fortran interface in a',
     &      'serial environment.'
         write(6, '(a//)') ' globalstrategy = INEXACT_NEWTON'
      endif

      do ii = 1,PROBSIZE
         ri = ii
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
     &      ' The resultant values of uu are:'

         do i = 1,PROBSIZE,4
            write(6,'(i4,4(1x,f10.6))') i,
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
            write (6,'(//a//)')
     &      'The Sensitivity values s are:'
            do i = 1,PROBSIZE,4
               write(6,'(i4,4(1x,f10.6))') i,
     &            ww(i)/pbar(ip), ww(i+1)/pbar(ip),
     &            ww(i+2)/pbar(ip), ww(i+3)/pbar(ip)
            enddo

         write(6,'(/6(1x,a,i4))')
     &      'nni=', iopt(4),  'nli=', iopt(11), 'nfe=',  iopt(5),
     &      'npe=', iopt(12), 'nps=', iopt(13), 'ncfl=', iopt(14)
         endif

      enddo


      call fsenskinfree

      stop
      end

c **************************************************************************
c     The function defining the system f(u) = 0. must be defined by a
c     FORTRAN function of the following form. To use this for a sensitivity
c     example, a single scaling parameter of nominal value 1.0 is
c     multiplied into the 'i^2' term.
c **************************************************************************

      subroutine KFUN(nloc, uu, fval, np, params)

      implicit none

      integer nloc, np
      double precision fval(*), uu(*), params(np)

      integer i

      integer PROBSIZE
      parameter(PROBSIZE=128)

      do i = 1,PROBSIZE
         fval(i) = uu(i)*uu(i) - params(1)*(i*i)
      enddo

      return
      end


c **************************************************************************
c   The routine kpreco is the preconditioner setup routine. It is required
c   that the specific name be used in order that the C code can find and
c   link to it. The argument list must also be as illustrated below:
c **************************************************************************

      subroutine kpreco(neq, udata, uscale, fdata, fscale,
     1                  vtemp1, vtemp2, uround, nfe, ier)

      integer neq, nfe, ier

      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)
      double precision uround

      integer i

c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      integer PROBSIZE
      parameter(PROBSIZE=128)
      double precision pp
      common /pcom/ pp(PROBSIZE)
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do i = 1,PROBSIZE
         pp(i) = 0.5/(udata(i)+5.)
      enddo

      return
      end


c **************************************************************************
c   The routine kpsol is the preconditioner solve routine. It is required
c   that the specific name be used in order that the C code can find and
c   link to it. The argument list must also be as illustrated below:
c **************************************************************************

      subroutine kpsol(neq, udata, uscale, fdata, fscale,
     &   vtem, ftem, uround, nfe, ier)

      implicit none

      integer neq, nfe, ier

      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtem(*), ftem(*)
      double precision uround

      integer i

c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *
      integer PROBSIZE
      parameter(PROBSIZE=128)
      double precision pp
      common /pcom/ pp(PROBSIZE)
c     common pcom * * * * * * * * * * * * * * * * * * * * * * * * * * *

      ier = 0

      do i = 1,PROBSIZE
         vtem(i) = vtem(i) * pp(i)
      enddo

      return
      end
