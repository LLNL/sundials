c     ----------------------------------------------------------------
c     $Revision: 1.2 $
c     $Date: 2005-06-21 19:12:56 $
c     ----------------------------------------------------------------
c     This simple example problem for FIDA, due to Robertson, 
c     is from chemical kinetics, and consists of the following three 
c     equations:
c
c          dy1/dt = -.04*y1 + 1.e4*y2*y3
c          dy2/dt =  .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c             0   = y1 + y2 + y3 - 1
c
c     on the interval from t = 0.0 to t = 4.e10, with initial
c     conditions: y1 = 1, y2 = y3 = 0.
c
c     The problem is solved using a dense linear solver, with a
c     user-supplied Jacobian. Output is printed at
c     t = .4, 4, 40, ..., 4e10.
c     ----------------------------------------------------------------
c
      program irobxf
c
      implicit none
c
      integer ier
      integer*4 iopt(40)
      double precision ropt(40)
c
      integer iatol, inopt, nout, iout, itask
      integer nst, kused
      integer*4 neq, i
      double precision t0, t1, rtol, tout, tret
      double precision hused
      double precision y(3), yp(3), atol(3), id(3), constr(3)
c
      data nst/15/, kused/19/, hused/16/
c
c Initialize variables
c
      neq = 3
      nout = 12
      rtol = 1.0d-4
      t0 = 0.0d0
      t1 = 0.4d0
      iatol = 2
      inopt = 0
      itask = 1
c
      y(1) = 1.0d0
      y(2) = 0.0d0
      y(3) = 0.0d0
c
      yp(1) = -0.04d0
      yp(2) = 0.04d0
      yp(3) = 0.0d0
c
      atol(1) = 1.0d-6
      atol(2) = 1.0d-10
      atol(3) = 1.0d-6
c
c Initialize IDA vector environment
c
      call fnvinits(2, neq, ier)
      if (ier .ne. 0) then
         write(6,10) ier
 10      format(///' SUNDIALS_ERROR: FNVINITS returned IER = ', i5)
         stop
      endif
c
      call fidamalloc(t0, y, yp, iatol, rtol, atol, id, constr, inopt,
     &                iopt, ropt, ier)
      if (ier .ne. 0) then
         write(6,20) ier
 20      format(///' SUNDIALS_ERROR: FIDAMALLOC returned IER = ', i5)
         stop
      endif
c
c Attach dense linear solver
c
      call fidadense(neq, ier)
      call fidadensesetjac(1, ier)
c
c Print header
c
      call prntintro(rtol, atol, y)
c
      tout = t1
      do 30 iout = 1, nout
c
         call fidasolve(tout, tret, y, yp, itask, ier)
c
         write(6,40) tret, (y(i), i = 1,3), iopt(nst), iopt(kused),
     &               ropt(hused)
 40       format(e8.2, 3(1x,e12.4), i5, i3, e12.4)
c
         if (ier .ne. 0) then
            write(6,50) ier
 50         format(///' SUNDIALS_ERROR: FIDASOLVE returned IER = ', i5)
            call fidafree
            stop
         endif
c
         tout = tout*10.0d0
 30   continue
c
c Print final statistics
c
      call prntstats(iopt)
c
c Free IDA memory
c
      call fidafree
c
      stop
      end
c
c ==========
c
      subroutine fidaresfun(tres, y, yp, res, reserr)
c
      implicit none
c
      integer reserr
      double precision tres
      double precision y(*), yp(*), res(*)
c
      res(1) = -0.04d0*y(1)+1.0d4*y(2)*y(3)
      res(2) = -res(1)-3.0d7*y(2)*y(2)-yp(2)
      res(1) = res(1)-yp(1)
      res(3) = y(1)+y(2)+y(3)-1.0d0
c
      return
      end
c
c ==========
c
      subroutine fidadjac(neq, t, y, yp, r, jac, cj, ewt, h,
     1                    wk1, wk2, wk3, djacerr)
c
      implicit none
c
      integer*4 neq
      integer djacerr
      double precision t, h, cj
      double precision y(*), yp(*), r(*), ewt(*), jac(neq,neq)
      double precision wk1(*), wk2(*), wk3(*)
c
      jac(1,1) = -0.04d0-cj
      jac(2,1) = 0.04d0
      jac(3,1) = 1.0d0
      jac(1,2) = 1.0d4*y(3)
      jac(2,2) = -1.0d4*y(3)-6.0d7*y(2)-cj
      jac(3,2) = 1.0d0
      jac(1,3) = 1.0d4*y(2)
      jac(2,3) = -1.0d4*y(2)
      jac(3,3) = 1.0d0
c
      return
      end
c
c ==========
c
      subroutine prntintro(rtol, atol, y)
c
      implicit none
c
      integer*4 i
      double precision rtol, atol(*), y(*)
c
      write(6,60) rtol, (atol(i), i = 1,3), (y(i), i = 1,3)
 60   format(/'irobxf: Robertson kinetics DAE serial example',
     &       'problem for IDA', /,'        Three equation chemical',
     &       'kinetics problem.', //,
     &       'Tolerance parameters:  rtol = ', e8.2,
     &       '   atol = ', 3(1x,e8.2), /,
     &       'Initial conditions y0 = (', 3(1x,e8.2), ')', /,
     &       'Constraints and id not used.', //,
     &       '  t          y1           y2           y3        nst',
     &       '  k    h')
c
      return
      end
c
c ==========
c
      subroutine prntstats(iopt)
c
      implicit none
c
      integer*4 iopt(40)
      integer nst, reseval, jaceval, nni, ncf, netf
c
      data nst/15/, reseval/16/, jaceval/28/, nni/21/, netf/18/,
     &     ncf/22/
c
      write(6,70) iopt(nst), iopt(reseval), iopt(jaceval),
     &             iopt(nni), iopt(netf), iopt(ncf)
 70   format(/'Final Run Statistics:', //,
     &         'Number of steps                    = ', i3, /,
     &         'Number of residual evaluations     = ', i3, /,
     &         'Number of Jacobian evaluations     = ', i3, /,
     &         'Number of nonlinear iterations     = ', i3, /,
     &         'Number of error test failures      = ', i3, /,
     &         'Number of nonlinear conv. failures = ', i3)
c
      return
      end
