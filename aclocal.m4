dnl sinclude(autoconf_macros/casc_fortran.m4)
sinclude(autoconf_macros/casc_libs_and_headers.m4)
sinclude(autoconf_macros/casc_mpi.m4)
sinclude(autoconf_macros/casc_misc.m4)

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_CC,
[
   case $ARCH in
      OSF1 | alpha)
        AC_MSG_CHECKING([for C compiler])
        CC=cc
        AC_MSG_RESULT([$CC])
      ;;
      AIX | rs6000)
        AC_MSG_CHECKING([for C compiler])
        CC=xlc
        AC_MSG_RESULT([$CC])
      ;;
      *)
        AC_PROG_CC
      ;;
   esac
]) 

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_F77,
[
   case $ARCH in
      OSF1 | alpha)
        AC_MSG_CHECKING([for Fortran compiler])
        F77=f77
        AC_MSG_RESULT([$F77])
      ;;
      AIX | rs6000)
        AC_MSG_CHECKING([for Fortran compiler])
        F77=xlf
        AC_MSG_RESULT([$F77])
      ;;
      *)
        AC_PROG_F77
      ;;
   esac
]) 

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_CFLAGS,
[
case $debug in
   yes)
      CFLAGS="-g"
   ;;
   no)
      CFLAGS="-O"
   ;;
   *)
      CFLAGS=$debug
   ;;
esac

case $ARCH in 
   AIX | rs6000)
	CFLAGS="$CFLAGS -qmaxmem=18432"
   ;;
   LINUX)
	CFLAGS="$CFLAGS -ffloat-store"
   ;;
esac

AC_MSG_CHECKING(for C compiler flags)
AC_MSG_RESULT($CFLAGS)
])

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_PROG_FFLAGS,
[
case $debug in
   yes)
      FFLAGS="-g"
   ;;
   no)
      FFLAGS="-O"
   ;;
   *)
      FFLAGS=$debug
   ;;
esac

AC_MSG_CHECKING(for Fortran compiler flags)
AC_MSG_RESULT($FFLAGS)
])

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_FIND_MPI,
[
case $ARCH in
  OSF1 | alpha)
	CASC_FIND_MPI
  ;;		
  *)
        CASC_FIND_MPI
  ;;
esac
])

dnl that's all