sinclude(autoconf_macros/casc_fortran.m4)
sinclude(autoconf_macros/casc_libs_and_headers.m4)
sinclude(autoconf_macros/casc_mpi.m4)
sinclude(/home/casc/software/autoconfig/macros/casc_cxx.m4)
sinclude(/home/casc/software/autoconfig/macros/casc_lang.m4)
sinclude(/home/casc/software/autoconfig/macros/casc_misc.m4)
sinclude(/home/casc/software/autoconfig/macros/casc_opt_debug.m4)
sinclude(/home/casc/software/autoconfig/macros/casc_specific_libs.m4)

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
  AIX | rs6000)
	CASC_FIND_MPI
  ;;
  sun4 | solaris)
	CASC_FIND_MPI
  ;;
  LINUX)
        CASC_FIND_MPI
  ;;
  OSF1 | alpha)
	CASC_FIND_MPI
  ;;		
esac
])

dnl ***********************************************************************

AC_DEFUN(SUNDIALS_FIND_MPIF77,
[
   AC_ARG_WITH(MPIF77,
       AC_HELP_STRING([--with-MPIF77=ARG], [manually set MPIF77 to ARG]),
       [AC_MSG_RESULT([setting MPIF77 to $withval]); MPIF77=$withval],
       [CASC_PROG_MPIF77])
])

dnl that's all