#---------------------------------------------------------------------------------------------
# 
#
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CC,
[
  case $target in

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      AC_MSG_CHECKING([for C compiler])
      CC=xlc
      AC_MSG_RESULT([$CC])

    ;;

    # All others
    *)

      AC_PROG_CC(cc gcc)

    ;;

  esac
])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CXX,
[
  case $target in

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      AC_MSG_CHECKING([for C++ compiler])
      CXX=xlCC
      AC_MSG_RESULT([$CXX])

    ;;

    # All others
    *)

      AC_PROG_CXX(CC g++ gcc c++ cxx)

    ;;

  esac
])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_F77,
[
  case $target in

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      AC_MSG_CHECKING([for Fortran compiler])
      F77=xlf
      AC_MSG_RESULT([$F77])

    ;;

    # All others

    *)

      AC_PROG_F77(f77 g77 f90 xlf90)

    ;;

  esac
])

#---------------------------------------------------------------------------------------------
#
#
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CFLAGS,
[

  case $target in 

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      CFLAGS="${CFLAGS} -qmaxmem=18432"

    ;;

    # Linux
    i686-pc-linux-gnu)

      if test "X${GCC}" = "Xyes"; then       
        CFLAGS="${CFLAGS} -ffloat-store"
      fi

    ;;

  esac

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CXXFLAGS,
[
  case $target in 

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      CXXFLAGS_MISC="-qnolm -qrtti"

    ;;

    # Linux
    i686-pc-linux-gnu)

      if test "X${GCC}" = "Xyes"; then       
        CXXFLAGS="${CXXFLAGS} -ffloat-store"
      fi

    ;;
  esac
])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_FFLAGS,
[
  case $target in

    # SGI/IRIX
    mips-sgi-irix* ) 

      FFLAGS_MISC="-64"

    ;;

  esac
])

#---------------------------------------------------------------------------------------------
# 
#
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_BUILD_MODULES_LIST,
[

if test "X${EXAMPLES}" = "Xyes"; then
 EX_MODULES=
else
 EX_MODULES=no
fi
SUNDIALS_MAKEFILES=Makefile

# SHARED module

MODULES=shared/source
SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} shared/source/Makefile"

# NVECTOR modules

NVEC_MODULES=

if test -d ${srcdir}/nvec_ser; then
  NVEC_MODULES="$NVEC_MODULES nvec_ser"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} nvec_ser/Makefile"
fi

if test -d ${srcdir}/nvec_par && test "X${MPICC_OK}" = "Xyes" ; then
  NVEC_MODULES="$NVEC_MODULES nvec_par"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} nvec_par/Makefile"
fi

# CVODE module

if test -d ${srcdir}/cvode; then

  MODULES="$MODULES cvode/source cvode/fcmix"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/source/Makefile cvode/fcmix/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvode/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvode/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/examples_par/Makefile"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvode/fcmix/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvode/fcmix/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/examples_par/Makefile"
  fi

fi

# CVODES module

if test -d ${srcdir}/cvodes; then

  MODULES="$MODULES cvodes/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvodes/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvodes/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/examples_par/Makefile"
  fi

fi

# IDA module

if test -d ${srcdir}/ida; then

  MODULES="$MODULES ida/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES ida/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES ida/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/examples_par/Makefile"
  fi

fi

# IDAS module

if test -d ${srcdir}/idas; then

  MODULES="$MODULES idas/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES idas/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES idas/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/examples_par/Makefile"
  fi

fi

# KINSOL module

if test -d ${srcdir}/kinsol; then

  MODULES="$MODULES kinsol/source kinsol/fcmix"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/source/Makefile kinsol/fcmix/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES kinsol/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES kinsol/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/examples_par/Makefile"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES kinsol/fcmix/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES kinsol/fcmix/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/examples_par/Makefile"
  fi

fi

])

#---------------------------------------------------------------------------------------------
# 
#
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_MPI_SPECIFY,
[

  # MPI compilers
  # -------------

  AC_ARG_WITH(mpi-comp,
  [AC_HELP_STRING([--with-mpi-comp=DIR], [use MPI compilers mpicc, mpif77, and mpicxx 
                                          (or mpiCC) in the specified path.])],
  [
    if test X${withval} != Xno; then
      USE_MPI_COMPILERS=yes
      if test X${withval} = Xyes; then
        MPI_CC=mpicc
        MPI_F77=mpif77
        MPI_CXX=mpiCC
      else
        MPI_CC=${withval}/mpicc
        MPI_F77=${withval}/mpif77
        MPI_TEMP_CXX=${withval}/mpicxx
        if test -f ${MPI_TEMP_CXX}; then
          MPI_CXX=${MPI_TEMP_CXX}
        else
          MPI_CXX=${withval}/mpiCC
        fi
      fi
    else
      USE_MPI_COMPILERS=no
    fi
  ],
  [
    USE_MPI_COMPILERS=yes
    MPI_CC=mpicc
    MPI_F77=mpif77
    MPI_CXX=mpiCC
  ])

  # MPI root directory
  # ------------------

  AC_ARG_WITH(mpi,
  [AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory])],
  [
    USE_MPI_COMPILERS=no
    MPI_DIR=${withval}
    AC_MSG_CHECKING(MPI directory)
    AC_MSG_RESULT([${MPI_DIR}])
  ])

  # MPI include directory
  # ---------------------

  AC_ARG_WITH(mpi-incdir,
  [AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@])],
  [
    USE_MPI_COMPILERS=no
    MPI_INC=${withval}
    AC_MSG_CHECKING(user-defined MPI includes)
    AC_MSG_RESULT([${MPI_INC}])
  ])

  # MPI library directory
  # ---------------------

  AC_ARG_WITH(mpi-libdir,
  [AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@])],
  [
    USE_MPI_COMPILERS=no
    MPI_LIBDIR=${withval}
    AC_MSG_CHECKING(user-defined MPI library directory)
    AC_MSG_RESULT([${MPI_LIBDIR}])
  ])

  # MPI libraries
  # -------------

  AC_ARG_WITH(mpi-libs,
  [AC_HELP_STRING([--with-mpi-libs],[MPI libraries @<:@"-lmpi"@:>@])],
  [ 
    USE_MPI_COMPILERS=no
    MPI_LIBS=${withval}
    AC_MSG_CHECKING(user-defined MPI libraries)
    AC_MSG_RESULT([${MPI_LIBS}])
  ])

])

#---------------------------------------------------------------------------------------------
# 
#
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPICC,
[

  # Test the MPICC compiler

  if test -n "${MPI_CC}"; then
    if test -f ${MPI_CC}; then
      MPI_CC_EXISTS=yes
    else
      AC_CHECK_PROG(MPI_CC_EXISTS, ${MPI_CC}, yes, no)
    fi

    if test "X${MPI_CC_EXISTS}" = "Xyes"; then
      MPICC_OK=yes
      MPICC=${MPI_CC}
    else
      MPICC_OK=no
      echo "-----"
      echo "Cannot find MPI C compiler ${MPI_CC}."
      echo "Specify a path to all mpi compilers with --with-mpi-comp=DIR"
      echo "or specify a C compiler using CC=<compiler>"
      echo "-----"
      AC_MSG_ERROR([MPI C compiler (${MPI_CC}) not found.])
    fi
  fi

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPICXX,
[

  # Test the MPICXX compiler

  if test -n "${MPI_CXX}"; then

    if test -f ${MPI_CXX}; then   
      MPI_CXX_EXISTS=yes
    else
      AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
    fi

    if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
      MPICXX_OK=yes
      MPICXX=${MPI_CXX}
    else
      MPICXX_OK=no
      echo "-----"
      echo "Cannot find MPI C++ compiler ${MPI_CXX}."
      echo "Specify with --with-mpi-cxx."
      echo "-----"
      AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
    fi

  fi

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPIF77,
[

  # Test the MPIF77 compiler

  if test -n "${MPI_F77}"; then
    if test -f ${MPI_F77}; then
      MPI_F77_EXISTS=yes
    else
      AC_CHECK_PROG(MPI_F77_EXISTS, ${MPI_F77}, yes, no)
    fi

    if test "X${MPI_F77_EXISTS}" = "Xyes"; then
      MPIF77_OK=yes
      MPIF77=${MPI_F77}
    else
      MPIF77_OK=no
      echo "-----"
      echo "Cannot find MPI Fortran compiler ${MPI_F77}."
      echo "Specify a path to all mpi compilers with --with-mpi-comp=DIR"
      echo "or specify a Fortran 77 compiler using F77=<compiler>"
        echo "-----"
      AC_MSG_ERROR([MPI Fortran 77 compiler (${MPI_F77}) not found.])
    fi
  fi

])

#---------------------------------------------------------------------------------------------
# 
#
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_CC_WITH_MPI,
[

  # Save copies of current CPPFLAGS, LDFLAGS, LIBS

  SAVED_CPPFLAGS=${CPPFLAGS}
  SAVED_LDFLAGS=${LDFLAGS}
  SAVED_LIBS=${LIBS}

  # Try a simple MPI C program

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    MPIINC="-I${MPI_INC}"
    CPPFLAGS="${CPPFLAGS} ${MPIINC}"
  fi

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    MPILIBDIR="-L${MPI_LIBDIR}"
    LDFLAGS="${LDFLAGS} ${MPILIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    MPILIBS=${MPI_LIBS}
    LIBS="${LIBS} ${MPILIBS}"
  fi

  AC_LANG([C]) 
  AC_MSG_CHECKING(whether MPI will link using the C compiler)
  AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[#include "mpi.h"]], [[int c; char** v; MPI_Init(&c,&v);]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)  
   echo "-----"
   echo "Cannot link simple MPI program."
   echo "Try --with-mpi-comp to specify MPI C compile script."
   echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
   echo "to specify all the specific MPI compile options."
   echo "-----"
   AC_MSG_ERROR(MPI cannot link)
  ])

  # Restore CPPFLAGS, LDFLAGS, LIBS

  CPPFLAGS=${SAVED_CPPFLAGS}
  LDFLAGS=${SAVED_LDFLAGS}
  LIBS=${SAVED_LIBS}

  # Everything worked fine
  MPICC_OK=yes
])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_CXX_WITH_MPI,
[

  # Save copies of current CPPFLAGS, LDFLAGS, LIBS

  SAVED_CPPFLAGS=${CPPFLAGS}
  SAVED_LDFLAGS=${LDFLAGS}
  SAVED_LIBS=${LIBS}

  # Test a simple MPI C++ program

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    MPIINC="-I${MPI_INC}"
    CPPFLAGS="${CPPFLAGS} ${MPIINC}"
  fi

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    MPILIBDIR="-L${MPI_LIBDIR}"
    LDFLAGS="${LDFLAGS} ${MPILIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    MPILIBS=${MPI_LIBS}
    LIBS="${LIBS} ${MPILIBS}"
  fi

  AC_LANG([C++]) 
  AC_MSG_CHECKING(whether MPI will link using the C++ compiler)
  AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[#include "mpi.h"]], [[int c; char** v; MPI_Init(&c,&v);]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)  
   echo "-----"
   echo "Cannot link simple MPI program."
   echo "Try --with-mpi-comp to specify MPI C++ compile script."
   echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
   echo "to specify all the specific MPI compile options."
   echo "-----"
   AC_MSG_ERROR(MPI cannot link)
  ])

  # Restore CPPFLAGS, LDFLAGS, LIBS

  CPPFLAGS=${SAVED_CPPFLAGS}
  LDFLAGS=${SAVED_LDFLAGS}
  LIBS=${SAVED_LIBS}

  # Everything worked fine
  MPICXX_OK=yes

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_F77_WITH_MPI,
[

  # Save copies of current FFLAGS, LDFLAGS, LIBS

  SAVED_FFLAGS=${FFLAGS}
  SAVED_LDFLAGS=${LDFLAGS}
  SAVED_LIBS=${LIBS}

  # Try a simple MPI F77 program

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    MPIINC="-I${MPI_INC}"
    FFLAGS="${FFLAGS} ${MPIINC}"
  fi

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    MPILIBDIR="-L${MPI_LIBDIR}"
    LDFLAGS="${LDFLAGS} ${MPILIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    MPILIBS=${MPI_LIBS}
    LIBS="${LIBS} ${MPILIBS}"
  fi

  AC_LANG([Fortran 77]) 
  AC_MSG_CHECKING(whether MPI will link using the Fortran compiler)
  AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], 
  [      
      INCLUDE "mpif.h"
      CALL MPI_INIT(IER)
  ])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)  
   echo "-----"
   echo "Cannot link simple MPI program."
   echo "Try --with-mpi-comp to specify MPI F77 compile script."
   echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
   echo "to specify all the specific MPI compile options."
   echo "-----"
   AC_MSG_ERROR(MPI cannot link)
  ])

  # Restore FFLAGS, LDFLAGS, LIBS

  FFLAGS=${SAVED_FFLAGS}
  LDFLAGS=${SAVED_LDFLAGS}
  LIBS=${SAVED_LIBS}

  # Everything worked fine
  MPIF77_OK=yes

])

