#---------------------------------------------------------------------------------------------
# 
#
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CC,
[
  case $target in

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      MY_CC=xlc

    ;;

  esac
])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CXX,
[
  case $target in

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      MY_CXX=xlCC

    ;;

  esac
])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_F77,
[
  case $target in

    # IBM
    rs6000-ibm-aix3.2.* | rs6000-ibm-aix4.* | powerpc-ibm-aix4.*)

      MY_F77=xlf

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

      MY_CFLAGS="${MY_CFLAGS} -qmaxmem=18432"

    ;;

    # Linux
    i686-pc-linux-gnu)

      if test "X${GCC}" = "Xyes"; then       
        MY_CFLAGS="${MY_CFLAGS} -ffloat-store"
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

      MY_CXXFLAGS="${MY_CXXFLAGS} -qnolm -qrtti"

    ;;

    # Linux
    i686-pc-linux-gnu)

      if test "X${GXX}" = "Xyes"; then       
        MY_CXXFLAGS="${MY_CXXFLAGS} -ffloat-store"
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

      MY_FFLAGS="-64"

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
if test -x ${MPI_CC}; then
  AC_MSG_CHECKING([for ${MPI_CC}])
  AC_MSG_RESULT([yes])
  MPI_CC_EXISTS=yes
else
  AC_CHECK_PROG(MPI_CC_EXISTS, ${MPI_CC}, yes, no)
fi

if test "X${MPI_CC_EXISTS}" = "Xyes"; then
  MPICC_OK=yes
  MPICC=${MPI_CC}
else
  MPICC_OK=no
  AC_MSG_WARN([MPI C compiler (${MPI_CC}) not found.])
  echo "
       Cannot find a working MPI C compiler.
       Try --with-mpi-comp to specify MPI C compiler script.
       Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
       to specify all the specific MPI compile options.

       Some modules will not be configured...
       "
fi

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPICXX,
[

# Test the MPICXX compiler

if test -x ${MPI_CXX}; then   
  AC_MSG_CHECKING([for ${MPI_CXX}])
  AC_MSG_RESULT([yes])
  MPI_CXX_EXISTS=yes
else
  AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
fi

if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
  MPICXX_OK=yes
  MPICXX=${MPI_CXX}
else
  MPICXX_OK=no
  AC_MSG_WARN([MPI C++ compiler (${MPI_CXX}) not found.])
  echo "
       Cannot find a working MPI C++ compiler.
       Try --with-mpi-comp to specify MPI C++ compiler script.
       Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
       to specify all the specific MPI compile options.

       Some modules will not be configured...
       "
fi

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPIF77,
[

# Test the MPIF77 compiler

if test -x ${MPI_F77}; then
  AC_MSG_CHECKING([for ${MPI_F77}])
  AC_MSG_RESULT([yes])
  MPI_F77_EXISTS=yes
else
  AC_CHECK_PROG(MPI_F77_EXISTS, ${MPI_F77}, yes, no)
fi

if test "X${MPI_F77_EXISTS}" = "Xyes"; then
  MPIF77_OK=yes
  MPIF77=${MPI_F77}
else
  MPIF77_OK=no
  AC_MSG_WARN([MPI F77 compiler (${MPI_F77}) not found.])
  echo "
       Cannot find a working MPI F77 compiler.
       Try --with-mpi-comp to specify MPI F77 compiler script.
       Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
       to specify all the specific MPI compile options.

       Some modules will not be configured...
       "
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
[AC_MSG_RESULT(yes)
    MPICC_OK=yes 
],
[AC_MSG_RESULT(no)
    MPICC_OK=no
    AC_MSG_WARN([MPI cannot link with ${CC}])
    echo "
         Cannot build a simple MPI program.
         Try --with-mpi-comp to specify MPI C compiler script.
         Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
         to specify all the specific MPI compile options.

         Some modules will not be configured...
         "
])

# Restore CPPFLAGS, LDFLAGS, LIBS

CPPFLAGS=${SAVED_CPPFLAGS}
LDFLAGS=${SAVED_LDFLAGS}
LIBS=${SAVED_LIBS}

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
[AC_MSG_RESULT(yes)
    MPICXX_OK=yes
],
[AC_MSG_RESULT(no)
    MPICXX_OK=no  
    AC_MSG_WARN([MPI cannot link with ${CXX}])
    echo "
         Cannot build a simple MPI program.
         Try --with-mpi-comp to specify MPI C++ compiler script.
         Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
         to specify all the specific MPI compile options.

         Some modules will not be configured...
         "
])

# Restore CPPFLAGS, LDFLAGS, LIBS

CPPFLAGS=${SAVED_CPPFLAGS}
LDFLAGS=${SAVED_LDFLAGS}
LIBS=${SAVED_LIBS}

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
[AC_MSG_RESULT(yes)
    MPIF77_OK=yes
],
[AC_MSG_RESULT(no)
    MPIF77_OK=no  
    AC_MSG_WARN([MPI cannot link with ${F77}])
    echo "
         Cannot build a simple MPI program.
         Try --with-mpi-comp to specify MPI F77 compiler script.
         Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
         to specify all the specific MPI compile options.

         Some modules will not be configured...
         "
])

# Restore FFLAGS, LDFLAGS, LIBS

FFLAGS=${SAVED_FFLAGS}
LDFLAGS=${SAVED_LDFLAGS}
LIBS=${SAVED_LIBS}

])

