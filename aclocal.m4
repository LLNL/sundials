# ------------------------------------------------------------------------
# $Revision: 1.15 $
# $Date: 2004-05-04 21:48:23 $
# ------------------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# ------------------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California
# Produced at the Lawrence Livermore National Laboratory
# All rights reserved
# ------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
# Perform misc. initializations
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_INITIALIZE,
[

# Say Hi!
echo "----------------------------------"
echo "Running SUNDIALS Configure Script"
echo "----------------------------------"

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.
AC_CONFIG_SRCDIR(/shared/source/nvector.c)

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.
AC_CONFIG_AUX_DIR(config)

# If --host was not specified, guess it
if test -z "${host}"; then
   AC_CANONICAL_BUILD
   AC_CANONICAL_HOST
fi

# Only systems running Tru64 Unix should compile the CVODE and KINSOL fcmix dev examples.
case $host in

  # Compaq/Tru64
  *-*-osf* ) SUNDIALS_DEV_EXAMPLES="yes" ;;

  # other
  * ) SUNDIALS_DEV_EXAMPLES="no" ;;

esac

# Overwrite the default installation path
DEFAULT_PREFIX=`pwd`
AC_PREFIX_DEFAULT(`pwd`)

# If make predefines the Make variable MAKE, define output variable 
# SET_MAKE to be empty. Otherwise, define SET_MAKE to contain `MAKE=make'. 
AC_PROG_MAKE_SET

])

#---------------------------------------------------------------------------------------------
# Function to check for user overrides (CC, F77, CFLAGS, FFLAGS, CXX, CXXFLAGS)
# After the call to this function, MY_VAR is 'none' if VAR is not set, or has
# the value of VAR if VAR is set (in which case VAR is unset)
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_OVERRIDES,
[
MY_CC=none
MY_CXX=none
MY_F77=none
MY_CFLAGS=none
MY_CXXFLAGS=none
MY_FFLAGS=none

if test -n "${CC}"; then
  MY_CC=${CC}
  unset CC
fi

if test -n "${CXX}"; then
  MY_CXX=${CXX}
  unset CXX
fi

if test -n "${F77}"; then
  MY_F77=${F77}
  unset F77
fi

if test -n "${CFLAGS}"; then
  MY_CFLAGS=${CFLAGS}
  unset CFLAGS
fi

if test -n "${CXXFLAGS}"; then
  MY_CXXFLAGS=${CXXFLAGS}
  unset CXXFLAGS
fi

if test -n "${FFLAGS}"; then
  MY_FFLAGS=${FFLAGS}
  unset FFLAGS
fi

])

#---------------------------------------------------------------------------------------------
# Find an archiver
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_AR,
[
AC_CHECK_PROG(AR,ar,ar,no)
if test "X${AR}" = "Xno"; then
  AC_MSG_ERROR(cannot find an archiver)
else
  AR="${AR} rc"
fi
])

#---------------------------------------------------------------------------------------------
# Function to set C compilation
# Check if the user requested additional preprocessor and compiler flags
# and/or additional libraries
# The following are set here:
#   - CPPFLAGS
#   - CFLAGS
#   - CC
#   - LDFLAGS
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_SET_CC,
[
AC_ARG_WITH([],[],[])

EXTRA_CPPFLAGS=""
AC_ARG_WITH(single-prec,
[AC_HELP_STRING([--with-single-prec],[use single-precision arithmetic @<:@double prec@:>@])],
[
  if test "X${withval}" != "Xno"; then
    EXTRA_CPPFLAGS="-DSUNDIALS_SINGLE_PRECISION"
  fi
])

AC_ARG_WITH([],[ ],[])

AC_ARG_WITH(cppflags,
[AC_HELP_STRING([--with-cppflags],[add extra preprocessor flags])],
[EXTRA_CPPFLAGS="${EXTRA_CPPFLAGS} ${withval}"]
)

AC_ARG_WITH(cflags,
[AC_HELP_STRING([--with-cflags],[add extra C compiler flags])],
[EXTRA_CFLAGS=${withval}]
)

AC_ARG_WITH(ldflags,
[AC_HELP_STRING([--with-ldflags],[add extra linker flags])],
[EXTRA_LDFLAGS=${withval}]
)

AC_ARG_WITH(libs,
[AC_HELP_STRING([--with-libs],[add extra libraries])],
[EXTRA_LIBS=${withval}]
)
LIBS="${LIBS} ${EXTRA_LIBS}"

CC_OK=no

# Check for preprocessor flags
AC_MSG_CHECKING([for extra preprocessor flags])
if test -z "${EXTRA_CPPFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CPPFLAGS}])
  CPPFLAGS="${CPPFLAGS} ${EXTRA_CPPFLAGS}"
fi

# C compiler and flags
if test "X${MY_CC}" = "Xnone"; then
  MY_CC=cc
fi
AC_PROG_CC(${MY_CC} cc gcc)
test -n "${CC}" && CC_OK=yes

AC_MSG_CHECKING([for C compiler flags])
if test "X${MY_CFLAGS}" = "Xnone"; then
  MY_CFLAGS="-O"
  SUNDIALS_PROG_CFLAGS
fi
CFLAGS=${MY_CFLAGS}
AC_MSG_RESULT([${CFLAGS}])

AC_MSG_CHECKING([for extra C compiler flags])
if test -z "${EXTRA_CFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CFLAGS}])
  CFLAGS="${CFLAGS} ${EXTRA_CFLAGS}"
fi

AC_LANG([C])
AC_HEADER_STDC

if test "X${EXAMPLES}" = "Xyes"; then

  AC_LANG([C])
  AC_CHECK_HEADERS(
    [math.h],
    break,
    AC_MSG_WARN([You do not appear to have math.h!])
  )

  AC_LANG([C])
  AC_SEARCH_LIBS(pow,[m],,AC_MSG_WARN([Cannot find math library]))
  AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_WARN([Cannot find math library]))

fi

# Linker flags
AC_MSG_CHECKING([for extra linker flags])
if test -z "${EXTRA_LDFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_LDFLAGS}])
  LDFLAGS="${LDFLAGS} ${EXTRA_LDFLAGS}"
fi


])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CFLAGS,
[

  case $host in 

    # Linux
    i*-pc-linux-gnu)

      if test "X${GCC}" = "Xyes"; then       
        MY_CFLAGS="${MY_CFLAGS} -ffloat-store"
      fi

    ;;

  esac

])


#---------------------------------------------------------------------------------------------
# Function to set Fortran compilation
# It first tests whether a Fortran compiler is needed
# Unless the user explicitely disables it, Fortran support is enabled if
# either the cvode or the kinsol module is present
# The following are set here:
#   - F77
#   - FFLAGS
#   - FCFLAGS
#   - LDFLAGS 
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_SET_F77,
[
AC_ARG_WITH([],[  ],[])

F77_ENABLED=no
if test -d ${srcdir}/cvode || test -d ${srcdir}/kinsol ; then
  F77_ENABLED=yes
fi
AC_ARG_WITH([f77],
[AC_HELP_STRING([--without-f77], [skip F77 support])],
[
  if test "X${withval}" = "Xno"; then
    F77_ENABLED=no
  fi
])

AC_ARG_WITH(fflags,
[AC_HELP_STRING([--with-fflags],[add extra Fortran compiler flags])],
[EXTRA_FFLAGS=${withval}]
)

AC_ARG_WITH(f77underscore, 
[AC_HELP_STRING([--with-f77underscore[[[[=ARG]]]]],[specify number of underscores to append to function names (none/one/two) [ARG=one]],[])],
[
  if test "X${withval}" = "Xnone"; then
    FCFLAGS="-DSUNDIALS_UNDERSCORE_NONE"
  elif test "X${withval}" = "Xtwo"; then
    FCFLAGS="-DSUNDIALS_UNDERSCORE_TWO"
  fi
])

F77_OK=no

if test "X${F77_ENABLED}" = "Xyes"; then

  if test "X${MY_F77}" = "Xnone"; then
    MY_F77=f77
  fi
  AC_PROG_F77(${MY_F77} f77 g77)

  if test -n "${F77}"; then 

    F77_OK=yes

    AC_MSG_CHECKING([for Fortran compiler flags])
    if test "X${MY_FFLAGS}" = "Xnone"; then
      MY_FFLAGS="-O"
      SUNDIALS_PROG_FFLAGS
    fi
    FFLAGS=${MY_FFLAGS}
    AC_MSG_RESULT([${FFLAGS}])

    AC_MSG_CHECKING([for extra Fortran compiler flags])
    if test -z "${EXTRA_FFLAGS}"; then
      AC_MSG_RESULT(none)
    else
      AC_MSG_RESULT([${EXTRA_FFLAGS}])
      FFLAGS="${FFLAGS} ${EXTRA_FFLAGS}"
    fi

    AC_F77_LIBRARY_LDFLAGS

  else

    AC_MSG_WARN([no acceptable F77 compiler found])
    echo "
         Cannot find a working Fortran compiler.
         Some modules will not be configured...
         "
  fi

fi

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_FFLAGS,
[
  case $host in

    # Linux
    i*-pc-linux-gnu)

      if test "X${G77}" = "Xyes"; then       
        MY_FFLAGS="${MY_FFLAGS} -ffloat-store"
      fi

    ;;

    # SGI/IRIX
    mips-sgi-irix* ) 

      MY_FFLAGS="-64"

    ;;

    # Compaq/Tru64
    alpha*-dec-osf*)

      if test "X${F77}" = "Xf77"; then
        MY_FFLAGS="-O1"
      fi

    ;;

  esac
])

#---------------------------------------------------------------------------------------------
# Function to set C++ compilation
# The following are set here:
#   - CXXFLAGS
#   - CXX
# It first tests whether a C++ compiler is needed
# Unless the user explicitely disables it, C++ support is enabled if
# the xs4c module is present
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_SET_CXX,
[
CXX_ENABLED=no
test -d ${srcdir}/xs4c && CXX_ENABLED=yes
AC_ARG_WITH([cxx],
[AC_HELP_STRING([--without-cxx], [skip C++ support])],
[
  if test "X${withval}" = "Xno"; then
    CXX_ENABLED=no
  fi
])

AC_ARG_WITH(cxxflags,
[AC_HELP_STRING([--with-cxxflags],[add extra C++ compiler flags])],
[EXTRA_CXXFLAGS=${withval}]
)

CXX_OK=no

if test "X${CXX_ENABLED}" = "Xyes"; then

  if test "X${MY_CXX}" = "Xnone"; then
    MY_CXX=g++
  fi
  AC_PROG_CXX(${MY_CXX} g++ CC gcc c++ cxx)

  if test -n "${CXX}"; then 

    CXX_OK=yes

    AC_MSG_CHECKING([for C++ compiler flags])
    if test "X${MY_CXXFLAGS}" = "Xnone"; then
      MY_CXXFLAGS="-O"
      SUNDIALS_PROG_CXXFLAGS
    fi
    CXXFLAGS=${MY_CXXFLAGS}
    AC_MSG_RESULT([${CXXFLAGS}])

    AC_MSG_CHECKING([for extra C++ compiler flags])
    if test -z "${EXTRA_CXXFLAGS}"; then
      AC_MSG_RESULT(none)
    else
      AC_MSG_RESULT([${EXTRA_CXXFLAGS}])
      CXXFLAGS="${CXXFLAGS} ${EXTRA_CXXFLAGS}" 
    fi

    AC_LANG([C++])
    AC_CHECK_HEADERS(
      [complex],
      break,
      AC_MSG_NOTICE([You do not appear to have complex.h!])
    )

  else

    AC_MSG_WARN([no acceptable C++ compiler found])
    echo "
         Cannot find a working C++ compiler.
         Some modules will not be configured...
         "
  fi
  
fi

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_PROG_CXXFLAGS,
[
  case $host in 

    # Linux
    i686-pc-linux-gnu)

      if test "X${GXX}" = "Xyes"; then       
        MY_CXXFLAGS="${MY_CXXFLAGS} -ffloat-store"
      fi

    ;;
  esac
])

#---------------------------------------------------------------------------------------------
# Function to set MPI support for C
# Unless the user explicitely disables it, MPI support is enabled if nvec_par is present
# Decide first whether we'll use an MPI script (mpicc) and then test it
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_SET_MPICC,
[
AC_ARG_WITH([],[   ],[])

MPI_ENABLED=no
test -d ${srcdir}/nvec_par && MPI_ENABLED=yes
AC_ARG_WITH([mpi],
[AC_HELP_STRING([--without-mpi], [skip MPI parallel support])],
[
  if test "X${withval}" = "Xno"; then
    MPI_ENABLED=no
  fi
])

if test "X${MPI_ENABLED}" = "Xyes"; then

  # MPI root directory
  # ------------------

  AC_ARG_WITH(mpi-root,
  [AC_HELP_STRING([--with-mpi-root=MPIROOT],[use MPI root directory])],
  [
    MPI_DIR=${withval}
    AC_MSG_CHECKING(MPI directory)
    AC_MSG_RESULT([${MPI_DIR}])
  ])

  # MPI include directory
  # ---------------------

  AC_ARG_WITH(mpi-incdir,
  [AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@])],
  [
    MPI_INC=${withval}
    AC_MSG_CHECKING(user-defined MPI includes)
    AC_MSG_RESULT([${MPI_INC}])
  ])

  # MPI library directory
  # ---------------------

  AC_ARG_WITH(mpi-libdir,
  [AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@])],
  [
    MPI_LIBDIR=${withval}
    AC_MSG_CHECKING(user-defined MPI library directory)
    AC_MSG_RESULT([${MPI_LIBDIR}])
  ])

  # MPI libraries
  # -------------

  AC_ARG_WITH(mpi-libs,
  [AC_HELP_STRING([--with-mpi-libs],[MPI libraries @<:@"-lmpi"@:>@])],
  [ 
    MPI_LIBS=${withval}
    AC_MSG_CHECKING(user-defined MPI libraries)
    AC_MSG_RESULT([${MPI_LIBS}])
  ])

  # MPICC compiler
  # --------------

  AC_ARG_WITH(mpicc,
  [AC_HELP_STRING([--with-mpicc[[[[=ARG]]]]],[specify MPI C compiler to use @<:@mpicc@:>@],
                  [                                ])],
  [
    if [ test "X${withval}" != "Xno" ]; then
      USE_MPICC="yes"
      MPI_CC="${withval}"
    else
      USE_MPICC="no"
    fi
  ],
  [
    USE_MPICC="yes"
    MPI_CC="mpicc"
  ])

  # Test MPI
  # --------

  if test "X${USE_MPICC}" = "Xyes"; then
    SUNDIALS_CHECK_MPICC
  else
    MPICC=${CC}
    SUNDIALS_CHECK_CC_WITH_MPI
  fi


fi

])

#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPICC,
[

# Test the MPICC compiler

if test -x ${MPI_CC}; then
  MPI_CC_EXISTS=yes
else
  TEMP_MPI_CC=`basename "${MPI_CC}"`
  AC_CHECK_PROG(MPI_CC_EXISTS, ${TEMP_MPI_CC}, yes, no)
fi

if test "X${MPI_CC_EXISTS}" = "Xyes"; then
  MPICC_OK=yes
  MPICC=${MPI_CC}
else
  MPICC_OK=no
  AC_MSG_WARN([MPI C compiler (${MPI_CC}) not found.])
  echo "
       Cannot find a working MPI C compiler.
       Try --with-mpicc to specify MPI C compiler script.
       Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
       to specify all the specific MPI compile options.

       Some modules will not be configured...
       "
fi

])

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
         Try --with-mpicc to specify a MPI C compiler script.
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
# Function to set MPI support for Fortran
# Unless the user explicitely disables it, MPI support is enabled if nvec_par is present
# Decide first whether we'll use an MPI script (mpif77) and then test it
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_SET_MPIF77,
[

if test "X${MPI_ENABLED}" = "Xyes" && test "X${F77_ENABLED}" = "Xyes" ; then

  AC_ARG_WITH(mpif77,
  [AC_HELP_STRING([--with-mpif77[[[[=ARG]]]]],[specify MPI Fortran compiler to use @<:@mpif77@:>@],
                  [                                ])],
  [
    if [ test "X${withval}" != "Xno" ]; then
      USE_MPIF77="yes"
      MPI_F77="${withval}"
    else
      USE_MPIF77="no"
    fi
  ],
  [
    USE_MPIF77="yes"
    MPI_F77="mpif77"
  ])

  if test "X${USE_MPIF77}" = "Xyes"; then
    SUNDIALS_CHECK_MPIF77
  else
    MPIF77=${F77}
    SUNDIALS_CHECK_F77_WITH_MPI
  fi

fi

])


#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPIF77,
[

# Test the MPIF77 compiler
if test -x ${MPI_F77}; then
  MPI_F77_EXISTS=yes
else
  TEMP_MPI_F77=`basename "${MPI_F77}"`
  AC_CHECK_PROG(MPI_F77_EXISTS, ${TEMP_MPI_F77}, yes, no)
fi

if test "X${MPI_F77_EXISTS}" = "Xyes"; then
  MPIF77_OK=yes
  MPIF77=${MPI_F77}
else
  MPIF77_OK=no
  AC_MSG_WARN([MPI F77 compiler (${MPI_F77}) not found.])
  echo "
       Cannot find a working MPI F77 compiler.
       Try --with-mpif77 to specify a MPI F77 compiler script.
       Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
       to specify all the specific MPI compile options.

       Some modules will not be configured...
       "
fi

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
         Try --with-mpif77 to specify MPI F77 compiler script.
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



#---------------------------------------------------------------------------------------------
# Function to set MPI support for C++
# Unless the user explicitely disables it, MPI support is enabled if nvec_par is present
# Decide first whether we'll use an MPI script (mpicxx) and then test it
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_SET_MPICXX,
[

if test "X${MPI_ENABLED}" = "Xyes" && test "X${CXX_ENABLED}" = "Xyes"; then

  AC_ARG_WITH(mpicxx,
  [AC_HELP_STRING([--with-mpicxx[[[[=ARG]]]]],[specify MPI C++ compiler to use @<:@mpiCC@:>@],
                  [                                ])],
  [
    if [ test "X${withval}" != "Xno" ]; then
      USE_MPICXX="yes"
      MPI_CXX="${withval}"
      MPICXX_FLAG_USED=yes
    else
      USE_MPICXX="no"
      MPICXX_FLAG_USED=yes
    fi
  ],
  [
    USE_MPICXX="yes"
    MPI_CXX="mpiCC"
  ])

  if test "X${USE_MPICXX}" = "Xyes"; then
    SUNDIALS_CHECK_MPICXX
  else
     MPICXX=${CXX}
    SUNDIALS_CHECK_CXX_WITH_MPI
  fi

fi
])


#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_CHECK_MPICXX,
[

# Test the MPICXX compiler

if test -x ${MPI_CXX}; then
  MPI_CXX_EXISTS=yes
else 
  TEMP_MPI_CXX=`basename "${MPI_CXX}"`
  AC_CHECK_PROG(MPI_CXX_EXISTS, ${TEMP_MPI_CXX}, yes, no)
fi

if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
  MPICXX_OK=yes
  MPICXX=${MPI_CXX}
else
  MPICXX_OK=no
  AC_MSG_WARN([MPI C++ compiler (${MPI_CXX}) not found.])
  echo "
       Cannot find a working MPI C++ compiler.
       Try --with-mpicxx to specify a MPI C++ compiler script.
       Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir
       to specify all the specific MPI compile options.

       Some modules will not be configured...
       "
fi

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
         Try --with-mpicxx to specify MPI C++ compiler script.
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
# Function to decide which modules should be configured
# This function looks through the list of released modules
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_ENABLE_MODULES,
[
AC_ARG_ENABLE([],[],[])

CVODE_ENABLED=no
test -d ${srcdir}/cvode && CVODE_ENABLED=yes
AC_ARG_ENABLE(cvode,
[AC_HELP_STRING([--disable-cvode],[disable configuration of cvode])],
[
  if test "X${enableval}" = "Xno"; then
    CVODE_ENABLED=no
  fi
])

CVODES_ENABLED=no
test -d ${srcdir}/cvodes && CVODES_ENABLED=yes
AC_ARG_ENABLE(cvodes,
[AC_HELP_STRING([--disable-cvodes],[disable configuration of cvodes])],
[
  if test "X${enableval}" = "Xno"; then
    CVODES_ENABLED=no
  fi
])

IDA_ENABLED=no
test -d ${srcdir}/ida && IDA_ENABLED=yes
AC_ARG_ENABLE(ida,
[AC_HELP_STRING([--disable-ida],[disable configuration of ida])],
[
  if test "X${enableval}" = "Xno"; then
    IDA_ENABLED=no
  fi
])

KINSOL_ENABLED=no
test -d ${srcdir}/kinsol && KINSOL_ENABLED=yes
AC_ARG_ENABLE(kinsol,
[AC_HELP_STRING([--disable-kinsol],[disable configuration of kinsol])],
[
  if test "X${enableval}" = "Xno"; then
    KINSOL_ENABLED=no
  fi
])

])

#---------------------------------------------------------------------------------------------
# Function to decide which modules should be configured
# This function looks through the list of development modules
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_ENABLE_DEV_MODULES,
[
IDAS_ENABLED=no
test -d ${srcdir}/idas && IDAS_ENABLED=yes
AC_ARG_ENABLE(idas,
[AC_HELP_STRING([--disable-idas],[disable configuration of idas])],
[
  if test "X${enableval}" = "Xno"; then
    IDAS_ENABLED=no
  fi
])
])

#---------------------------------------------------------------------------------------------
# Function to decide if and what examples will be built
# Unless the user explicitely disables example building, all C examples will
# be built, C++ examples are built only if either cvodes or idas are present,
# and Fortran examples are built if either cvode or kinsol are present.
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_ENABLE_EXAMPLES,
[
AC_ARG_ENABLE([],[ ],[])

EXAMPLES=yes
AC_ARG_ENABLE(examples,
[AC_HELP_STRING([--disable-examples],[disable configuration of examples])],
[
  if test "X${enableval}" = "Xno"; then
    EXAMPLES=no
  fi
])

CXX_EXAMPLES=no
if test "X${EXAMPLES}" = "Xyes" && test "X${CXX_ENABLED}" = "Xyes"; then
  if test "X${CVODES_ENABLED}" = "Xyes" || test "X${IDAS_ENABLED}" = "Xyes"; then
    CXX_EXAMPLES=yes
  fi
fi

F77_EXAMPLES=no;
if test "X${EXAMPLES}" = "Xyes" && test "X${F77_ENABLED}" = "Xyes"; then
  if test "X${CVODE_ENABLED}" = "Xyes" || test "X${KINSOL_ENABLED}" = "Xyes"; then
    F77_EXAMPLES=yes
  fi
fi

if test "X${F77_EXAMPLES}" = "Xno"; then
  F77_ENABLED=no
fi

# Test whether examples can be built
# ----------------------------------

SERIAL_C_EXAMPLES=no
SERIAL_CXX_EXAMPLES=no
SERIAL_F77_EXAMPLES=no
PARALLEL_C_EXAMPLES=no
PARALLEL_CXX_EXAMPLES=no
PARALLEL_F77_EXAMPLES=no

# C examples

if test "X${EXAMPLES}" = "Xyes"; then

    AC_MSG_CHECKING([whether we can build the serial C examples])
    test "X${CC_OK}" = "Xyes" && SERIAL_C_EXAMPLES=yes
    AC_MSG_RESULT([${SERIAL_C_EXAMPLES}])

  if test "X${MPI_ENABLED}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel C examples])
    test "X${MPICC_OK}" = "Xyes" && PARALLEL_C_EXAMPLES=yes
    AC_MSG_RESULT([${PARALLEL_C_EXAMPLES}])
  fi

fi

# C++ examples

if test "X${CXX_EXAMPLES}" = "Xyes"; then

    AC_MSG_CHECKING([whether we can build the serial C++ examples])
    test "X${CXX_OK}" = "Xyes" && SERIAL_CXX_EXAMPLES=yes
    AC_MSG_RESULT([${SERIAL_CXX_EXAMPLES}])

  if test "X${MPI_ENABLED}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel C++ examples])
    test "X${MPICXX_OK}" = "Xyes" && PARALLEL_CXX_EXAMPLES=yes
    AC_MSG_RESULT([${PARALLEL_CXX_EXAMPLES}])
  fi

fi

# Fortran examples

if test "X${F77_EXAMPLES}" = "Xyes"; then

    AC_MSG_CHECKING([whether we can build the serial Fortran examples])
    test "X${F77_OK}" = "Xyes" && SERIAL_F77_EXAMPLES=yes
    AC_MSG_RESULT([${SERIAL_F77_EXAMPLES}])

  if test "X${MPI_ENABLED}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel Fortran examples])
    test "X${MPIF77_OK}" = "Xyes" && PARALLEL_F77_EXAMPLES=yes
    AC_MSG_RESULT([${PARALLEL_F77_EXAMPLES}])
  fi

fi

])

#---------------------------------------------------------------------------------------------
# Function to decide if and what dev examples will be built
# Unless the user explicitely disables example building, all C dev examples will
# be built and Fortran dev examples are built if either cvode or kinsol are present.
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_ENABLE_DEV_EXAMPLES,
[
AC_ARG_ENABLE([],[ ],[])

# Test whether dev examples can be built
# ----------------------------------

SERIAL_DEV_F77_EXAMPLES=no
PARALLEL_DEV_F77_EXAMPLES=no
PARALLEL_DEV_C_EXAMPLES=no

# C examples

if test "X${EXAMPLES}" = "Xyes"; then

  if test "X${MPI_ENABLED}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel dev C examples])
    test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && PARALLEL_DEV_C_EXAMPLES=yes
    AC_MSG_RESULT([${PARALLEL_DEV_C_EXAMPLES}])
  fi

fi

# Fortran examples

if test "X${SUNDIALS_DEV_EXAMPLES}" = "Xyes"; then

  if test "X${F77_EXAMPLES}" = "Xyes"; then

    AC_MSG_CHECKING([whether we can build the serial dev Fortran examples])
    test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && SERIAL_DEV_F77_EXAMPLES=yes
    AC_MSG_RESULT([${SERIAL_DEV_F77_EXAMPLES}])

    if test "X${MPI_ENABLED}" = "Xyes"; then
      AC_MSG_CHECKING([whether we can build the parallel dev Fortran examples])
      test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && PARALLEL_DEV_F77_EXAMPLES=yes
      AC_MSG_RESULT([${PARALLEL_DEV_F77_EXAMPLES}])
    fi

  fi

fi

])

#---------------------------------------------------------------------------------------------
# Function to collect the list of modules to be built
# SUNDIALS_MAKEFILES will contain the list of all makefiles to be generated
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

if test -d ${srcdir}/nvec_par && test "X${MPICC_OK}" = "Xyes"; then
  NVEC_MODULES="$NVEC_MODULES nvec_par"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} nvec_par/Makefile"
fi

# CVODE module

if test "X${CVODE_ENABLED}" = "Xyes"; then

  MODULES="$MODULES cvode/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/source/Makefile"

  MODULES="$MODULES cvode/fcmix"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvode/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvode/examples_par; then
    EX_MODULES="$EX_MODULES cvode/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/examples_par/Makefile"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvode/fcmix/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/examples_ser/Makefile"
  fi

  if test "X${SUNDIALS_DEV_EXAMPLES}" = "Xyes"; then
    if test "X${SERIAL_DEV_F77_EXAMPLES}" = "Xyes"; then
      EX_MODULES="$EX_MODULES cvode/fcmix/test_examples_ser"
      SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/test_examples_ser/Makefile"
    fi
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvode/fcmix/examples_par; then
    EX_MODULES="$EX_MODULES cvode/fcmix/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/examples_par/Makefile"
  fi

  if test "X${SUNDIALS_DEV_EXAMPLES}" = "Xyes"; then
    if test "X${PARALLEL_DEV_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvode/fcmix/test_examples_par; then
      EX_MODULES="$EX_MODULES cvode/fcmix/test_examples_par"
      SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/test_examples_par/Makefile"
    fi
  fi

fi

# CVODES module

if test "X${CVODES_ENABLED}" = "Xyes"; then

  MODULES="$MODULES cvodes/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES cvodes/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvodes/examples_par; then
    EX_MODULES="$EX_MODULES cvodes/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/examples_par/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvodes/test_examples_par; then
    EX_MODULES="$EX_MODULES cvodes/test_examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/test_examples_par/Makefile"
  fi

fi

# IDA module

if test "X${IDA_ENABLED}" = "Xyes"; then

  MODULES="$MODULES ida/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES ida/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/ida/examples_par; then
    EX_MODULES="$EX_MODULES ida/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/examples_par/Makefile"
  fi

fi

# IDAS module

if test "X${IDAS_ENABLED}" = "Xyes"; then

  MODULES="$MODULES idas/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES idas/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/idas/examples_par; then
    EX_MODULES="$EX_MODULES idas/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/examples_par/Makefile"
  fi

fi

# KINSOL module

if test "X${KINSOL_ENABLED}" = "Xyes"; then

  MODULES="$MODULES kinsol/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/source/Makefile"

  MODULES="$MODULES kinsol/fcmix"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES kinsol/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/kinsol/examples_par; then
    EX_MODULES="$EX_MODULES kinsol/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/examples_par/Makefile"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes"; then
    EX_MODULES="$EX_MODULES kinsol/fcmix/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/examples_ser/Makefile"
  fi

  if test "X${SUNDIALS_DEV_EXAMPLES}" = "Xyes"; then
    if test "X${SERIAL_DEV_F77_EXAMPLES}" = "Xyes"; then
      EX_MODULES="$EX_MODULES kinsol/fcmix/test_examples_ser"
      SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/test_examples_ser/Makefile"
    fi
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/kinsol/fcmix/examples_par; then
    EX_MODULES="$EX_MODULES kinsol/fcmix/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/examples_par/Makefile"
  fi

  if test "X${SUNDIALS_DEV_EXAMPLES}" = "Xyes"; then
    if test "X${PARALLEL_DEV_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/kinsol/fcmix/test_examples_par; then
      EX_MODULES="$EX_MODULES kinsol/fcmix/test_examples_par"
      SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/test_examples_par/Makefile"
    fi
  fi

fi

])

#---------------------------------------------------------------------------------------------
# Function to decide the install path
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_INSTALL_PATH,
[
if test "${prefix}" != "NONE"; then
 my_includedir="${prefix}/include"
else
 my_includedir="${DEFAULT_PREFIX}/include"
fi

if test "${exec_prefix}" != "NONE"; then
 my_libdir="${exec_prefix}/lib"
else
 if test "${prefix}" != "NONE"; then
  my_libdir="${prefix}/lib"
 else
  my_libdir="${DEFAULT_PREFIX}/lib"
 fi
fi

AC_MSG_CHECKING([for lib directory])
AC_MSG_RESULT([${my_libdir}])

AC_MSG_CHECKING([for include directory])
AC_MSG_RESULT([${my_includedir}])

# Create the include and lib directories if they don't exist

test ! -d $my_includedir && mkdir -p $my_includedir
test ! -d $my_libdir && mkdir -p $my_libdir
])

#---------------------------------------------------------------------------------------------
# Report what we're going to do
#---------------------------------------------------------------------------------------------

AC_DEFUN(SUNDIALS_REPORT,
[
echo "
Configuration:
--------------

  Host System Type:           $host
  Source code location:       $srcdir
  Install path (include):     $my_includedir
  Install path (lib):         $my_libdir

  C/C++ Preprocessor:         $CPP $CPPFLAGS
  C Compiler:	              $CC $CFLAGS
  C Linker:                   $CC $LDFLAGS $LIBS"

if test "X${CXX_OK}" = "Xyes"; then
echo "  
  C++ Compiler:               $CXX $CXXFLAGS"
fi

if test "X${F77_OK}" = "Xyes"; then
echo "
  F77 Compiler:               $F77 $FFLAGS"
fi

if test "X${MPICC_OK}" = "Xyes"; then
echo "  
  MPICC:                      $MPICC"
fi

if test "X${MPICXX_OK}" = "Xyes"; then
echo "  
  MPICXX:                     $MPICXX"
fi

if test "X${MPIF77_OK}" = "Xyes"; then
echo "  
  MPIF77:                     $MPIF77"
fi

if test "X${MPICC_OK}" = "Xyes"; then
echo "
  MPI include directory:      $MPIINC
  MPI lib directory:          $MPILIBDIR
  MPI libraries               $MPILIBS"
fi

echo "  
  Type 'make' and then 'make install' to build and install ${PACKAGE_STRING}"

echo "
Modules:
--------
"
if test "X${CVODE_ENABLED}" = "Xyes"; then
echo "  cvode"
fi
if test "X${CVODES_ENABLED}" = "Xyes"; then
echo "  cvodes"
fi
if test "X${IDA_ENABLED}" = "Xyes"; then
echo "  ida"
fi
if test "X${IDAS_ENABLED}" = "Xyes"; then
echo "  idas"
fi
if test "X${KINSOL_ENABLED}" = "Xyes"; then
echo "  kinsol"
fi

if test "X${EXAMPLES}" = "Xyes"; then
echo "
Examples:
---------
  
  Type 'make examples' to build the following examples:
"
if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
echo "  Serial C examples"
fi
if test "X${PARALLEL_C_EXAMPLES}" = "Xyes"; then
if test "X${MPICC_OK}" = "Xyes"; then
echo "  Parallel C examples"
fi
if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes"; then
if test "X${MPICC_OK}" = "Xyes"; then
echo "  Parallel dev C examples"
fi
fi
fi
if test "X${SERIAL_CXX_EXAMPLES}" = "Xyes"; then
if test "X${CXX_OK}" = "Xyes"; then
echo "  Serial C++ examples"
fi
fi
if test "X${PARALLEL_CXX_EXAMPLES}" = "Xyes"; then
if test "X${MPICXX_OK}" = "Xyes"; then
echo "  Parallel C++ examples"
fi
fi
if test "X${SERIAL_F77_EXAMPLES}" = "Xyes"; then
if test "X${F77_OK}" = "Xyes"; then
echo "  Serial Fortran examples"
fi
fi
if test "X${SUNDIALS_DEV_EXAMPLES}" = "Xyes"; then
if test "X${SERIAL_DEV_F77_EXAMPLES}" = "Xyes"; then
if test "X${F77_OK}" = "Xyes"; then
echo "  Serial dev Fortran examples"
fi
fi
fi
if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
if test "X${MPIF77_OK}" = "Xyes"; then
echo "  Parallel Fortran examples"
fi
fi
if test "X${SUNDIALS_DEV_EXAMPLES}" = "Xyes"; then
if test "X${PARALLEL_DEV_F77_EXAMPLES}" = "Xyes"; then
if test "X${MPIF77_OK}" = "Xyes"; then
echo "  Parallel dev Fortran examples"
fi
fi
fi

fi


echo "
----------------------------------
Finished SUNDIALS Configure Script
----------------------------------
"

])
