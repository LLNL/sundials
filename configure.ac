# ------------------------------------------------------------------------
# File: configure.ac
# Programmers: Radu Serban @ LLNL
# ------------------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California
# Produced at the Lawrence Livermore National Laboratory
# All rights reserved
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# This must be the first line in configure.ac.
# Optional 3rd argument is email address for bugs.

AC_INIT(Sundials, 1.1, radu@llnl.gov)

# Say Hi!
USERNAME=`whoami`
echo "----------------------------------"
echo "Hi ${USERNAME}!"
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
   AC_CANONICAL_TARGET
fi

# Overwrite the default installation path
DEFAULT_PREFIX=`pwd`
AC_PREFIX_DEFAULT(`pwd`)

# If make predefines the Make variable MAKE, define output variable 
# SET_MAKE to be empty. Otherwise, define SET_MAKE to contain `MAKE=make'. 

AC_PROG_MAKE_SET

# ------------------------------------------------------------------------
# Check for user overrides (some to be used later)
# ------------------------------------------------------------------------
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

# ------------------------------------------------------------------------
# Checks for extra flags, etc
# ------------------------------------------------------------------------

AC_ARG_WITH(cppflags,
[AC_HELP_STRING([--with-cppflags],[add extra preprocessor flags])],
[EXTRA_CPPFLAGS=${withval}]
)

AC_ARG_WITH(cflags,
[AC_HELP_STRING([--with-cflags],[add extra C compiler flags])],
[EXTRA_CFLAGS=${withval}]
)

AC_ARG_WITH(cxxflags,
[AC_HELP_STRING([--with-cxxflags],[add extra C++ compiler flags])],
[EXTRA_CXXFLAGS=${withval}]
)

AC_ARG_WITH(fflags,
[AC_HELP_STRING([--with-fflags],[add extra Fortran compiler flags])],
[EXTRA_FFLAGS=${withval}]
)

AC_ARG_WITH(ldflags,
[AC_HELP_STRING([--with-ldflags],[add extra linker flags])],
[EXTRA_LDFLAGS=${withval}]
)

# ------------------------------------------------------------------------
# Archiver
# ------------------------------------------------------------------------

AC_ARG_WITH(ar,
[AC_HELP_STRING([--with-ar],[specify archiver command @<:@"ar rc"@:>@])],
[
  AC_MSG_CHECKING(user-defined archiver program)
  AC_MSG_RESULT([${withval}])
  AR=${withval}
],
[
  AC_CHECK_PROG(AR,ar,ar,no)
  if test "X${AR}" = "Xno"; then
    AC_MSG_ERROR(cannot find an archiver)
  else
    AR="${AR} rc"
  fi
])

# ------------------------------------------------------------------------
# See whether we need a C++ compiler and MPI support 
# ------------------------------------------------------------------------

if test -d ${srcdir}/nvec_par; then
  PARALLEL=yes
else
  PARALLEL=no
fi


if test -d ${srcdir}/nvec_ser; then
  SERIAL=yes
else
  SERIAL=no
fi

if test -d ${srcdir}/xs4c; then
  COMPLEX=yes
else
  COMPLEX=no
fi

# ------------------------------------------------------------------------
# See whether we will compile examples
# ------------------------------------------------------------------------

AC_ARG_ENABLE(examples,
   AC_HELP_STRING([--enable-examples],[compile examples @<:@no@:>@]),
   EXAMPLES=yes,
   EXAMPLES=no)

CXX_EXAMPLES=no
if test "X${EXAMPLES}" = "Xyes" && test "X${COMPLEX}" = "Xyes"; then
  if test -d ${srcdir}/cvodes || test -d ${srcdir}/idas; then
    CXX_EXAMPLES=yes
  fi
fi

F77_EXAMPLES=no;
if test "X${EXAMPLES}" = "Xyes"; then
  if test -d ${srcdir}/cvode || test -d ${srcdir}/kinsol; then
    F77_EXAMPLES=yes
  fi
fi

# ------------------------------------------------------------------------
# some initializations
# ------------------------------------------------------------------------
CC_OK=no
CXX_OK=no
F77_OK=no
MPICC_OK=no
MPICXX_OK=no
MPIF77_OK=no
# ------------------------------------------------------------------------
# Check for preprocessor flags
# ------------------------------------------------------------------------

AC_MSG_CHECKING([for extra preprocessor flags])
if test -z "${EXTRA_CPPFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CPPFLAGS}])
  CPPFLAGS="${CPPFLAGS} ${EXTRA_CPPFLAGS}"
fi

# ------------------------------------------------------------------------
# C compiler and flags
# ------------------------------------------------------------------------

if test "X${MY_CC}" = "Xnone"; then
  SUNDIALS_PROG_CC
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

# ------------------------------------------------------------------------
# C++ compiler and flags
# ------------------------------------------------------------------------

if test "X${COMPLEX}" = "Xyes"; then

  if test "X${MY_CXX}" = "Xnone"; then
    SUNDIALS_PROG_CXX
  fi
  AC_PROG_CXX(${MY_CXX} CC g++ gcc c++ cxx)

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

# ------------------------------------------------------------------------
# Fortran compiler and flags
# ------------------------------------------------------------------------

if test "X${F77_EXAMPLES}" = "Xyes"; then

  if test "X${MY_F77}" = "Xnone"; then
    SUNDIALS_PROG_F77
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

# ------------------------------------------------------------------------
# Linker flags
# ------------------------------------------------------------------------

AC_MSG_CHECKING([for extra linker flags])
if test -z "${EXTRA_LDFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_LDFLAGS}])
  LDFLAGS="${LDFLAGS} ${EXTRA_LDFLAGS}"
fi

# ------------------------------------------------------------------------
# Check for MPI
# ------------------------------------------------------------------------

if test "X${PARALLEL}" = "Xyes"; then

  # Decide whether we use MPI compilers (mpicc, mpicxx, mpif77)
  # or if we use the compiler with special flags

  SUNDIALS_MPI_SPECIFY

  # Test if we can use MPI

  if test "X${USE_MPI_COMPILERS}" = "Xyes"; then

    SUNDIALS_CHECK_MPICC

    if test "X${COMPLEX}" = "Xyes"; then
      SUNDIALS_CHECK_MPICXX
    fi

    if test "X${F77_EXAMPLES}" = "Xyes"; then
      SUNDIALS_CHECK_MPIF77
    fi

  else

    MPICC=${CC}
    SUNDIALS_CHECK_CC_WITH_MPI

    if test "X${COMPLEX}" = "Xyes"; then
      MPICXX=${CXX}
      SUNDIALS_CHECK_CXX_WITH_MPI
    fi

    if test "X${F77_EXAMPLES}" = "Xyes"; then
      MPIF77=${F77}
      SUNDIALS_CHECK_F77_WITH_MPI
    fi

  fi

fi

# ------------------------------------------------------------------------
# Check whether we have everything we need to build the examples
# ------------------------------------------------------------------------

SERIAL_C_EXAMPLES=no
SERIAL_CXX_EXAMPLES=no
SERIAL_F77_EXAMPLES=no
PARALLEL_C_EXAMPLES=no
PARALLEL_CXX_EXAMPLES=no
PARALLEL_F77_EXAMPLES=no

# C examples

if test "X${EXAMPLES}" = "Xyes"; then

  if test "X${SERIAL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the serial C examples])
    test "X${CC_OK}" = "Xyes" && SERIAL_C_EXAMPLES=yes
    AC_MSG_RESULT([${SERIAL_C_EXAMPLES}])
  fi

  if test "X${PARALLEL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel C examples])
    test "X${MPICC_OK}" = "Xyes" && PARALLEL_C_EXAMPLES=yes
    AC_MSG_RESULT([${PARALLEL_C_EXAMPLES}])
  fi

fi

# C++ examples

if test "X${CXX_EXAMPLES}" = "Xyes"; then

  if test "X${SERIAL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the serial C++ examples])
    test "X${CXX_OK}" = "Xyes" && SERIAL_CXX_EXAMPLES=yes
    AC_MSG_RESULT([${SERIAL_CXX_EXAMPLES}])
  fi

  if test "X${PARALLEL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel C++ examples])
    test "X${MPICXX_OK}" = "Xyes" && PARALLEL_CXX_EXAMPLES=yes
    AC_MSG_RESULT([${PARALLEL_CXX_EXAMPLES}])
  fi

fi

# Fortran examples

if test "X${F77_EXAMPLES}" = "Xyes"; then

  if test "X${SERIAL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the serial Fortran examples])
    test "X${F77_OK}" = "Xyes" && SERIAL_F77_EXAMPLES=yes
    AC_MSG_RESULT([${SERIAL_F77_EXAMPLES}])
  fi

  if test "X${PARALLEL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel Fortran examples])
    test "X${MPIF77_OK}" = "Xyes" && PARALLEL_F77_EXAMPLES=yes
    AC_MSG_RESULT([${PARALLEL_F77_EXAMPLES}])
  fi

fi

# ------------------------------------------------------------------------ 
# Check where we are going to install
# ------------------------------------------------------------------------ 

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

# ------------------------------------------------------------------------
# Set the lists containing the solver, nvector, and examples modules
# ------------------------------------------------------------------------

SUNDIALS_BUILD_MODULES_LIST

# ------------------------------------------------------------------------
# Substitute variables
# ------------------------------------------------------------------------

AC_SUBST(FLIBS)

AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(MPIF77)

AC_SUBST(MPIINC)
AC_SUBST(MPILIBDIR)
AC_SUBST(MPILIBS)

AC_SUBST(MODULES)
AC_SUBST(NVEC_MODULES)
AC_SUBST(EX_MODULES)

# ------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------

AC_CONFIG_FILES([${SUNDIALS_MAKEFILES}])

AC_OUTPUT()

# Say Bye!

echo "
Configuration:
--------------

  Host System Type:           $host
  Source code location:       $srcdir
  Install path (include):     $my_includedir
  Install path (lib):         $my_libdir

  C Preprocessor:             $CPP $CPPFLAGS
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
  MPICC:                      $MPICC
  MPI include directory:      $MPIINC
  MPI lib directory:          $MPILIBDIR
  MPI libraries               $MPILIBS"
fi

if test "X${MPICXX_OK}" = "Xyes"; then
echo "  
  MPICXX:                     $MPICXX"
fi

if test "X${MPIF77_OK}" = "Xyes"; then
echo "  
  MPIF77:                     $MPIF77"
fi

echo "  
  Type 'make' and then 'make install' to build and install ${PACKAGE_STRING}"

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
if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
if test "X${MPIF77_OK}" = "Xyes"; then
echo "  Parallel Fortran examples"
fi
fi

fi


echo "
----------------------------------
Finished SUNDIALS Configure Script
----------------------------------
"

#that's all