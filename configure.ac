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

# Define $build, $host, $target, etc
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

# Overwrite the default installation path
DEFAULT_PREFIX=`pwd`
AC_PREFIX_DEFAULT(`pwd`)

# If make predefines the Make variable MAKE, define output variable 
# SET_MAKE to be empty. Otherwise, define SET_MAKE to contain `MAKE=make'. 

AC_PROG_MAKE_SET

# ------------------------------------------------------------------------
# Check for user overrides (some to be used later)
# ------------------------------------------------------------------------

if test -n "${CC}"; then
  USERDEF_CC=yes
fi

if test -n "${CXX}"; then
  USERDEF_CXX=yes
fi

if test -n "${F77}"; then
  USERDEF_F77=yes
fi

if test -n "${CPPFLAGS}"; then
  USERDEF_CPPFLAGS=yes
fi

if test -n "${CFLAGS}"; then
  USERDEF_CFLAGS=yes
fi

if test -n "${CXXFLAGS}"; then
  USERDEF_CXXFLAGS=yes
fi

if test -n "${FFLAGS}"; then
  USERDEF_FFLAGS=yes
fi

if test -n "${LDFLAGS}"; then
  USERDEF_LDFLAGS=yes
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

if test "X${USERDEF_CC}" != "Xyes"; then
  SUNDIALS_PROG_CC
fi
AC_CHECK_PROG(CC_OK, ${CC}, yes, no)

AC_MSG_CHECKING([for C compiler flags])
if test "X${USERDEF_CFLAGS}" = "Xyes"; then
  AC_MSG_RESULT([${CFLAGS} (user-specified via CFLAGS)])
else
  CFLAGS="-O"
  SUNDIALS_PROG_CFLAGS
  AC_MSG_RESULT([${CFLAGS}])
fi

AC_MSG_CHECKING([for extra C compiler flags])
if test -z "${EXTRA_CFLAGS}"; then
  AC_MSG_RESULT(none)
else
  AC_MSG_RESULT([${EXTRA_CFLAGS}])
  CFLAGS="${CFLAGS} ${EXTRA_CFLAGS}"
fi

# ------------------------------------------------------------------------
# C++ compiler and flags
# ------------------------------------------------------------------------

if test "X${COMPLEX}" = "Xyes"; then

  if test "X${USERDEF_CXX}" != "Xyes"; then
    SUNDIALS_PROG_CXX
  fi
  AC_CHECK_PROG(CXX_OK, ${CXX}, yes, no)

  AC_MSG_CHECKING([for C++ compiler flags])
  if test "X${USERDEF_CXXFLAGS}" = "Xyes"; then
    AC_MSG_RESULT([${CXXFLAGS} (user-specified via CXXFLAGS)])
  else
    CXXFLAGS="-O"
    SUNDIALS_PROG_CXXFLAGS
    AC_MSG_RESULT([${CXXFLAGS}])
  fi

  AC_MSG_CHECKING([for extra C++ compiler flags])
  if test -z "${EXTRA_CXXFLAGS}"; then
    AC_MSG_RESULT(none)
  else
    AC_MSG_RESULT([${EXTRA_CXXFLAGS}])
    CXXFLAGS="${CXXFLAGS} ${EXTRA_CXXFLAGS}" 
  fi

fi

# ------------------------------------------------------------------------
# Fortran compiler and flags
# ------------------------------------------------------------------------

if test "X${F77_EXAMPLES}" = "Xyes"; then

    if test "X${USERDEF_F77}" != "Xyes"; then
      SUNDIALS_PROG_F77
    fi
    AC_CHECK_PROG(F77_OK, ${F77}, yes, no)

    AC_MSG_CHECKING([for Fortran compiler flags])
    if test "X${USERDEF_FFLAGS}" = "Xyes"; then
      AC_MSG_RESULT([${FFLAGS} (user-specified via FFLAGS)])
    else
      FFLAGS="-O"
      SUNDIALS_PROG_FFLAGS
      AC_MSG_RESULT([${FFLAGS}])
    fi

    AC_MSG_CHECKING([for extra Fortran compiler flags])
    if test -z "${EXTRA_FFLAGS}"; then
      AC_MSG_RESULT(none)
    else
      AC_MSG_RESULT([${EXTRA_FFLAGS}])
      FFLAGS="${FFLAGS} ${EXTRA_FFLAGS}"
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
# Check for required header files
# ------------------------------------------------------------------------

AC_LANG([C])
AC_CHECK_HEADERS(
  [stdio.h],
  break,
  AC_MSG_NOTICE([You do not appear to have stdio.h!])
)

if test "X${EXAMPLES}" = "Xyes"; then

  AC_LANG([C])
  AC_CHECK_HEADERS(
    [stdlib.h],
    break,
    AC_MSG_NOTICE([You do not appear to have stdlib.h!])
  )

  AC_LANG([C])
  AC_CHECK_HEADERS(
    [math.h],
    break,
    AC_MSG_NOTICE([You do not appear to have math.h!])
  )

fi

if test "X${COMPLEX}" = "Xyes"; then

  AC_LANG([C++])
  AC_CHECK_HEADERS(
    [complex.h],
    break,
    AC_MSG_NOTICE([You do not appear to have complex.h!])
  )

fi

# ------------------------------------------------------------------------
# Check for required libraries
# ------------------------------------------------------------------------

if test "X${EXAMPLES}" = "Xyes"; then

  AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
  AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

  if test "X${F77_EXAMPLES}" = "Xyes"; then
    AC_F77_LIBRARY_LDFLAGS
  fi  

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

    SUNDIALS_CHECK_MPICXX

    if test "X${F77_EXAMPLES}" = "Xyes"; then
      SUNDIALS_CHECK_MPIF77
    fi

  else

    MPICC=${CC}
    SUNDIALS_CHECK_CC_WITH_MPI

    MPICXX=${CXX}
    SUNDIALS_CHECK_CXX_WITH_MPI

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
    if test "X${CC_OK}" = "Xyes"; then
      SERIAL_C_EXAMPLES=yes
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no (C compiler problem)])
    fi
  fi

  if test "X${PARALLEL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel C examples])
    if test "X${MPICC_OK}" = "Xyes"; then
      PARALLEL_C_EXAMPLES=yes
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no (C compiler with MPI problem)])
    fi
  fi

fi

# C++ examples

if test "X${CXX_EXAMPLES}" = "Xyes"; then

  if test "X${SERIAL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the serial C++ examples])
    if test "X${CXX_OK}" = "Xyes"; then
      SERIAL_CXX_EXAMPLES=yes
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no (C++ compiler problem)])
    fi
  fi

  if test "X${PARALLEL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel C++ examples])
    if test "X${MPICXX_OK}" = "Xyes"; then
      PARALLEL_CXX_EXAMPLES=yes
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no (C++ compiler with MPI problem)])
    fi
  fi

fi

# Fortran examples

if test "X${F77_EXAMPLES}" = "Xyes"; then

  if test "X${SERIAL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the serial Fortran examples])
    if test "X${F77_OK}" = "Xyes"; then
      SERIAL_F77_EXAMPLES=yes
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no (F77 compiler problem)])
    fi
  fi

  if test "X${PARALLEL}" = "Xyes"; then
    AC_MSG_CHECKING([whether we can build the parallel Fortran examples])
    if test "X${MPIF77_OK}" = "Xyes"; then
      PARALLEL_F77_EXAMPLES=yes
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no (Fortran compiler with MPI problem)])
    fi
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
echo "----------------------------------"
echo "Finished SUNDIALS Configure Script"
echo "----------------------------------"

#that's all