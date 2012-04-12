# ------------------------------------------------------------------------
# $Revision: 1.53 $
# $Date: 2012-04-12 00:40:13 $
# ------------------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# ------------------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# Initialize configure
AC_INIT(SUNDIALS, 2.5.0, radu@llnl.gov)

# Display greeting
SUNDIALS_GREETING

# Specify directory containing auxillary build tools and M4 files
AC_CONFIG_AUX_DIR([config])

# Miscellaneous SUNDIALS initializations
echo "Initialization"
echo "--------------"
echo ""
SUNDIALS_INITIALIZE

# Test enable/disable features
SUNDIALS_ENABLES

# Set C compilation (Required)
echo ""
echo "C Compiler Settings"
echo "-------------------"
echo ""
AC_PROG_CC(cc gcc)
SUNDIALS_SET_CC

# Set Fortran support
if test "X${FCMIX_ENABLED}" = "Xyes" || test "X${LAPACK_ENABLED}" = "Xyes"; then
  echo ""
  echo "Fortran Settings"
  echo "----------------"
  echo ""
  SUNDIALS_F77_SUPPORT
fi

# Set MPI support (Optional)
if test "X${MPI_ENABLED}" = "Xyes"; then
  echo ""
  echo "MPI-C Settings"
  echo "--------------"
  echo ""
  SUNDIALS_SET_MPICC
fi
if test "X${MPI_C_COMP_OK}" = "Xyes"; then
  if test "X${FCMIX_ENABLED}" = "Xyes"; then
    echo ""
    echo "MPI-2 Functionality"
    echo "-------------------"
    echo ""
    SUNDIALS_CHECK_MPI2
  fi
  if test "X${F77_EXAMPLES_ENABLED}" = "Xyes"; then
    echo ""
    echo "MPI-Fortran Settings"
    echo "--------------------"
    echo ""
    SUNDIALS_SET_MPIF77
  fi
fi

# Set examples modules
if test "X${EXAMPLES_ENABLED}" = "Xyes"; then
  echo ""
  echo "Examples"
  echo "--------"
  echo ""
  SUNDIALS_SET_EXAMPLES
fi

# Run libtool checks
echo ""
echo "Libtool Settings"
echo "----------------"
echo ""
AC_DISABLE_SHARED
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

# Do we build shared libraries?
# (enable_shared is set by the function AC_ENABLE_SHARED in libtool.m4)
SHARED_LIBS="${enable_shared}"

# Print a few more lines for configure --help
SUNDIALS_MORE_HELP

# Substitute variables
AC_SUBST(SHARED_LIBS)
AC_SUBST(MPICC)
AC_SUBST(MPIF77)
AC_SUBST(MPI_INC_DIR)
AC_SUBST(MPI_LIB_DIR)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_FLAGS)
AC_SUBST(FCMIX_ENABLED)
AC_SUBST(FLOAT_TYPE)
AC_SUBST(LIBTOOL_DEPS)
AC_SUBST(F77_MANGLE_MACRO1)
AC_SUBST(F77_MANGLE_MACRO2)
AC_SUBST(F77_CASE)
AC_SUBST(F77_UNDERSCORES)
AC_SUBST(PRECISION_LEVEL)
AC_SUBST(GENERIC_MATH_LIB)
AC_SUBST(BLAS_LAPACK_MACRO)
AC_SUBST(SUNDIALS_EXPORT)
AC_SUBST(F77_MPI_COMM_F2C)
AC_SUBST(F77_LNKR)
AC_SUBST(F77_LIBS)
AC_SUBST(F77_LDFLAGS)
AC_SUBST(LAPACK_ENABLED)
AC_SUBST(BLAS_LAPACK_LIBS)
AC_SUBST(MPIF77_LNKR)
AC_SUBST(OBJEXT)
AC_SUBST(EXEEXT)
AC_SUBST(SLV_MODULES)
AC_SUBST(EXS_MODULES)
AC_SUBST(EXS_INSTDIR)

# Generate makefiles
echo ""
echo "Create Makefiles and configuration files"
echo "----------------------------------------"
echo ""
SUNDIALS_BUILD_MODULES_LIST

# Set-up postprocessing commands
if test "X${EXAMPLES_ENABLED}" = "Xyes" && test "X${EXS_INSTDIR}" != "Xno"; then
  AC_CONFIG_COMMANDS([postprocessing],
  [
    echo ""
    echo "Postprocess exported Makefiles"
    echo "------------------------------"
    echo ""
  ])
  SUNDIALS_POST_PROCESSING
fi

# Specify files to be configured and call AC_OUTPUT
AC_CONFIG_FILES([${SUNDIALS_MAKEFILES}])
AC_CONFIG_FILES([${SUNDIALS_CONFIGFILES}])
AC_OUTPUT

# Display summary
SUNDIALS_REPORT
