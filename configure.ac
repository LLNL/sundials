# ------------------------------------------------------------------------
# $Revision: 1.11 $
# $Date: 2004-10-21 20:18:23 $
# ------------------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# ------------------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# Initialize configure
AC_INIT(SUNDIALS, 2.0, radu@llnl.gov)

# Display greeting
SUNDIALS_GREETING

# Miscellaneous SUNDIALS initializations
SUNDIALS_INITIALIZE

# Set C compilation
SUNDIALS_SET_CC

# Set Fortran compilation
SUNDIALS_SET_F77

# Set C++ compilation
SUNDIALS_SET_CXX

# Set MPI support
SUNDIALS_SET_MPICC
SUNDIALS_SET_MPIF77
SUNDIALS_SET_MPICXX

# Set modules to build
SUNDIALS_ENABLE_MODULES
SUNDIALS_ENABLE_DEV_MODULES

# Set examples modules
SUNDIALS_ENABLE_EXAMPLES
SUNDIALS_ENABLE_DEV_EXAMPLES

# Run libtool checks
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# Substitute variables
AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(MPIF77)
AC_SUBST(MPI_INC_DIR)
AC_SUBST(MPI_LIB_DIR)
AC_SUBST(MPI_LIBS)
AC_SUBST(MODULES)
AC_SUBST(NVEC_MODULES)
AC_SUBST(EXAMPLE_MODULES)
AC_SUBST(F77_FLOAT_TYPE)
AC_SUBST(LIBTOOL_DEPS)

# Set installation path
SUNDIALS_INSTALL_PATH

# Generate makefiles
SUNDIALS_BUILD_MODULES_LIST
AC_CONFIG_FILES([${SUNDIALS_MAKEFILES}])
AC_OUTPUT()

# Display summary
SUNDIALS_REPORT
