# ------------------------------------------------------------------------
# $Revision: 1.9 $
# $Date: 2004-03-30 19:36:03 $
# ------------------------------------------------------------------------
# Programmers: Radu Serban @ LLNL
# ------------------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California
# Produced at the Lawrence Livermore National Laboratory
# All rights reserved
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# Initialize configure
AC_INIT(Sundials, 1.1, radu@llnl.gov)

# Misc. sundials initializations
SUNDIALS_INITIALIZE

# Check for user overrides (some to be used later)
SUNDIALS_CHECK_OVERRIDES

# Archiver
SUNDIALS_CHECK_AR

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

# Set modules to configure
SUNDIALS_ENABLE_MODULES
SUNDIALS_ENABLE_DEV_MODULES

# Set example modules
SUNDIALS_ENABLE_EXAMPLES

# Substitute variables
AC_SUBST(FLIBS)
AC_SUBST(FCFLAGS)
AC_SUBST(MPICC)
AC_SUBST(MPICXX)
AC_SUBST(MPIF77)
AC_SUBST(MPIINC)
AC_SUBST(MPILIBDIR)
AC_SUBST(MPILIBS)
AC_SUBST(MODULES)
AC_SUBST(NVEC_MODULES)
AC_SUBST(EX_MODULES)

# Set install path
SUNDIALS_INSTALL_PATH

# Set the lists of modules and generate makefiles
SUNDIALS_BUILD_MODULES_LIST
AC_CONFIG_FILES([${SUNDIALS_MAKEFILES}])
AC_OUTPUT()

# Report what we've decide to do
SUNDIALS_REPORT

#that's all
