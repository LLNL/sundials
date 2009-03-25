# -----------------------------------------------------------------
# $Revision: 1.57 $
# $Date: 2009-03-25 18:32:37 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
#
# SUNDIALS autoconf macros
#
# The functions defined here fall into the following categories:
#
# (1) Initializations:
#     SUNDIALS_GREETING
#     SUNDIALS_INITIALIZE
#     SUNDIALS_ENABLES
#
# (2) C compiler tests
#     SUNDIALS_SET_CC
#     SUNDIALS_CC_CHECK
#     SUNDIALS_CPLUSPLUS_CHECK
#
# (3) Fortran support
#     SUNDIALS_F77_SUPPORT
#     SUNDIALS_F77_CHECK
#     SUNDIALS_F77_LNKR_CHECK
#     SUNDIALS_F77_NAME_MANGLING
#     SUNDIALS_F77_LAPACK_SET
#
# (4) Parallel support
#     SUNDIALS_SET_MPICC
#     SUNDIALS_CHECK_MPICC
#     SUNDIALS_CC_WITH_MPI_CHECK
#     SUNDIALS_SET_MPIF77
#     SUNDIALS_CHECK_MPIF77
#     SUNDIALS_MPIF77_LNKR_CHECK
#     SUNDIALS_F77_WITH_MPI_CHECK
#     SUNDIALS_CHECK_MPI2
#
# (5) Finalizations:
#     SUNDIALS_MORE_HELP
#     SUNDIALS_SET_EXAMPLES
#     SUNDIALS_BUILD_MODULES_LIST
#     SUNDIALS_POST_PROCESSING
#     SUNDIALS_REPORT
#
# -----------------------------------------------------------------


#=================================================================#
#                                                                 #
#                                                                 #
#               I N I T I A L I Z A T I O N S                     #
#                                                                 #
#                                                                 #
#==================================================================

#------------------------------------------------------------------
# GREETING
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_GREETING],
[

# Say Hi!
echo "
---------------------------------
Running SUNDIALS Configure Script
---------------------------------
"

]) dnl END SUNDIALS_GREETING

#------------------------------------------------------------------
# PERFORM INITIALIZATIONS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_INITIALIZE],
[

# Reference custom macros
m4_include([config/mod_fortran.m4])
m4_include([config/mod_c.m4])
m4_include([config/cust_general.m4])

# Make input filename DOS compatible (change config.h.in to config.hin)
AC_CONFIG_HEADERS([config.h:config.hin])

# Make user aware of copyright notice (input COPYRIGHT information)
AC_COPYRIGHT(
[
Copyright (c) 2002, The Regents of the University of California.
Produced at the Lawrence Livermore National Laboratory.
All rights reserved.
For details, see the LICENSE file.
])

# Specify root of source tree
# Given file is guaranteed to exist in all SUNDIALS packages
AC_CONFIG_SRCDIR([/src/sundials/sundials_nvector.c])

# Get host information
# AC_CANONICAL_BUILD defines the following variables: build, build_cpu,
# build_vendor, and build_os
AC_CANONICAL_BUILD
# AC_CANONICAL_HOST defines the following variables: host, host_cpu,
# host_vendor, and host_os
AC_CANONICAL_HOST

# Set MAKE if necessary
# Must include @SET_MAKE@ in each Makefile.in file
# AC_SUBST is called automatically for SET_MAKE
AC_PROG_MAKE_SET

# Defines INSTALL (sets to path of "install" program)
# Also sets INSTALL_PROGRAM and INSTALL_SCRIPT
AC_PROG_INSTALL

# Set defaults for config/sundials_config.in file
F77_MANGLE_MACRO1=""
F77_MANGLE_MACRO2=""
PRECISION_LEVEL=""
GENERIC_MATH_LIB=""
BLAS_LAPACK_MACRO=""
F77_MPI_COMM_F2C=""
SUNDIALS_EXPORT="#define SUNDIALS_EXPORT"

# Initialize enable status of various modules, options, and features
# to their default values
#
# NOTE: when CPODES is released, change its default to enabled.
#
CVODE_ENABLED="yes"
CVODES_ENABLED="yes"
IDA_ENABLED="yes"
IDAS_ENABLED="yes"
KINSOL_ENABLED="yes"
LAPACK_ENABLED="yes"
FCMIX_ENABLED="yes"
MPI_ENABLED="yes"
#
CPODES_ENABLED="no"
#
EXAMPLES_ENABLED="no"
F77_EXAMPLES_ENABLED="no"

# Initialize variables that may NOT necessarily be initialized
# during normal execution. Should NOT use uninitialized variables
F77_OK="no"
LAPACK_OK="no"
MPI_C_COMP_OK="no"
MPI_F77_COMP_OK="no"

# This variable is set to "yes" if an AC_MSG_WARN statement
# was executed
SUNDIALS_WARN_FLAG="no"

]) dnl END SUNDIALS_INITIALIZE

#------------------------------------------------------------------
# TEST ENABLES
#
# The following variables may be changed here (default value in []):
#
#   CVODE_ENABLED    - enable CVODE module [yes]
#   CVODES_ENABLED   - enable CVODES module [yes]
#   IDA_ENABLED      - enable IDA module [yes]
#   IDAS_ENABLED     - enable IDAS module [yes]
#   KINSOL_ENABLED   - enable KINSOL module [yes]
#   FCMIX_ENABLED    - enable Fortran-C interfaces [yes]
#   LAPACK_ENABLED   - enable Lapack support [yes]
#   MPI_ENABLED      - enable parallel support [yes]
#   EXAMPLES_ENABLED - enable example programs [no]
#   F77_EXAMPLES_ENABLED - enable Fortran example programs [no]
#  
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_ENABLES],
[

# Check if user wants to disable CVODE module
# If not, then make certain source directory actually exists
AC_ARG_ENABLE(cvode,
[AC_HELP_STRING([--disable-cvode],[disable configuration of CVODE])],
[
if test "X${enableval}" = "Xno"; then
  CVODE_ENABLED="no"
fi
],
[
if test -d ${srcdir}/src/cvode ; then
  CVODE_ENABLED="yes"
else
  CVODE_ENABLED="no"
fi
])

# Check if user wants to disable CVODES module
# If not, then make certain source directory actually exists
AC_ARG_ENABLE(cvodes,
[AC_HELP_STRING([--disable-cvodes],[disable configuration of CVODES])],
[
if test "X${enableval}" = "Xno"; then
  CVODES_ENABLED="no"
fi
],
[
if test -d ${srcdir}/src/cvodes ; then
  CVODES_ENABLED="yes"
else
  CVODES_ENABLED="no"
fi
])

# Check if user wants to disable IDA module
# If not, then make certain source directory actually exists
AC_ARG_ENABLE(ida,
[AC_HELP_STRING([--disable-ida],[disable configuration of IDA])],
[
if test "X${enableval}" = "Xno"; then
  IDA_ENABLED="no"
fi
],
[
if test -d ${srcdir}/src/ida  ; then
  IDA_ENABLED="yes"
else
  IDA_ENABLED="no"
fi
])

# Check if user wants to disable IDAS module
# If not, then make certain source directory actually exists
AC_ARG_ENABLE(idas,
[AC_HELP_STRING([--disable-idas],[disable configuration of IDAS])],
[
if test "X${enableval}" = "Xno"; then
  IDAS_ENABLED="no"
fi
],
[
if test -d ${srcdir}/src/idas  ; then
  IDAS_ENABLED="yes"
else
  IDAS_ENABLED="no"
fi
])

# Check if user wants to disable KINSOL MODULE
# If not, then make certain source directory actually exists
AC_ARG_ENABLE(kinsol,
[AC_HELP_STRING([--disable-kinsol],[disable configuration of KINSOL])],
[
if test "X${enableval}" = "Xno"; then
  KINSOL_ENABLED="no"
fi
],
[
if test -d ${srcdir}/src/kinsol ; then
  KINSOL_ENABLED="yes"
else
  KINSOL_ENABLED="no"
fi
])

# Check if user wants to disable CPODES module
# If not, then make certain source directory actually exists
AC_ARG_ENABLE(cpodes,
[AC_HELP_STRING([--disable-cpodes],[disable configuration of CPODES])],
[
if test "X${enableval}" = "Xno"; then
  CPODES_ENABLED="no"
fi
],
[
if test -d ${srcdir}/src/cpodes ; then
  CPODES_ENABLED="yes"
else
  CPODES_ENABLED="no"
fi
])

# Check if user wants to disable Fortran support (FCMIX components).
AC_ARG_ENABLE([fcmix],
[AC_HELP_STRING([--disable-fcmix], [disable Fortran-C support])],
[
if test "X${enableval}" = "Xno"; then
  FCMIX_ENABLED="no"
fi
],
[
if test "X${CVODE_ENABLED}" = "Xno" && test "X${KINSOL_ENABLED}" = "Xno" && test "X${IDA_ENABLED}" = "Xno"; then
  FCMIX_ENABLED="no"
fi
])

# Check if user wants to disable Lapack support.
AC_ARG_ENABLE([lapack],
[AC_HELP_STRING([--disable-lapack], [disable Lapack support])],
[
if test "X${enableval}" = "Xno"; then
  LAPACK_ENABLED="no"
fi
])

# Check if user wants to disable support for MPI.
# If not, set the default based on whetehr certain source directories exist
AC_ARG_ENABLE([mpi],
[AC_HELP_STRING([--disable-mpi],[disable MPI support])],
[
if test "X${enableval}" = "Xno"; then
  MPI_ENABLED="no"
fi
],
[
if test -d ${srcdir}/src/nvec_par || test -d ${srcdir}/src/nvec_spcpar; then
  MPI_ENABLED="yes"
else
  MPI_ENABLED="no"
fi
])

# Check if user wants to enable all examples.
# Examples are NOT built by default
AC_ARG_ENABLE(examples,
[AC_HELP_STRING([--enable-examples],[enable configuration of examples])],
[
if test "X${enableval}" = "Xno"; then
  EXAMPLES_ENABLED="no"
else
  EXAMPLES_ENABLED="yes"
fi
])

# Fortran examples are enabled only if both FCMIX and EXAMPLES are enabled
if test "X${FCMIX_ENABLED}" = "Xyes" && test "X${EXAMPLES_ENABLED}" = "Xyes"; then
  F77_EXAMPLES_ENABLED="yes"
fi

]) dnl END SUNDIALS_ENABLES


#=================================================================#
#                                                                 #
#                                                                 #
#              C    C O M P I L E R     T E S T S                 #
#                                                                 #
#                                                                 #
#==================================================================


#------------------------------------------------------------------
# CHECK C COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_CC],
[


if test "X${CC}" = "X"; then

  echo ""
  echo "   Unable to find a working C compiler"
  echo ""
  echo "   Try using CC to explicitly specify a C compiler"
  echo ""

  AC_MSG_ERROR([cannot find a C compiler])

else

  SUNDIALS_CC_CHECK

fi

]) dnl END SUNDIALS_SET_CC
 

AC_DEFUN([SUNDIALS_CC_CHECK],
[

# Default is C programming language (initialize language stack)
AC_LANG([C])

AC_ARG_WITH([],[],[])

# Set floating-point precision: single   [C type 'float']
#                               double   [C type 'double'] (default)
#                               extended [C type 'long double']
# Provide variable description templates for config.hin and config.h files
# Required by autoheader utility
AH_TEMPLATE([SUNDIALS_SINGLE_PRECISION],
            [Define SUNDIALS data type 'realtype' as 'float'])
AH_TEMPLATE([SUNDIALS_DOUBLE_PRECISION],
            [Define SUNDIALS data type 'realtype' as 'double'])
AH_TEMPLATE([SUNDIALS_EXTENDED_PRECISION],
            [Define SUNDIALS data type 'realtype' as 'long double'])

AC_MSG_CHECKING([floating-point data type to use])
AC_ARG_WITH(precision,
[AC_HELP_STRING([--with-precision=ARG],
[specify floating-point precision (single/double/extended) [double]])],
[
if test "X${withval}" = "Xsingle"; then
  AC_MSG_RESULT([float])
  AC_DEFINE([SUNDIALS_SINGLE_PRECISION],[1],[])
  FLOAT_TYPE="single"
  PRECISION_LEVEL="#define SUNDIALS_SINGLE_PRECISION 1"
elif test "X${withval}" = "Xdouble"; then
  AC_MSG_RESULT([double])
  AC_DEFINE([SUNDIALS_DOUBLE_PRECISION],[1],[])
  FLOAT_TYPE="double"
  PRECISION_LEVEL="#define SUNDIALS_DOUBLE_PRECISION 1"
elif test "X${withval}" = "Xextended"; then
  AC_MSG_RESULT([long double])
  AC_DEFINE([SUNDIALS_EXTENDED_PRECISION],[1],[])
  FLOAT_TYPE="extended"
  PRECISION_LEVEL="#define SUNDIALS_EXTENDED_PRECISION 1"
else
  AC_MSG_ERROR([invalid input])
fi
],
[
# Use 'double' by default
AC_MSG_RESULT([double])
AC_DEFINE([SUNDIALS_DOUBLE_PRECISION],[1],[])
FLOAT_TYPE="double"
PRECISION_LEVEL="#define SUNDIALS_DOUBLE_PRECISION 1"
])

AC_ARG_WITH([],[ ],[])

# Overwrite CFLAGS
AC_MSG_CHECKING([for C compiler flags])
AC_ARG_WITH(cflags,
[AC_HELP_STRING([--with-cflags=ARG],[specify C compiler flags (CFLAGS will be overridden)])],
[
AC_MSG_RESULT([${withval}])
CFLAGS="${withval}"
],
[
AC_MSG_RESULT([none])
])

# Set CPP to command that runs C preprocessor
AC_PROG_CPP

# Overwrite CPPFLAGS
AC_MSG_CHECKING([for C/C++ preprocessor flags])
AC_ARG_WITH(cppflags,
[AC_HELP_STRING([--with-cppflags=ARG],[specify C/C++ preprocessor flags (CPPFLAGS will be overridden)])],
[
AC_MSG_RESULT([${withval}])
CPPFLAGS="${withval}"
],
[
AC_MSG_RESULT([none])
])

# Overwrite LDFLAGS
AC_MSG_CHECKING([for linker flags])
AC_ARG_WITH(ldflags,
[AC_HELP_STRING([--with-ldflags=ARG],[specify linker flags (LDFLAGS will be overridden)])],
[
AC_MSG_RESULT([${withval}])
LDFLAGS="${withval}"
],
[
AC_MSG_RESULT([none])
])

# Add any additional libraries
AC_MSG_CHECKING([for extra libraries])
AC_ARG_WITH(libs,
[AC_HELP_STRING([--with-libs=ARG],[add extra libraries])],
[
AC_MSG_RESULT([${withval}])
if test "X${LIBS}" = "X"; then
  LIBS="${withval}"
else
  LIBS="${LIBS} ${withval}"
fi
],
[
AC_MSG_RESULT([none])
])

# Defines STDC_HEADERS if the following header files are found: stdlib.h,
# stdarg.h, string.h, and float.h
# We really only need stdlib.h and float.h
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h float.h math.h])

# Set flag indicating if generic function names should be used
# Provide variable description template for config.hin and config.h files
# Required by autoheader utility
AH_TEMPLATE([SUNDIALS_USE_GENERIC_MATH],
            [Use generic math functions])

# Check if math library contains abs(), fabs(), pow(), and sqrt() functions (required)
# May update LIBS (meaning add additional library, namely libm)
MATH_FABS_OK="yes"
MATH_POW_OK="yes"
MATH_SQRT_OK="yes"
# Save copy of LIBS variable and unset LIBS
SAVED_LIBS="${LIBS}"
LIBS=""
# The abs routine is defined for an integer argument, so check for it regardless of
# the level of precision chosen
AC_CHECK_LIB([m],abs,[],[AC_MSG_ERROR([cannot find abs function])])
TEMP_MATH_LIB="${LIBS}"
LIBS=""
# Check for single-precision math routines
if test "X${FLOAT_TYPE}" = "Xsingle"; then
  AC_CHECK_LIB([m],fabsf,[],[MATH_FABS_OK="no"])
  AC_CHECK_LIB([m],powf,[],[MATH_POW_OK="no"])
  AC_CHECK_LIB([m],sqrtf,[],[MATH_SQRT_OK="no"])
# Check for extended-precision math routines
elif test "X${FLOAT_TYPE}" = "Xextended"; then
  AC_CHECK_LIB([m],fabsl,[],[MATH_FABS_OK="no"])
  AC_CHECK_LIB([m],powl,[],[MATH_POW_OK="no"])
  AC_CHECK_LIB([m],sqrtl,[],[MATH_SQRT_OK="no"])
# Check for (generic) double-precision math routines
elif test "X${FLOAT_TYPE}" = "Xdouble"; then
  AC_CHECK_LIB([m],fabs,[],[AC_MSG_ERROR([cannot find fabs function])])
  AC_CHECK_LIB([m],pow,[],[AC_MSG_ERROR([cannot find pow function])])
  AC_CHECK_LIB([m],sqrt,[],[AC_MSG_ERROR([cannot find sqrt function])])
fi
# If cannot find precision-specific implementations, then check for generic versions
if test "X${MATH_FABS_OK}" = "Xno" || test "X${MATH_POW_OK}" = "Xno" || test "X${MATH_SQRT_OK}" = "Xno"; then
  AC_CHECK_LIB([m],fabs,[],[AC_MSG_ERROR([cannot find fabs function])])
  AC_CHECK_LIB([m],pow,[],[AC_MSG_ERROR([cannot find pow function])])
  AC_CHECK_LIB([m],sqrt,[],[AC_MSG_ERROR([cannot find sqrt function])])
  # If all generic math routines are available, then set SUNDIALS_USE_GENERIC_MATH flag
  # for use by sundials_math.c file (preprocessor macros)
  AC_DEFINE([SUNDIALS_USE_GENERIC_MATH],[1],[])
  GENERIC_MATH_LIB="#define SUNDIALS_USE_GENERIC_MATH 1"
# If found all precision-specific routines, then set SUNDIALS_USE_GENERIC_MATH only if
# building SUNDIALS libraries with double-precision
else
  if test "X${FLOAT_TYPE}" = "Xdouble"; then
    AC_DEFINE([SUNDIALS_USE_GENERIC_MATH],[1],[])
    GENERIC_MATH_LIB="#define SUNDIALS_USE_GENERIC_MATH 1"
  else
    AC_DEFINE([SUNDIALS_USE_GENERIC_MATH],[0],[])
  fi
fi

# Add math library to LIBS environment variable
LIBS="${TEMP_MATH_LIB}"
AC_MSG_CHECKING([for additional required C libraries])
if test "X${LIBS}" = "X"; then
  if test "X${SAVED_LIBS}" = "X"; then
    LIBS=""
  else
    LIBS="${SAVED_LIBS}"
  fi
  AC_MSG_RESULT([none])
else
  AC_MSG_RESULT([${LIBS}])
  if test "X${SAVED_LIBS}" = "X"; then
    LIBS="${LIBS}"
  else
    LIBS="${LIBS} ${SAVED_LIBS}"
  fi
fi

# Check sizeof(int) - used to modify Fortran examples
AC_CHECK_SIZEOF(int)

# Check sizeof(long int) - used to modify Fortran examples
AC_CHECK_SIZEOF(long int)

# Check sizeof(realtype), where realtype is either float, double
# or long double - used to modify Fortran examples
if test "X${FLOAT_TYPE}" = "Xsingle"; then
  AC_CHECK_SIZEOF(float)
elif test "X${FLOAT_TYPE}" = "Xdouble"; then
  AC_CHECK_SIZEOF(double)
elif test "X${FLOAT_TYPE}" = "Xextended"; then
  AC_CHECK_SIZEOF(long double)
fi

# Defines EGREP and exports via AC_SUBST - used by FCMIX Makefile's
AC_PROG_EGREP

# Defines FGREP and exports via AC_SUBST - used by FCMIX Makefile's
AC_PROG_FGREP

# Check if CC is a C++ compiler
# Note: If CC is a C++ compiler and MPI is enabled, then we will
# check for "mpiCC" instead of "mpicc" if an MPI compiler was NOT specified
SUNDIALS_CPLUSPLUS_CHECK([${CC}])

]) dnl END SUNDIALS_SET_CC

#------------------------------------------------------------------
# CHECK IF COMPILER IS A C++ COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CPLUSPLUS_CHECK],
[

# Rename argument
COMP_NAME="$1"

# Update the language stack
AC_LANG_PUSH([C])

# Check if using a C++ compiler
AC_MSG_CHECKING([if ${COMP_NAME} is a C++ compiler])
AC_RUN_IFELSE(
[AC_LANG_PROGRAM([[]],
[[
#ifdef __cplusplus
  return(0);
#else
  return(1);
#endif
]])],
[
AC_MSG_RESULT([yes])
# COMP_NAME is a C++ compiler
USING_CPLUSPLUS_COMP="yes"
],
[
AC_MSG_RESULT([no])
# COMP_NAMPE is NOT a C++ compiler
USING_CPLUSPLUS_COMP="no"
])

# Revert back to previous language
AC_LANG_POP([C])

]) dnl END SUNDIALS_CPLUSPLUS_CHECK




#=================================================================#
#                                                                 #
#                                                                 #
#             F O R T R A N     S U P P O R T                     #
#                                                                 #
#                                                                 #
#==================================================================



#------------------------------------------------------------------
# FORTRAN SUPPORT
#
# Fortran support is required if FCMIX is enabled OR if LAPACK
# is enabled. In either case, we need a working F77 compiler in
# order to determine the Fortran name-mangling scheme.
#
# If we do need Fortran support, we first find and test a F77
# compiler, determine the mangling scheme, then we find the 
# libraries required to link C and Fortran.
#
# Throughout this function we use the control variable F77_OK
# which was initialized to "no".
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_SUPPORT],
[

F77_OK="yes"

# Look for a F77 compiler
# If unsuccessful, disable all Fortran support

AC_PROG_F77(f77 g77)

if test "X${F77}" = "X"; then

  F77_OK="no"
  SUNDIALS_WARN_FLAG="yes"

  echo ""
  echo "   Unable to find a working Fortran compiler"
  echo ""
  echo "   Try using F77 to explicitly specify a C compiler"
  echo ""
  if test "X${FCMIX_ENABLED}" = "Xyes"; then
    echo "   Disabling compilation of Fortran-C interfaces..."
  fi
  if test "X${LAPACK_ENABLED}" = "Xyes"; then
    echo "   Disabling compilation of Blas/Lapack interfaces..."
  fi
  echo ""

  FCMIX_ENABLED="no"
  LAPACK_ENABLED="no"
  F77_EXAMPLES_ENABLED="no"

fi

# Check Fortran compiler
# If unsuccessful, disable all Fortran support

if test "X${F77_OK}" = "Xyes"; then

  SUNDIALS_F77_CHECK

  if test "X${F77_OK}" = "Xno"; then

    SUNDIALS_WARN_FLAG="yes"

    echo ""
    echo "   Unable to compile test program using given Fortran compiler."
    echo ""
    if test "X${FCMIX_ENABLED}" = "Xyes"; then
      echo "   Disabling compilation of Fortran-C interfaces..."
    fi
    if test "X${LAPACK_ENABLED}" = "Xyes"; then
      echo "   Disabling compilation of Blas/Lapack interfaces..."
    fi
    echo ""

    FCMIX_ENABLED="no"
    LAPACK_ENABLED="no"
    F77_EXAMPLES_ENABLED="no"
    
  fi

fi


# Determine the Fortran name-mangling scheme
# If successfull, provide variable description templates for config.hin 
# and config.h files required by autoheader utility
# Otherwise, disable all Fortran support.

if test "X${F77_OK}" = "Xyes"; then

  SUNDIALS_F77_NAME_MANGLING

  AH_TEMPLATE([SUNDIALS_F77_FUNC], [FCMIX: Define name-mangling macro for C identifiers])
  AH_TEMPLATE([SUNDIALS_F77_FUNC_], [FCMIX: Define name-mangling macro for C identifiers with underscores])

  if test "X${F77_OK}" = "Xno"; then

    SUNDIALS_WARN_FLAG="yes"

    echo ""
    echo "   Unable to determine Fortran name-mangling scheme."
    echo ""
    if test "X${FCMIX_ENABLED}" = "Xyes"; then
      echo "   Disabling compilation of Fortran-C interfaces..."
    fi
    if test "X${LAPACK_ENABLED}" = "Xyes"; then
      echo "   Disabling compilation of Blas/Lapack interfaces..."
    fi
    echo ""

    F77_EXAMPLES_ENABLED="no"
    FCMIX_ENABLED="no"
    LAPACK_ENABLED="no"

  fi

fi


# If LAPACK is enabled, determine the proper library linkage
# If successful, set the libaries
# Otherwise, disable all Blas/Lapack support.

if test "X${LAPACK_ENABLED}" = "Xyes" && test "X${F77_OK}" = "Xyes"; then


  SUNDIALS_F77_LAPACK_SET

  if test "X${LAPACK_OK}" = "Xyes"; then

    AC_MSG_CHECKING([for Blas/Lapack library linkage])
    BLAS_LAPACK_LIBS="${LAPACK_LIBS} ${BLAS_LIBS} ${LIBS} ${FLIBS}"
    AC_MSG_RESULT([${LAPACK_LIBS} ${BLAS_LIBS}])

  else

    SUNDIALS_WARN_FLAG="yes"
    AC_MSG_CHECKING([for Blas/Lapack library linkage])
    AC_MSG_RESULT("no")
    echo ""
    echo "   Unable to determine Blas/Lapack library linkage."
    echo ""
    echo "   Try using --with-blas and --with-lapack."
    echo ""
    echo "   Disabling compilation of Blas/Lapack interfaces..."

    LAPACK_ENABLED="no"

  fi

fi


# Set the macro BLAS_LAPACK_MACRO for expansion in sundials_config.h

AH_TEMPLATE([SUNDIALS_BLAS_LAPACK], [Availability of Blas/Lapack libraries])

if test "X${LAPACK_ENABLED}" = "Xyes"; then
  AC_DEFINE([SUNDIALS_BLAS_LAPACK],[1],[])
  BLAS_LAPACK_MACRO="#define SUNDIALS_BLAS_LAPACK 1"
else
  AC_DEFINE([SUNDIALS_BLAS_LAPACK],[0],[])
  BLAS_LAPACK_MACRO="#define SUNDIALS_BLAS_LAPACK 0"
fi

]) dnl SUNDIALS_F77_SUPPORT

#------------------------------------------------------------------
# CHECK FORTRAN COMPILER
#
# Test the Fortran compiler by attempting to compile and link a 
# simple Fortran program. If the test succeeds, set F77_OK=yes.
# If the test fails, set F77_OK="no"
#
# Finally, check if we must use a Fortran compiler to link the
# Fortran codes (default is to use CC).
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_CHECK],
[

AC_LANG_PUSH([Fortran 77])

# Add any additional FFLAGS
AC_MSG_CHECKING([for extra Fortran compiler flags])
AC_ARG_WITH(fflags,
[AC_HELP_STRING([--with-fflags=ARG],[add extra Fortran compiler flags])],
[
AC_MSG_RESULT([${withval}])
FFLAGS="${FFLAGS} ${withval}"
],
[
AC_MSG_RESULT([none])
])
  
# Add any required linker flags to FLIBS
# Note: if FLIBS is defined, it is left unchanged
AC_F77_LIBRARY_LDFLAGS

# Try to compile a simple Fortran program (no linking)
AC_COMPILE_IFELSE(
[AC_LANG_SOURCE(
[[
	SUBROUTINE SUNDIALS()
	RETURN
	END
]])],
[F77_OK="yes"],
[F77_OK="no"])

# If CC is a C++ compiler (decided in SUNDIALS_CPLUSPLUS_CHECK), we must use 
# it to link the Fortran examples. In this case, test if that is successful.
# Otherwise, simply use F77 as the linker

if test "X${F77_OK}" = "Xyes"; then
  AC_MSG_CHECKING([which linker to use])
  if test "X${USING_CPLUSPLUS_COMP}" = "Xyes"; then
    SUNDIALS_F77_LNKR_CHECK
  else
    F77_LNKR="${F77}"
  fi
  AC_MSG_RESULT([${F77_LNKR}])
fi

# Reset language (remove 'Fortran 77' from stack)
AC_LANG_POP([Fortran 77])

]) dnl END SUNDIALS_SET_F77


#------------------------------------------------------------------
# F77 LINKER CHECK
# Check if the C++ compiler CC can be used to link a Fortran program.
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_LNKR_CHECK],
[

F77_LNKR_CHECK_OK="no"

# Compile simple Fortran example, but do NOT link
# Note: result stored as conftest.${ac_objext}
AC_COMPILE_IFELSE(
[AC_LANG_SOURCE(
[[
	PROGRAM SUNDIALS
	WRITE(*,*)'TEST'
	END
]])],
[

# Temporarily reset LIBS environment variable to perform test
SAVED_LIBS="${LIBS}"
LIBS="${LIBS} ${FLIBS}"

# Switch working language to C for next test
AC_LANG_PUSH([C])

# Check if CC can link Fortran example
# Note: AC_LINKONLY_IFELSE is a custom macro (modifications made to
# general.m4 and c.m4) (see config/cust_general.m4 and config/mod_c.m4)
AC_LINKONLY_IFELSE([],[F77_LNKR_CHECK_OK="yes"],[F77_LNKR_CHECK_OK="no"])

# Revert back to previous language (Fortran 77)
AC_LANG_POP([C])

# Set LIBS environment variable back to original value
LIBS="${SAVED_LIBS}"

])

# If either the compilation or the linking failed, we should
# disable building the Fortran examples
# For now, use F77 as the linker...
if test "X${F77_LNKR_CHECK_OK}" = "Xyes"; then
  F77_LNKR="${CC}"
else
  F77_LNKR="${F77}"
fi

]) dnl SUNDIALS_F77_LNKR_CHECK


#------------------------------------------------------------------
# DETERMINE FORTRAN NAME-MANGLING SCHEME
#
# Compiling a simple Fortran example and link it using a C compiler.
# Interpret results to infer name-mangling scheme.
#------------------------------------------------------------------


AC_DEFUN([SUNDIALS_F77_NAME_MANGLING],
[

AC_LANG_PUSH([Fortran 77])

# (1) Compile a dummy Fortran subroutine named SUNDIALS

FNAME_STATUS="none"

AC_COMPILE_IFELSE(
  [AC_LANG_SOURCE(
  [[
	SUBROUTINE SUNDIALS()
	RETURN
	END
  ]])],
  [

  mv conftest.${ac_objext} f77_wrapper_check.${ac_objext}

  # Temporarily reset LIBS environment variable to perform test
  SAVED_LIBS="${LIBS}"
  LIBS="f77_wrapper_check.${ac_objext} ${LIBS} ${FLIBS}"

  AC_LANG_PUSH([C])

  for i in "sundials" "SUNDIALS"
  do
    for j in "" "_" "__"
    do
      F77_MANGLED_NAME="${i}${j}"
      AC_LINK_IFELSE([AC_LANG_CALL([],[${F77_MANGLED_NAME}])],[FNAME_STATUS="set" ; break 2])
    done
  done

  AC_LANG_POP([C])

  # If test succeeded, then set the F77_MANGLE_MACRO1 macro

  if test "X${FNAME_STATUS}" = "Xset"; then

    if test "X${i}" = "Xsundials"; then

      FNAME_MSG="lower case "

      if test "X${j}" = "X"; then
        FNAME_MSG="${FNAME_MSG} + no underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC(name,NAME)],[name],[])
        F77_MANGLE_MACRO1="#define SUNDIALS_F77_FUNC(name,NAME) name"
        dgemm="dgemm"
        dgetrf="dgetrf"
      elif test "X${j}" = "X_"; then
        FNAME_MSG="${FNAME_MSG} + one underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC(name,NAME)],[name ## _],[])
        F77_MANGLE_MACRO1="#define SUNDIALS_F77_FUNC(name,NAME) name ## _"
        dgemm="dgemm_"
        dgetrf="dgetrf_"
      else
        FNAME_MSG="${FNAME_MSG} + two underscores"
        AC_DEFINE([SUNDIALS_F77_FUNC(name,NAME)],[name ## __],[])
        F77_MANGLE_MACRO1="#define SUNDIALS_F77_FUNC(name,NAME) name ## __"
        dgemm="dgemm__"
        dgetrf="dgetrf__"
      fi

    else

      FNAME_MSG="upper case "

      if test "X${j}" = "X"; then
        FNAME_MSG="${FNAME_MSG} + no underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC(name,NAME)],[name],[])
        F77_MANGLE_MACRO1="#define SUNDIALS_F77_FUNC(name,NAME) NAME"
        dgemm="DGEMM"
        dgetrf="DGETRF"
      elif test "X${j}" = "X_"; then
        FNAME_MSG="${FNAME_MSG} + one underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC(name,NAME)],[name ## _],[])
        F77_MANGLE_MACRO1="#define SUNDIALS_F77_FUNC(name,NAME) NAME ## _"
        dgemm="DGEMM_"
        dgetrf="DGETRF_"
      else
        FNAME_MSG="${FNAME_MSG} + two underscores"
        AC_DEFINE([SUNDIALS_F77_FUNC(name,NAME)],[name ## __],[])
        F77_MANGLE_MACRO1="#define SUNDIALS_F77_FUNC(name,NAME) NAME ## __"
        dgemm="DGEMM__"
        dgetrf="DGETRF__"
      fi

    fi

    AC_MSG_CHECKING([for Fortran name-mangling scheme of C identifiers])
    AC_MSG_RESULT([${FNAME_MSG}])

  else

    F77_OK="no"

  fi

  # Set LIBS environment variable back to original value
  LIBS="${SAVED_LIBS}"

  ])

# Remove temporary file
rm -f f77_wrapper_check.${ac_objext}


# (2) Compile a dummy Fortran subroutine named SUN_DIALS

FNAME_STATUS="none"

AC_COMPILE_IFELSE(
  [AC_LANG_SOURCE(
  [[
	SUBROUTINE SUN_DIALS()
	RETURN
	END
  ]])],
  [

  mv conftest.${ac_objext} f77_wrapper_check.${ac_objext}

  # Temporarily reset LIBS environment variable to perform test
  SAVED_LIBS="${LIBS}"
  LIBS="f77_wrapper_check.${ac_objext} ${LIBS} ${FLIBS}"

  AC_LANG_PUSH([C])

  for i in "sun_dials" "SUN_DIALS"
  do
    for j in "" "_" "__"
    do
      F77_MANGLED_NAME="${i}${j}"
      AC_LINK_IFELSE([AC_LANG_CALL([],[${F77_MANGLED_NAME}])],[FNAME_STATUS="set" ; break 2])
    done
  done

  AC_LANG_POP([C])

  # If test succeeded, then set the F77_MANGLE_MACRO2 macro

  if test "X${FNAME_STATUS}" = "Xset"; then

    if test "X${i}" = "Xsun_dials"; then

      FNAME_MSG="lower case "

      if test "X${j}" = "X"; then
        FNAME_MSG="${FNAME_MSG} + no underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC_(name,NAME)],[name],[])
        F77_MANGLE_MACRO2="#define SUNDIALS_F77_FUNC_(name,NAME) name"
      elif test "X${j}" = "X_"; then
        FNAME_MSG="${FNAME_MSG} + one underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC_(name,NAME)],[name ## _],[])
        F77_MANGLE_MACRO2="#define SUNDIALS_F77_FUNC_(name,NAME) name ## _"
      else
        FNAME_MSG="${FNAME_MSG} + two underscores"
        AC_DEFINE([SUNDIALS_F77_FUNC_(name,NAME)],[name ## __],[])
        F77_MANGLE_MACRO2="#define SUNDIALS_F77_FUNC_(name,NAME) name ## __"
      fi

    else

      FNAME_MSG="upper case "

      if test "X${j}" = "X"; then
        FNAME_MSG="${FNAME_MSG} + no underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC_(name,NAME)],[name],[])
        F77_MANGLE_MACRO2="#define SUNDIALS_F77_FUNC_(name,NAME) NAME"
      elif test "X${j}" = "X_"; then
        FNAME_MSG="${FNAME_MSG} + one underscore"
        AC_DEFINE([SUNDIALS_F77_FUNC_(name,NAME)],[name ## _],[])
        F77_MANGLE_MACRO2="#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## _"
      else
        FNAME_MSG="${FNAME_MSG} + two underscores"
        AC_DEFINE([SUNDIALS_F77_FUNC_(name,NAME)],[name ## __],[])
        F77_MANGLE_MACRO2="#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## __"
      fi

    fi

    AC_MSG_CHECKING([for Fortran name-mangling scheme of C identifiers with underscores])
    AC_MSG_RESULT([${FNAME_MSG}])

  else

    F77_OK="no"

  fi

  # Set LIBS environment variable back to original value
  LIBS="${SAVED_LIBS}"

  ])

# Remove temporary file
rm -f f77_wrapper_check.${ac_objext}


AC_LANG_POP([Fortran 77])

]) dnl END SUNDIALS_SET_FNAME


#------------------------------------------------------------------
# DETERMINE BLAS/LAPACK LIBRARY LINKAGE SCHEME
#
# If successful, this function sets LAPACK_OK="yes".
# Otherwise, it sets LAPACK_OK="no" 
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_LAPACK_SET],
[

# Check if the user specifies Blas libraries
AC_ARG_WITH(blas,
[AC_HELP_STRING([--with-blas=ARG],[specify Blas library])],
[
  case $withval in
    -* | */* | *.a | *.so | *.so.* | *.o) 
        BLAS_LIBS="$withval" 
        ;;
    *) 
        BLAS_LIBS="-l$withval" 
        ;;
  esac
])

# Check if the user specifies Lapack libraries
AC_ARG_WITH(lapack,
[AC_HELP_STRING([--with-lapack=ARG],[specify Lapack library])],
[
  case $withval in
    -* | */* | *.a | *.so | *.so.* | *.o) 
        LAPACK_LIBS="$withval" 
        ;;
    *) 
        LAPACK_LIBS="-l$withval" 
        ;;
  esac
])

acx_blas_ok=no
acx_lapack_ok=no

# BLAS_LIBS
# ---------

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test "x$BLAS_LIBS" != x; then
  save_LIBS="$LIBS"
  LIBS="$BLAS_LIBS $LIBS"
  AC_MSG_CHECKING([aha for $dgemm in $BLAS_LIBS])
  AC_TRY_LINK_FUNC($dgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
  AC_MSG_RESULT($acx_blas_ok)
  LIBS="$save_LIBS"
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS"
  AC_CHECK_FUNC($dgemm, [acx_blas_ok=yes])
  LIBS="$save_LIBS"
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(cxml, $dgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(dxml, $dgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
  if test "x$GCC" != xyes; then # only works with Sun CC
    AC_CHECK_LIB(sunmath, acosp,
                 [AC_CHECK_LIB(sunperf, $dgemm,
                               [BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                acx_blas_ok=yes],[],[-lsunmath])])
  fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(scs, $dgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(complib.sgimath, $dgemm,
               [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
 AC_CHECK_LIB(blas, $dgemm,
              [AC_CHECK_LIB(essl, $dgemm,
                        [acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
                        [], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(blas, $dgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

LIBS="$acx_blas_save_LIBS"

# LAPACK
# ------

# If we didn't find a Blas implementation, disable tests for Lapack
if test $acx_blas_ok = no; then
  acx_lapack_ok=disabled
fi

# Check LAPACK_LIBS environment variable
if test $acx_lapack_ok = no; then
if test "x$LAPACK_LIBS" != x; then
  save_LIBS="$LIBS"; 
  LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
  AC_MSG_CHECKING([for $dgetrf in $LAPACK_LIBS])
  AC_TRY_LINK_FUNC($dgetrf, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
  AC_MSG_RESULT($acx_lapack_ok)
  LIBS="$save_LIBS"
  if test acx_lapack_ok = no; then
     LAPACK_LIBS=""
  fi
fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
  save_LIBS="$LIBS"
  LIBS="$LIBS $BLAS_LIBS $FLIBS"
  AC_CHECK_FUNC($dgetrf, [acx_lapack_ok=yes])
  LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
  if test $acx_lapack_ok = no; then
    save_LIBS="$LIBS"
      LIBS="$BLAS_LIBS $LIBS"
      AC_CHECK_LIB($lapack, $dgetrf,
                   [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
      LIBS="$save_LIBS"
  fi
done

# If we have both libraries, set LAPACK_OK to yes
# -----------------------------------------------------

if test $acx_blas_ok = yes && test $acx_lapack_ok = yes; then
  LAPACK_OK="yes"
else
  LAPACK_OK="no"
fi

]) dnl SUNDIALS_F77_LAPACK_SET





#=================================================================#
#                                                                 #
#                                                                 #
#           P A R A L L E L     S U P P O R T                     #
#                                                                 #
#                                                                 #
#==================================================================


#------------------------------------------------------------------
# CHECK MPI-C COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_MPICC],
[

AC_ARG_WITH([],[    ],[])

# MPI root directory
AC_ARG_WITH(mpi-root,
[AC_HELP_STRING([--with-mpi-root=MPIROOT],[use MPI root directory])],
[
MPI_ROOT_DIR="${withval}"
],
[
MPI_ROOT_DIR=""
])

# MPI include directory
AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@])],
[
MPI_INC_DIR="${withval}"
],
[
MPI_INC_DIR=""
])

# MPI library directory
AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@])],
[
MPI_LIB_DIR="${withval}"
],
[
MPI_LIB_DIR=""
])

# MPI libraries
AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs=ARG],[MPI libraries])],
[
MPI_LIBS="${withval}"
],
[
MPI_LIBS=""
])

# MPI flags
AC_ARG_WITH(mpi-flags,
[AC_HELP_STRING([--with-mpi-flags=ARG],[MPI-specific flags])],
[
MPI_FLAGS="${withval}"
MPI_FLAGS_OK="yes"
],
[
MPI_FLAGS=""
MPI_FLAGS_OK="no"
])

# MPI-C compiler
MPICC_COMP_GIVEN="yes"
AC_MSG_CHECKING([if using MPI-C script])
AC_ARG_WITH(mpicc,
[AC_HELP_STRING([--with-mpicc[[[[=ARG]]]]],[specify MPI-C compiler to use @<:@mpicc@:>@])],
[
if test "X${withval}" = "Xno"; then
  USE_MPICC_SCRIPT="no"
else
  USE_MPICC_SCRIPT="yes"
  MPICC_COMP="${withval}"
fi
],
[
  USE_MPICC_SCRIPT="yes"
  MPICC_COMP="mpicc"
  MPICC_COMP_GIVEN="no"
])
AC_MSG_RESULT([${USE_MPICC_SCRIPT}])

# If CC is a C++ compiler, then we certainly do NOT want to use an MPI-C script
# Note: USING_CPLUSPLUS_COMP was defined by a call to SUNDIALS_CPLUSPLUS_CHECK
# in SUNDIALS_SET_CC
# Note: If the user specified an MPI-C script, then we will NOT do anything for now
if test "X${MPICC_COMP_GIVEN}" = "Xno" && test "X${USING_CPLUSPLUS_COMP}" = "Xyes"; then
  MPICC_COMP="mpiCC"
fi

# Check MPI-C compiler (either MPI compiler script or regular C compiler)
if test "X${USE_MPICC_SCRIPT}" = "Xyes"; then
  SUNDIALS_CHECK_MPICC
else
  MPICC_COMP="${CC}"
  MPICC="${CC}"
  SUNDIALS_CC_WITH_MPI_CHECK
fi

]) dnl END SUNDIALS_SET_MPICC

#------------------------------------------------------------------
# TEST MPI-C COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_MPICC],
[

# Test MPI-C compiler (meaning test MPICC_COMP)
# Check if MPI-C compiler can be found

AC_MSG_CHECKING([if absolute path to ${MPICC_COMP} was given])

# CASE 1: MPICC_COMP was found (cannot check if executable because the
# "-x" flag is NOT portable)
if test -f ${MPICC_COMP} ; then

  AC_MSG_RESULT([yes])
  MPICC_COMP_EXISTS="yes"
  # Determine MPI_INC_DIR and MPI_LIB_DIR for use by Makefile
  MPI_BASE_DIR=`AS_DIRNAME(["${MPICC_COMP}"])`
  TMP_MPI_INC_DIR="${MPI_BASE_DIR}/../include"
  TMP_MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"

# CASE 2: MPICC_COMP could NOT be found and MPI_ROOT_DIR was NOT specified,
# so search in PATH
else

  AC_MSG_RESULT([no])

  if test "X${MPI_ROOT_DIR}" = "X"; then
    # Try to find location of executable (perhaps directory was entered
    # incorrectly)
    TEMP_MPICC_COMP=`basename "${MPICC_COMP}"`
    AC_PATH_PROG([MPICC_COMP],[${TEMP_MPICC_COMP}],[none])
    # Cannot find executable in PATH
    if test "X${MPICC_COMP}" = "Xnone"; then
      MPICC_COMP_EXISTS="no"
      MPICC_COMP=""
    # Found executable and set MPICC_COMP to absolute pathname
    else
      MPICC_COMP_EXISTS="yes"
      MPI_BASE_DIR=`AS_DIRNAME(["${MPICC_COMP}"])`
      TMP_MPI_INC_DIR="${MPI_BASE_DIR}/../include"
      TMP_MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
    fi

  # CASE 3: MPICC_COMP could NOT be found, but MPI_ROOT_DIR was specified
  else

    AC_MSG_CHECKING([if ${MPICC_COMP} exists in ${MPI_ROOT_DIR}/bin])
    # MPICC_COMP should really only contain an executable name
    # Found location of MPICC_COMP
    if test -f ${MPI_ROOT_DIR}/bin/${MPICC_COMP} ; then
      AC_MSG_RESULT([yes])
      MPICC_COMP_EXISTS="yes"
      MPICC_COMP="${MPI_ROOT_DIR}/bin/${MPICC_COMP}"
      TMP_MPI_INC_DIR="${MPI_ROOT_DIR}/include"
      TMP_MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    # Could NOT find MPICC_COMP anywhere
    else
      AC_MSG_RESULT([no])
      MPICC_COMP_EXISTS="no"
      MPICC_COMP=""
    fi

  fi

fi

# If MPICC_COMP exists, set MPICC and (conditionally) set MPI_INC_DIR
# and MPI_LIB_DIR so that we do not end up with empty -I options.
# Otherwise, issue warning message
if test "X${MPICC_COMP_EXISTS}" = "Xyes"; then

  MPICC="${MPICC_COMP}"
  MPI_C_COMP_OK="yes"

  # If MPI_INC_DIR is empty, set it to TMP_MPI_INC_DIR
  if test "X${MPI_INC_DIR}" = "X"; then
    MPI_INC_DIR="$TMP_MPI_INC_DIR"
  fi

  # If MPI_LIB_DIR is empty, set it to TMP_MPI_LIB_DIR
  if test "X${MPI_LIB_DIR}" = "X"; then
    MPI_LIB_DIR="$TMP_MPI_LIB_DIR"
  fi

else

  AC_MSG_WARN([cannot find MPI-C compiler])
  echo ""
  echo "   Unable to find a functional MPI-C compiler."
  echo ""
  echo "   Try using --with-mpicc to specify a MPI-C compiler script,"
  echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
  echo "   to specify the locations of all relevant MPI files, or"
  echo "   --with-mpi-root to specify the base installation directory"
  echo "   of the MPI implementation to be used."
  echo ""
  echo "   Disabling the parallel NVECTOR module and all parallel examples..."
  echo ""
  MPICC=""
  MPI_C_COMP_OK="no"
  SUNDIALS_WARN_FLAG="yes"

fi

]) dnl END SUNDIALS_CHECK_MPICC

#------------------------------------------------------------------
# TEST C COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CC_WITH_MPI_CHECK],
[

# Test if we can compile MPI programs using the CC compiler
# and current MPI settings

AC_MSG_NOTICE([Testing CC with MPI settings])

# Save copies of CPPFLAGS, LDFLAGS and LIBS (preserve information)
# Temporarily overwritten so we can test MPI implementation
SAVED_CPPFLAGS="${CPPFLAGS}"
SAVED_LDFLAGS="${LDFLAGS}"
SAVED_LIBS="${LIBS}"

# Determine location of MPI header files (find MPI include directory)
MPI_EXISTS="yes"

AC_MSG_CHECKING([for location of MPI implementation])

# If MPI include directory was NOT explicitly specified, check if MPI root
# directory was given by user
if test "X${MPI_INC_DIR}" = "X"; then
  # If MPI root directory was NOT given so issue a warning message
  if test "X${MPI_ROOT_DIR}" = "X"; then
    AC_MSG_RESULT([not found])
    MPI_EXISTS="no"
    AC_MSG_WARN([cannot find MPI implementation files])
    echo ""
    echo "   Unable to find MPI implementation files."
    echo ""
    echo "   Try using --with-mpicc to specify a MPI-C compiler script,"
    echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
    echo "   to specify the locations of all relevant MPI files, or"
    echo "   --with-mpi-root to specify the base installation directory"
    echo "   of the MPI implementation to be used."
    echo ""
    echo "   Disabling the parallel NVECTOR module and all parallel examples..."
    echo ""
    SUNDIALS_WARN_FLAG="yes"
  # MPI root directory was given so set MPI_INC_DIR accordingly
  # Update CPPFLAGS
  else
    MPI_INC_DIR="${MPI_ROOT_DIR}/include"
    AC_MSG_RESULT([${MPI_INC_DIR}])
    if test "X${CPPFLAGS}" = "X"; then
      CPPFLAGS="-I${MPI_INC_DIR}"
    else
      CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
    fi
    # Add MPI_FLAGS if non-empty
    if test "X${MPI_FLAGS}" = "X"; then
      CPPFLAGS="${CPPFLAGS}"
    else
      CPPFLAGS="${CPPFLAGS} ${MPI_FLAGS}"
    fi
  fi
# MPI include directory was specified so update CPPFLAGS
else
  AC_MSG_RESULT([${MPI_INC_DIR}])
  if test "X${CPPFLAGS}" = "X"; then
    CPPFLAGS="-I${MPI_INC_DIR}"
  else
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
  fi
  # Add MPI_FLAGS if non-empty
  if test "X${MPI_FLAGS}" = "X"; then
    CPPFLAGS="${CPPFLAGS}"
  else
    CPPFLAGS="${CPPFLAGS} ${MPI_FLAGS}"
  fi
fi

# Only continue if found an MPI implementation
if test "X${MPI_EXISTS}" = "Xyes"; then

  AC_MSG_CHECKING([for location of MPI libraries])

  # Determine location of MPI libraries
  # MPI library directory was NOT specified by user so set based upon MPI_ROOT_DIR
  # Update LDFLAGS
  if test "X${MPI_LIB_DIR}" = "X"; then
    MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    AC_MSG_RESULT([${MPI_LIB_DIR}])
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  # MPI library directory was specified so update LDFLAGS
  else
    AC_MSG_RESULT([${MPI_LIB_DIR}])
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  fi

  # Check if user specified which MPI libraries must be included
  # If no libraries are given, then issue a warning message
  AC_MSG_CHECKING([for MPI libraries])
  if test "X${MPI_LIBS}" = "X"; then
    AC_MSG_RESULT([none])
    AC_MSG_WARN([no MPI libraries were given])
    echo ""
    echo "   Unable to compile MPI program using C compiler because"
    echo "   MPI libraries were not specified."
    echo ""
    echo "   Try using --with-mpi-libdir and --with-mpi-libs to"
    echo "   specify the location and names of the MPI libraries."
    echo ""
    echo "   Disabling the parallel NVECTOR module and all parallel examples..."
    echo ""
    MPI_C_COMP_OK="no"
    SUNDIALS_WARN_FLAG="yes"
  # MPI libraries were specified so update LIBS
  else
    AC_MSG_RESULT([${MPI_LIBS}])
    if test "X${LIBS}" = "X"; then
      LIBS="${MPI_LIBS}"
    else
      LIBS="${LIBS} ${MPI_LIBS}"
    fi
    # Set the MPI_C_COMP_OK variable to NULL so we can conditionally execute
    # the next test
    MPI_C_COMP_OK=""
  fi

  if test "X${MPI_C_COMP_OK}" = "X"; then
    AC_MSG_CHECKING([if C compiler can compile MPI programs])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include "mpi.h"]],[[int c; char **v; MPI_Init(&c,&v);]])],
    [AC_MSG_RESULT([yes])
     MPI_C_COMP_OK="yes"],
    [AC_MSG_RESULT([no])
     AC_MSG_WARN([C compiler cannot compile MPI programs])
     echo ""
     echo "   Unable to compile MPI program using C compiler."
     echo ""
     echo "   Try using --with-mpicc to specify a MPI-C compiler script,"
     echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
     echo "   to specify the locations of all relevant MPI files, or"
     echo "   --with-mpi-root to specify the base installation directory"
     echo "   of the MPI implementation to be used."
     echo ""
     echo "   Disabling the parallel NVECTOR module and all parallel examples..."
     echo ""
     MPI_C_COMP_OK="no"
     SUNDIALS_WARN_FLAG="yes"])
  fi
else
  MPI_C_COMP_OK="no"
fi
  
# Restore CPPFLAGS, LDFLAGS and LIBS
CPPFLAGS="${SAVED_CPPFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

]) dnl END SUNDIALS_CC_WITH_MPI_CHECK

#------------------------------------------------------------------
# SET MPI-F77 COMPILER
#
# These tests are done only if all of the following are still true:
#  - MPI is enabled
#  - F77 examples are enabled
#  - F77 works
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_MPIF77],
[

AC_MSG_CHECKING([if using MPI-Fortran script])
AC_ARG_WITH(mpif77,
[AC_HELP_STRING([--with-mpif77[[[[=ARG]]]]],[specify MPI-Fortran compiler to use @<:@mpif77@:>@])],
[
if test "X${withval}" = "Xno"; then
  USE_MPIF77_SCRIPT="no"
else
  USE_MPIF77_SCRIPT="yes"
  MPIF77_COMP="${withval}"
fi
],
[
USE_MPIF77_SCRIPT="yes"
MPIF77_COMP="mpif77"
])
AC_MSG_RESULT([${USE_MPIF77_SCRIPT}])

# Check MPI-Fortran compiler (either MPI compiler script or regular Fortran compiler)
if test "X${USE_MPIF77_SCRIPT}" = "Xyes"; then
  SUNDIALS_CHECK_MPIF77
else
  MPIF77_COMP="${F77}"
  MPIF77="${F77}"
  SUNDIALS_F77_WITH_MPI_CHECK
fi

]) dnl END SUNDIALS_SET_MPIF77

#------------------------------------------------------------------
# TEST MPIF77 COMPILER SCRIPT
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_MPIF77],
[

# Test the MPI-Fortran compiler (meaning test MPIF77_COMP)
# Check if MPI-Fortran compiler can be found

AC_MSG_CHECKING([if absolute path to ${MPIF77_COMP} was given])

# CASE 1: MPIF77_COMP was found (cannot check if executable because the
# "-x" flag is NOT portable)
if test -f ${MPIF77_COMP} ; then

  AC_MSG_RESULT([yes])
  MPIF77_COMP_EXISTS="yes"
  # Determine MPI_INC_DIR and MPI_LIB_DIR for use by Makefile
  MPI_BASE_DIR=`AS_DIRNAME(["${MPIF77_COMP}"])`

# CASE 2: MPIF77_COMP could NOT be found and MPI_ROOT_DIR was NOT specified,
# so search in PATH
else

  AC_MSG_RESULT([no])

  if test "X${MPI_ROOT_DIR}" = "X"; then

    # Try to find location of executable (perhaps directory was entered incorrectly)
    TEMP_MPIF77_COMP=`basename "${MPIF77_COMP}"`
    AC_PATH_PROG([MPIF77_COMP],[${TEMP_MPIF77_COMP}],[none])
    # Cannot find executable in PATH
    if test "X${MPIF77_COMP}" = "Xnone"; then
      MPIF77_COMP_EXISTS="no"
      MPIF77_COMP=""
    # Found executable and set MPIF77_COMP to absolute pathname
    else
      MPIF77_COMP_EXISTS="yes"
      MPI_BASE_DIR=`AS_DIRNAME(["${MPIF77_COMP}"])`
    fi

  # CASE 3: MPIF77_COMP could NOT be found, but MPI_ROOT_DIR was specified
  else

    AC_MSG_CHECKING([if ${MPIF77_COMP} exists in ${MPI_ROOT_DIR}/bin])
    # MPIF77_COMP should really only contain an executable name
    # Found location of MPIF77_COMP
    if test -f ${MPI_ROOT_DIR}/bin/${MPIF77_COMP} ; then
      AC_MSG_RESULT([yes])
      MPIF77_COMP_EXISTS="yes"
      MPIF77_COMP="${MPI_ROOT_DIR}/bin/${MPIF77_COMP}"
    # Could NOT find MPIF77_COMP anywhere
    else
      AC_MSG_RESULT([no])
      MPIF77_COMP_EXISTS="no"
      MPIF77_COMP=""
    fi

  fi

fi

# Issue warning message if MPIF77_COMP does NOT exist, else set MPIF77
if test "X${MPIF77_COMP_EXISTS}" = "Xyes"; then

  MPIF77="${MPIF77_COMP}"
  MPI_F77_COMP_OK="yes"

  # Note that we do not have to worry about empty MPI_INC_DIR and MPI_LIB_DIR
  # here as they were set in SUNDIALS_CHECK_MPICC

  # Check if we must use the MPI-Fortran compiler script (MPIF77) to link
  # the Fortran examples (default is to use MPICC)
  SUNDIALS_MPIF77_LNKR_CHECK

else

  AC_MSG_WARN([cannot find MPI-Fortran compiler])
  echo ""
  echo "   Unable to find a functional MPI-Fortran compiler."
  echo ""
  echo "   Try using --with-mpif77 to specify a MPI-Fortran compiler script,"
  echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
  echo "   to specify the locations of all relevant MPI files, or"
  echo "   --with-mpi-root to specify the base installation directory"
  echo "   of the MPI implementation to be used."
  echo ""
  echo "   Disabling parallel Fortran examples...."
  echo ""
  MPIF77=""
  MPI_F77_COMP_OK="no"
  SUNDIALS_WARN_FLAG="yes"

fi

]) dnl END SUNDIALS_CHECK_MPIF77

#------------------------------------------------------------------
# DETERMINE MPI-FORTRAN LINKER IF USING MPIF77 SCRIPT
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_MPIF77_LNKR_CHECK],
[

# If we are NOT using an MPI script, then MPICC_COMP == CC and we do NOT need
# to check again if CC is a C++ compiler as we already know the answer
if test "X${USE_MPICC_SCRIPT}" = "Xyes"; then

  # Check if using a C++ compiler (meaning MPI-C++ script)
  # Save result from CC check
  SAVED_USING_CPLUSPLUS_COMP="${USING_CPLUSPLUS_COMP}"
  SUNDIALS_CPLUSPLUS_CHECK([${MPICC_COMP}])
  # MPICC uses a C++ compiler so run the next test
  if test "X${USING_CPLUSPLUS_COMP}" = "Xyes" && test "X${SAVED_USING_CPLUSPLUS_COMP}" = "Xyes"; then
    RUN_MPIF77_LNKR_CHECK="yes"
  # ERROR
  elif test "X${USING_CPLUSPLUS_COMP}" = "Xyes" && test "X${SAVED_USING_CPLUSPLUS_COMP}" = "Xno"; then
    AC_MSG_ERROR([${MPICC_COMP} is a C++ compiler but ${CC} is a C compiler])
  # MPICC uses a C compiler so skip the next test
  elif test "X${USING_CPLUSPLUS_COMP}" = "Xno" && test "X${SAVED_USING_CPLUSPLUS_COMP}" = "Xno" ; then
    RUN_MPIF77_LNKR_CHECK="no"
  # ERROR
  elif test "X${USING_CPLUSPLUS_COMP}" = "Xno" && test "X${SAVED_USING_CPLUSPLUS_COMP}" = "Xyes" ; then
    AC_MSG_ERROR([${MPICC_COMP} is a C compiler but ${CC} is a C++ compiler])
  fi
  # Restore result from CC check
  USING_CPLUSPLUS_COMP="${SAVED_USING_CPLUSPLUS_COMP}"

else

  AC_MSG_CHECKING([if ${MPICC_COMP} is a C++ compiler])
  if test "X${USING_CPLUSPLUS_COMP}" = "Xyes"; then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
  fi

fi

AC_MSG_CHECKING([which linker to use])
# Perform the next test only if using a C++ compiler to build NVECTOR_PARALLEL
if test "X${RUN_MPIF77_LNKR_CHECK}" = "Xyes"; then

  MPIF77_LNKR_CHECK_OK="no"

  # Switch language to "Fortran 77"
  AC_LANG_PUSH([Fortran 77])

  # Temporarily reset F77 environment variable to perform test
  SAVED_F77="${F77}"
  F77="${MPIF77_COMP}"

  # Compile simple Fortran example, but do NOT link
  # Note: result stored as conftest.${ac_objext}
  AC_COMPILE_IFELSE(
  [AC_LANG_SOURCE(
  [[
	PROGRAM SUNDIALS
	INTEGER IER
	CALL MPI_INIT(IER)
	END
  ]])],
  [

  # Reset F77 to original value
  F77="${SAVED_F77}"

  # Revert to previous language
  AC_LANG_POP([Fortran 77])

  # Temporarily reset LIBS environment variable to perform test
  SAVED_LIBS="${LIBS}"
  LIBS="${LIBS} ${FLIBS}"

  # Switch working language to C for next test
  AC_LANG_PUSH([C])

  # Temporarily reset CC environment variable to perform next test
  SAVED_CC="${CC}"
  CC="${MPICC_COMP}"

  # Check if MPICC_COMP can link Fortran example
  # Note: AC_LINKONLY_IFELSE is a custom macro (modifications made to
  # general.m4 and c.m4)
  AC_LINKONLY_IFELSE([],[MPIF77_LNKR_CHECK_OK="yes"],[MPIF77_LNKR_CHECK_OK="no"])

  # Reset CC to original value
  CC="${SAVED_CC}"

  # Revert back to previous language (Fortran 77)
  AC_LANG_POP([C])

  # Set LIBS environment variable back to original value
  LIBS="${SAVED_LIBS}"

  ])

  # If either the compilation or the linking failed, we should
  # disable building the parallel Fortran examples
  # For now, use MPIF77 as the linker...

  if test "X${MPIF77_LNKR_CHECK_OK}" = "Xyes"; then
    MPIF77_LNKR="${MPICC}"
  else
    MPIF77_LNKR="${MPIF77}"
  fi

else

  # Using a C compiler so use MPIF77 to link parallel Fortran examples
  MPIF77_LNKR="${MPIF77}"

fi
AC_MSG_RESULT([${MPIF77_LNKR}])


]) dnl SUNDIALS_MPIF77_LNKR_CHECK

#------------------------------------------------------------------
# TEST FORTRAN COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_WITH_MPI_CHECK],
[

# Test if we can compile MPI programs using the F77 compiler
# and current MPI settings

AC_MSG_NOTICE([Testing F77 with MPI settings])

AC_LANG_PUSH([Fortran 77])

# Save copies of FFLAGS, LDFLAGS and LIBS (preserve information)
# Temporarily overwritten so we can test MPI implementation
SAVED_FFLAGS="${FFLAGS}"
SAVED_LDFLAGS="${LDFLAGS}"
SAVED_LIBS="${LIBS}"

# This may seem redundant, but we are not guaranteed that
# SUNDIALS_CC_WITH_MPI_CHECK has been executed
# Determine location of MPI header files (find MPI include directory)
MPI_EXISTS="yes"

AC_MSG_CHECKING([for location of MPI implementation])

# If MPI include directory was NOT explicitly specified so check if MPI root
# directory was given by user
if test "X${MPI_INC_DIR}" = "X"; then
  # If MPI root directory was NOT given so issue a warning message
  if test "X${MPI_ROOT_DIR}" = "X"; then
    AC_MSG_RESULT([not found])
    MPI_EXISTS="no"
    AC_MSG_WARN([cannot find MPI implementation files])
    echo ""
    echo "   Unable to find MPI implementation files."
    echo ""
    echo "   Try using --with-mpif77 to specify a MPI-Fortran compiler script,"
    echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
    echo "   to specify the locations of all relevant MPI files, or"
    echo "   --with-mpi-root to specify the base installation directory"
    echo "   of the MPI implementation to be used."
    echo ""
    echo "   Disabling all parallel Fortran examples..."
    echo ""
    SUNDIALS_WARN_FLAG="yes"
  # MPI root directory was given so set MPI_INC_DIR accordingly
  # Update FFLAGS
  else
    MPI_INC_DIR="${MPI_ROOT_DIR}/include"
    AC_MSG_RESULT([${MPI_INC_DIR}])
    if test "X${FFLAGS}" = "X"; then
      FFLAGS="-I${MPI_INC_DIR}"
    else
      FFLAGS="${FFLAGS} -I${MPI_INC_DIR}"
    fi
  fi
# MPI include directory was specified so update FFLAGS
else
  AC_MSG_RESULT([${MPI_INC_DIR}])
  if test "X${FFLAGS}" = "X"; then
    FFLAGS="-I${MPI_INC_DIR}"
  else
    FFLAGS="${FFLAGS} -I${MPI_INC_DIR}"
  fi
fi

# Only continue if found an MPI implementation
if test "X${MPI_EXISTS}" = "Xyes"; then

  AC_MSG_CHECKING([for location of MPI libraries])

  # Determine location of MPI libraries
  # MPI library directory was NOT specified by user so set based upon MPI_ROOT_DIR
  # Update LDFLAGS
  if test "X${MPI_LIB_DIR}" = "X"; then
    MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    AC_MSG_RESULT([${MPI_LIB_DIR}])
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  # MPI library directory was specified so update LDFLAGS
  else
    AC_MSG_RESULT([${MPI_LIB_DIR}])
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  fi

  # Check if user specified which MPI libraries must be included
  # If no libraries are given, then issue a warning message
  AC_MSG_CHECKING([for MPI libraries])
  if test "X${MPI_LIBS}" = "X"; then
    AC_MSG_RESULT([none])
    echo ""
    echo "   Unable to compile MPI program using Fortran compiler because"
    echo "   MPI libraries were not specified."
    echo ""
    echo "   Try using --with-mpi-libdir and --with-mpi-libs to"
    echo "   specify the location and names of the MPI libraries."
    echo ""
    echo "   Disabling all parallel Fortran examples..."
    echo ""
    MPI_F77_COMP_OK="no"
  # MPI libraries were specified so update LIBS
  else
    AC_MSG_RESULT([${MPI_LIBS}])
    if test "X${LIBS}" = "X"; then
      LIBS="${MPI_LIBS}"
    else
      LIBS="${LIBS} ${MPI_LIBS}"
    fi
    # Set the MPI_F77_COMP_OK variable to NULL so we can conditionally execute
    # the next test
    MPI_F77_COMP_OK=""
  fi

  if test "X${MPI_F77_COMP_OK}" = "X"; then
    AC_MSG_CHECKING([if Fortran compiler can compile MPI programs])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([],
    [      
        INCLUDE "mpif.h"
        CALL MPI_INIT(IER)
    ])],
    [AC_MSG_RESULT([yes])
     MPI_F77_COMP_OK="yes"],
    [AC_MSG_RESULT([no])
     AC_MSG_WARN([Fortran compiler cannot compile MPI programs])
     echo ""
     echo "   Unable to compile MPI program using Fortran compiler."
     echo ""
     echo "   Try using --with-mpif77 to specify a MPI-Fortran compiler script,"
     echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
     echo "   to specify the locations of all relevant MPI files, or"
     echo "   --with-mpi-root to specify the base installation directory"
     echo "   of the MPI implementation to be used."
     echo ""
     echo "   Disabling all parallel Fortran examples..."
     echo ""
     MPI_F77_COMP_OK="no"
     SUNDIALS_WARN_FLAG="yes"])

    # Set MPIF77_LNKR based on value of F77_LNKR
    # Note: setting MPIF77_LNKR is trivial if NOT using the MPI compiler script
    # since the SUNDIALS_F77_LNKR_CHECK macro already checked if CC or F77
    # should be used
    AC_MSG_CHECKING([which linker to use])
    if test "X${F77_LNKR}" = "X${CC}"; then
      MPIF77_LNKR="${MPICC}"
    elif test "X${F77_LNKR}" = "X${F77}"; then
      MPIF77_LNKR="${MPIF77}"
    fi
    AC_MSG_RESULT([${MPIF77_LNKR}])
  fi

else
  MPI_F77_COMP_OK="no"
fi

# Restore FFLAGS, LDFLAGS and LIBS
FFLAGS="${SAVED_FFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

AC_LANG_POP([Fortran 77])

]) dnl END SUNDIALS_F77_WITH_MPI_CHECK

#------------------------------------------------------------------
# TEST MPI-2 FUNCTIONALITY
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_MPI2],
[

# Determine if MPI implementation used to build SUNDIALS provides
# MPI-2 functionality.
#
# Test for MPI_Comm_f2c() function:
#   (1) NO  : FNVECTOR_PARALLEL module will NOT allow user to specify
#             an MPI communicator and MPI_COMM_WORLD will be used
#   (2) YES : FNVECTOR_PARALLEL module will allow user to specify
#             an MPI communicator
#
 
# Provide variable description templates for config.hin and config.h files
# Required by autoheader utility
AH_TEMPLATE([SUNDIALS_MPI_COMM_F2C],
            [FNVECTOR: Allow user to specify different MPI communicator])

# Save copies of CPPFLAGS, LDFLAGS and LIBS (preserve information)
# Temporarily overwritten so we can test MPI implementation
SAVED_CPPFLAGS="${CPPFLAGS}"
SAVED_LDFLAGS="${LDFLAGS}"
SAVED_LIBS="${LIBS}"

# Determine location of MPI header files (find MPI include directory)
MPI_EXISTS="yes"

# MPI include directory was NOT explicitly specified so check if MPI root
# directory was given by user
if test "X${MPI_INC_DIR}" = "X"; then
  # MPI root directory was NOT given so issue a warning message
  if test "X${MPI_ROOT_DIR}" = "X"; then
    MPI_EXISTS="no"
    AC_MSG_WARN([cannot find MPI implementation files])
    echo ""
    echo "   Unable to find MPI implementation files."
    echo ""
    echo "   Try using --with-mpicc to specify a MPI-C compiler script,"
    echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
    echo "   to specify the locations of all relevant MPI files, or"
    echo "   --with-mpi-root to specify the base installation directory"
    echo "   of the MPI implementation to be used."
    echo ""
    echo "   Disabling FNVECTOR_PARALLEL support for user-specified"
    echo "   MPI communicator..."
    echo ""
    SUNDIALS_WARN_FLAG="yes"
  # MPI root directory was given so set MPI_INC_DIR accordingly
  # Update CPPFLAGS
  else
    MPI_INC_DIR="${MPI_ROOT_DIR}/include"
    if test "X${CPPFLAGS}" = "X"; then
      CPPFLAGS="-I${MPI_INC_DIR}"
    else
      CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
    fi
    # Add MPI_FLAGS if non-empty
    if test "X${MPI_FLAGS}" = "X"; then
      CPPFLAGS="${CPPFLAGS}"
    else
      CPPFLAGS="${CPPFLAGS} ${MPI_FLAGS}"
    fi
  fi
# MPI include directory was specified so update CPPFLAGS
else
  if test "X${CPPFLAGS}" = "X"; then
    CPPFLAGS="-I${MPI_INC_DIR}"
  else
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
  fi
  # Add MPI_FLAGS if non-empty
  if test "X${MPI_FLAGS}" = "X"; then
    CPPFLAGS="${CPPFLAGS}"
  else
    CPPFLAGS="${CPPFLAGS} ${MPI_FLAGS}"
  fi
fi

# Only continue if found an MPI implementation
if test "X${MPI_EXISTS}" = "Xyes"; then

  # Determine location of MPI libraries
  # MPI library directory was NOT specified by user so set based upon MPI_ROOT_DIR
  # Update LDFLAGS
  if test "X${MPI_LIB_DIR}" = "X"; then
    MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  # MPI library directory was specified so update LDFLAGS
  else
    if test "X${LDFLAGS}" = "X"; then
      LDFLAGS="-L${MPI_LIB_DIR}"
    else
      LDFLAGS="${LDFLAGS} -L${MPI_LIB_DIR}"
    fi
  fi

  # Check if user specified which MPI libraries linker should be use
  if test "X${MPI_LIBS}" = "X"; then
    :
  # MPI libraries were specified so update LIBS
  else
    if test "X${LIBS}" = "X"; then
      LIBS="${MPI_LIBS}"
    else
      LIBS="${LIBS} ${MPI_LIBS}"
    fi
  fi

  # Since AC_LINK_IFELSE uses CC, set CC = MPICC if using
  # an MPI compiler script
  if test "X${USE_MPICC_SCRIPT}" = "Xyes"; then
    SAVED_CC="${CC}"
    CC="${MPICC_COMP}"
  fi

  # Check if MPI implementation supports MPI_Comm_f2c() from
  # MPI-2 specification
  if test "X${FCMIX_ENABLED}" = "Xyes"; then
    AC_MSG_CHECKING([for MPI_Comm_f2c() from MPI-2 specification])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include "mpi.h"]],
    [[
        int c;
        char **v;
        MPI_Comm C_comm;
        MPI_Init(&c, &v);
        C_comm = MPI_Comm_f2c((MPI_Fint) 1);
        MPI_Finalize();
    ]])],
    [AC_MSG_RESULT([yes])
     AC_DEFINE([SUNDIALS_MPI_COMM_F2C],[1],[])
     F77_MPI_COMM_F2C="#define SUNDIALS_MPI_COMM_F2C 1"],
    [AC_MSG_RESULT([no])
     AC_DEFINE([SUNDIALS_MPI_COMM_F2C],[0],[])
     F77_MPI_COMM_F2C="#define SUNDIALS_MPI_COMM_F2C 0"])
  fi

  # Reset CC if necessary
  if test "X${USE_MPICC_SCRIPT}" = "Xyes"; then
    CC="${SAVED_CC}"
  fi

else
  AC_DEFINE([SUNDIALS_MPI_COMM_F2C],[0],[])
  F77_MPI_COMM_F2C="#define SUNDIALS_MPI_COMM_F2C 0"
fi

# Restore CPPFLAGS, LDFLAGS and LIBS
CPPFLAGS="${SAVED_CPPFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

]) dnl END SUNDIALS_CHECK_MPI2




#=================================================================#
#                                                                 #
#                                                                 #
#                   F I N A L I Z A T I O N S                     #
#                                                                 #
#                                                                 #
#==================================================================


#------------------------------------------------------------------
# ADD SOME MORE STUFF TO configure --help
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_MORE_HELP],
[

AC_ARG_WITH([],[          ],[])
AC_ARG_WITH([],[NOTES],[])
AC_ARG_WITH([],[  It is legal to set --with-exinstdir to "no", in which case the examples],[])
AC_ARG_WITH([],[  are built but not installed.],[])
AC_ARG_WITH([],[  Enabling the compilation of the examples (--enable-examples) but disabling their],[])
AC_ARG_WITH([],[  installation (--with-exinstdir=no) can be used to test the SUNDIALS libraries.],[])

]) dnl END SUNDIALS_MORE_HELP

#------------------------------------------------------------------
# SET EXAMPLES
#
# Decide which examples can be built
#
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_EXAMPLES],
[

# Set proper object file extension
# Must export OBJ_EXT via AC_SUBST
OBJEXT=".lo"

# Check if serial C examples can actually be built
SERIAL_C_EXAMPLES="yes"

# Check if parallel C examples can actually be built
if test "X${MPI_ENABLED}" = "Xyes"; then
  if test "X${MPI_C_COMP_OK}" = "Xyes"; then
    PARALLEL_C_EXAMPLES="yes"
  else
    PARALLEL_C_EXAMPLES="no"
  fi
else
  PARALLEL_C_EXAMPLES="disabled"
fi

# Check if serial F77 examples can actually be built
if test "X${FCMIX_ENABLED}" = "Xyes"; then
  if test "X${F77_OK}" = "Xyes"; then
    SERIAL_F77_EXAMPLES="yes"
  else
    SERIAL_F77_EXAMPLES="no"
  fi
else
  SERIAL_F77_EXAMPLES="disabled"
fi

# Check if parallel F77 examples can actually be built
if test "X${FCMIX_ENABLED}" = "Xyes" &&  test "X${MPI_ENABLED}" = "Xyes"; then
  if test "X${MPI_F77_COMP_OK}" = "Xyes"; then
    PARALLEL_F77_EXAMPLES="yes"
  else
    PARALLEL_F77_EXAMPLES="no"
  fi
else
  PARALLEL_F77_EXAMPLES="disabled"
fi

# Notify user
AC_MSG_CHECKING([if we can build serial C examples])
AC_MSG_RESULT([${SERIAL_C_EXAMPLES}])
AC_MSG_CHECKING([if we can build parallel C examples])
AC_MSG_RESULT([${PARALLEL_C_EXAMPLES}])
AC_MSG_CHECKING([if we can build serial Fortran examples])
AC_MSG_RESULT([${SERIAL_F77_EXAMPLES}])
AC_MSG_CHECKING([if we can build parallel Fortran examples])
AC_MSG_RESULT([${PARALLEL_F77_EXAMPLES}])

# Check if the Fortran update script (bin/fortran-update.in) is needed
if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" || test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
  BUILD_F77_UPDATE_SCRIPT="yes";
else
  BUILD_F77_UPDATE_SCRIPT="no"
fi

# Where should we install the examples?
# Note: setting this to "no" will disable example installation!
AC_MSG_CHECKING([where to install the SUNDIALS examples])
AC_ARG_WITH([],[   ],[])
AC_ARG_WITH([exinstdir],
[AC_HELP_STRING([--with-exinstdir=DIR], [install SUNDIALS examples in DIR @<:@EPREFIX/examples@:>@])],
[
  EXS_INSTDIR="${withval}"
],
[
  if test "X${exec_prefix}" = "XNONE"; then
    if test "X${prefix}" = "XNONE"; then
      EXS_INSTDIR="\${exec_prefix}/examples"
    else
      EXS_INSTDIR="${prefix}/examples"
    fi
  else
    EXS_INSTDIR="${exec_prefix}/examples"
  fi
])
AC_MSG_RESULT([${EXS_INSTDIR}])

# Prepare substitution variables to create the exported example Makefiles

F77_LIBS="${FLIBS} ${LIBS}"
if test "X${F77_LNKR}" = "X${F77}"; then
  F77_LDFLAGS="${FFLAGS} ${LDFLAGS}"
else
  F77_LDFLAGS="${CFLAGS} ${LDFLAGS}"
fi

]) dnl END SUNDIALS_SET_EXAMPLES

#------------------------------------------------------------------
# BUILD MODULES LIST
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_BUILD_MODULES_LIST],
[

# Initialize the list of Makefiles to be created
SUNDIALS_MAKEFILES="Makefile"

# Initialize list of additional configure files to be created
SUNDIALS_CONFIGFILES="include/sundials/sundials_config.h:include/sundials/sundials_config.in"
SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} bin/sundials-config:bin/sundials-config.in"

# Initialize lists of solver modules and example modules
SLV_MODULES="src/sundials"
SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/sundials/Makefile"

EXS_MODULES=""

# NVECTOR modules
if test -d ${srcdir}/src/nvec_ser ; then
  SLV_MODULES="${SLV_MODULES} src/nvec_ser"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/nvec_ser/Makefile"
fi

if test -d ${srcdir}/src/nvec_par && test "X${MPI_C_COMP_OK}" = "Xyes"; then
  SLV_MODULES="${SLV_MODULES} src/nvec_par"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/nvec_par/Makefile"
fi

if test -d ${srcdir}/src/nvec_spcpar && test "X${MPI_C_COMP_OK}" = "Xyes"; then
  SLV_MODULES="${SLV_MODULES} src/nvec_spcpar"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/nvec_spcpar/Makefile"
fi

# CVODE module
if test "X${CVODE_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/cvode"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/cvode/Makefile"

  if test "X${FCMIX_ENABLED}" = "Xyes" && test -d ${srcdir}/src/cvode/fcmix ; then
    SLV_MODULES="${SLV_MODULES} src/cvode/fcmix"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/cvode/fcmix/Makefile"
  fi

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/serial/Makefile_ex:examples/templates/makefile_serial_C_ex.in"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/fcmix_serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/fcmix_serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/fcmix_serial/Makefile_ex:examples/templates/makefile_serial_F77_ex.in"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/parallel/Makefile_ex:examples/templates/makefile_parallel_C_ex.in"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/fcmix_parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/fcmix_parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/fcmix_parallel/Makefile_ex:examples/templates/makefile_parallel_F77_ex.in"
  fi

fi

# CVODES module
if test "X${CVODES_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/cvodes"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/cvodes/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvodes/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvodes/serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvodes/serial/Makefile_ex:examples/templates/makefile_serial_C_ex.in"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvodes/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvodes/parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvodes/parallel/Makefile_ex:examples/templates/makefile_parallel_C_ex.in"
  fi

fi

# IDA module
if test "X${IDA_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/ida"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/ida/Makefile"

  if test "X${FCMIX_ENABLED}" = "Xyes" && test -d ${srcdir}/src/ida/fcmix ; then
    SLV_MODULES="${SLV_MODULES} src/ida/fcmix"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/ida/fcmix/Makefile"
  fi

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/serial/Makefile_ex:examples/templates/makefile_serial_C_ex.in"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/fcmix_serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/fcmix_serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/fcmix_serial/Makefile_ex:examples/templates/makefile_serial_F77_ex.in"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/parallel/Makefile_ex:examples/templates/makefile_parallel_C_ex.in"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/fcmix_parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/fcmix_parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/fcmix_parallel/Makefile_ex:examples/templates/makefile_parallel_F77_ex.in"
  fi

fi

# IDAS module
if test "X${IDAS_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/idas"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/idas/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/idas/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/idas/serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/idas/serial/Makefile_ex:examples/templates/makefile_serial_C_ex.in"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/idas/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/idas/parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/idas/parallel/Makefile_ex:examples/templates/makefile_parallel_C_ex.in"
  fi

fi

# KINSOL module
if test "X${KINSOL_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/kinsol"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/kinsol/Makefile"

  if test "X${FCMIX_ENABLED}" = "Xyes" && test -d ${srcdir}/src/kinsol/fcmix ; then
    SLV_MODULES="${SLV_MODULES} src/kinsol/fcmix"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/kinsol/fcmix/Makefile"
  fi

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/serial/Makefile_ex:examples/templates/makefile_serial_C_ex.in"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/fcmix_serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/fcmix_serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/fcmix_serial/Makefile_ex:examples/templates/makefile_serial_F77_ex.in"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/parallel/Makefile_ex:examples/templates/makefile_parallel_C_ex.in"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/fcmix_parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/fcmix_parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/fcmix_parallel/Makefile_ex:examples/templates/makefile_parallel_F77_ex.in"
  fi

fi

# CPODES module
if test "X${CPODES_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/cpodes"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/cpodes/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cpodes/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cpodes/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cpodes/serial/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cpodes/serial/Makefile_ex:examples/templates/makefile_serial_C_ex.in"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cpodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cpodes/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cpodes/parallel/Makefile"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cpodes/parallel/Makefile_ex:examples/templates/makefile_parallel_C_ex.in"
  fi

fi

# Add Fortran update script to the list of additional files to be generated
if test "X${BUILD_F77_UPDATE_SCRIPT}" = "Xyes"; then
  SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} bin/fortran-update.sh:bin/fortran-update.in"
fi

# If needed, add Makefile update script to the list of additional files to be generated
if test "X${EXAMPLES_ENABLED}" = "Xyes" && test "X${EXS_INSTDIR}" != "Xno"; then
  SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} bin/makefile-update.sh:bin/makefile-update.in"
fi

]) dnl END SUNDIALS_BUILD_MODULES_LIST

#------------------------------------------------------------------
# POST PROCESSING OF EXAMPLE Makefiles for export
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_POST_PROCESSING],
[

# If installing examples, the Makefiles that will be exported must
# be post-processed to complete the substitution of all variables. 
# After config.status runs, each example subdirectory contains an 
# export makefile, named Makefile_ex, which was created from the 
# common template in examples/templates.
#
# The following variables are still to be substituted at this point:
#   SOLVER
#   EXAMPLES
#   EXAMPLES_BL
#   SOLVER_LIB SOLVER_FLIB
#   NVEC_LIB NVEC_FLIB
# 
# This function is called ONLY if examples are enabled AND examples will 
# be installed. If so, it sets up commands to be called after config.status
# has generated a first version of the Makefiles for export:
# 
# (1) For each solver, proceed ONLY if the solver is enabled.
# (2) For each type of examples, proceed ONLY if they can be compiled AND
#     the example directory exists.

# CVODE module
if test "X${CVODE_ENABLED}" = "Xyes"; then

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([cvode_ser_ex_bl],
       [   
       IN_FILE="examples/cvode/serial/Makefile_ex"
       SOLVER="CVODE"
       SOLVER_LIB="sundials_cvode"
       SOLVER_FLIB=""
       EXAMPLES="cvAdvDiff_bnd cvDirectDemo_ls cvDiurnal_kry_bp cvDiurnal_kry cvKrylovDemo_ls cvKrylovDemo_prec cvRoberts_dns cvRoberts_dns_uw"
       EXAMPLES_BL="cvAdvDiff_bndL cvRoberts_dnsL"
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([cvode_ser_ex],
       [   
       IN_FILE="examples/cvode/serial/Makefile_ex"
       SOLVER="CVODE"
       SOLVER_LIB="sundials_cvode"
       SOLVER_FLIB=""
       EXAMPLES="cvAdvDiff_bnd cvDirectDemo_ls cvDiurnal_kry_bp cvDiurnal_kry cvKrylovDemo_ls cvKrylovDemo_prec cvRoberts_dns cvRoberts_dns_uw"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([cvode_fser_ex_bl],
       [
       IN_FILE="examples/cvode/fcmix_serial/Makefile_ex"
       SOLVER="CVODE"
       SOLVER_LIB="sundials_cvode"
       SOLVER_FLIB="sundials_fcvode"
       EXAMPLES="fcvAdvDiff_bnd fcvDiurnal_kry_bp fcvDiurnal_kry fcvRoberts_dns"
       EXAMPLES_BL="fcvRoberts_dnsL"
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([cvode_fser_ex],
       [
       IN_FILE="examples/cvode/fcmix_serial/Makefile_ex"
       SOLVER="CVODE"
       SOLVER_LIB="sundials_cvode"
       SOLVER_FLIB="sundials_fcvode"
       EXAMPLES="fcvAdvDiff_bnd fcvDiurnal_kry_bp fcvDiurnal_kry fcvRoberts_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/parallel ; then
     AC_CONFIG_COMMANDS([cvode_par_ex],
     [
     IN_FILE="examples/cvode/parallel/Makefile_ex"
     SOLVER="CVODE"
     SOLVER_LIB="sundials_cvode"
     SOLVER_FLIB=""
     EXAMPLES="cvAdvDiff_non_p cvDiurnal_kry_bbd_p cvDiurnal_kry_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_parallel ; then
     AC_CONFIG_COMMANDS([cvode_fpar_ex],
     [
     IN_FILE="examples/cvode/fcmix_parallel/Makefile_ex"
     SOLVER="CVODE"
     SOLVER_LIB="sundials_cvode"
     SOLVER_FLIB="sundials_fcvode"
     EXAMPLES="fcvDiag_non_p fcvDiag_kry_bbd_p fcvDiag_kry_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

fi


# CVODES module
if test "X${CVODES_ENABLED}" = "Xyes"; then

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([cvodes_ser_ex_bl],
       [
       IN_FILE="examples/cvodes/serial/Makefile_ex"
       SOLVER="CVODES"
       SOLVER_LIB="sundials_cvodes"
       SOLVER_FLIB=""
       EXAMPLES="cvsAdvDiff_ASAi_bnd cvsAdvDiff_FSA_non cvsDiurnal_kry_bp cvsFoodWeb_ASAp_kry cvsKrylovDemo_prec cvsAdvDiff_bnd cvsDirectDemo_ls cvsDiurnal_kry cvsHessian_ASA_FSA cvsRoberts_ASAi_dns cvsRoberts_dns_uw  cvsDiurnal_FSA_kry cvsFoodWeb_ASAi_kry cvsKrylovDemo_ls cvsRoberts_dns cvsRoberts_FSA_dns"
       EXAMPLES_BL="cvsRoberts_dnsL cvsAdvDiff_bndL"
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([cvodes_ser_ex],
       [
       IN_FILE="examples/cvodes/serial/Makefile_ex"
       SOLVER="CVODES"
       SOLVER_LIB="sundials_cvodes"
       SOLVER_FLIB=""
       EXAMPLES="cvsAdvDiff_ASAi_bnd cvsAdvDiff_FSA_non cvsDiurnal_kry_bp cvsFoodWeb_ASAp_kry cvsKrylovDemo_prec cvsAdvDiff_bnd cvsDirectDemo_ls cvsDiurnal_kry cvsHessian_ASA_FSA cvsRoberts_ASAi_dns cvsRoberts_dns_uw  cvsDiurnal_FSA_kry cvsFoodWeb_ASAi_kry cvsKrylovDemo_ls cvsRoberts_dns cvsRoberts_FSA_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/parallel ; then
     AC_CONFIG_COMMANDS([cvodes_par_ex],
     [
     IN_FILE="examples/cvodes/parallel/Makefile_ex"
     SOLVER="CVODES"
     SOLVER_LIB="sundials_cvodes"
     SOLVER_FLIB=""
     EXAMPLES="cvsAdvDiff_ASAp_non_p cvsAdvDiff_non_p cvsDiurnal_FSA_kry_p cvsDiurnal_kry_p cvsAdvDiff_FSA_non_p cvsAtmDisp_ASAi_kry_bbd_p cvsDiurnal_kry_bbd_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

fi


# IDA module
if test "X${IDA_ENABLED}" = "Xyes"; then

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([ida_ser_ex_bl],
       [
       IN_FILE="examples/ida/serial/Makefile_ex"
       SOLVER="IDA"
       SOLVER_LIB="sundials_ida"
       SOLVER_FLIB=""
       EXAMPLES="idaFoodWeb_bnd idaHeat2D_bnd idaHeat2D_kry idaKrylovDemo_ls idaRoberts_dns idaSlCrank_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([ida_ser_ex],
       [
       IN_FILE="examples/ida/serial/Makefile_ex"
       SOLVER="IDA"
       SOLVER_LIB="sundials_ida"
       SOLVER_FLIB=""
       EXAMPLES="idaFoodWeb_bnd idaHeat2D_bnd idaHeat2D_kry idaKrylovDemo_ls idaRoberts_dns idaSlCrank_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([ida_fser_ex_bl],
       [
       IN_FILE="examples/ida/fcmix_serial/Makefile_ex"
       SOLVER="IDA"
       SOLVER_LIB="sundials_ida"
       SOLVER_FLIB="sundials_fida"
       EXAMPLES="fidaRoberts_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([ida_fser_ex],
       [
       IN_FILE="examples/ida/fcmix_serial/Makefile_ex"
       SOLVER="IDA"
       SOLVER_LIB="sundials_ida"
       SOLVER_FLIB="sundials_fida"
       EXAMPLES="fidaRoberts_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/parallel ; then
     AC_CONFIG_COMMANDS([ida_par_ex],
     [
     IN_FILE="examples/ida/parallel/Makefile_ex"
     SOLVER="IDA"
     SOLVER_LIB="sundials_ida"
     SOLVER_FLIB=""
     EXAMPLES="idaFoodWeb_kry_bbd_p idaFoodWeb_kry_p idaHeat2D_kry_bbd_p idaHeat2D_kry_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_parallel ; then
     AC_CONFIG_COMMANDS([ida_fpar_ex],
     [
     IN_FILE="examples/ida/fcmix_parallel/Makefile_ex"
     SOLVER="IDA"
     SOLVER_LIB="sundials_ida"
     SOLVER_FLIB="sundials_fida"
     EXAMPLES="fidaHeat2D_kry_bbd_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

fi


# IDAS module
if test "X${IDAS_ENABLED}" = "Xyes"; then

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([idas_ser_ex_bl],
       [
       IN_FILE="examples/idas/serial/Makefile_ex"
       SOLVER="IDAS"
       SOLVER_LIB="sundials_idas"
       SOLVER_FLIB=""
       EXAMPLES="idasAkzoNob_ASAi_dns idasFoodWeb_bnd idasHeat2D_kry idasKrylovDemo_ls idasRoberts_dns idasSlCrank_dns idasAkzoNob_dns idasHeat2D_bnd idasHessian_ASA_FSA idasRoberts_ASAi_dns idasRoberts_FSA_dns idasSlCrank_FSA_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([idas_ser_ex],
       [
       IN_FILE="examples/idas/serial/Makefile_ex"
       SOLVER="IDAS"
       SOLVER_LIB="sundials_idas"
       SOLVER_FLIB=""
       EXAMPLES="idasAkzoNob_ASAi_dns idasFoodWeb_bnd idasHeat2D_kry idasKrylovDemo_ls idasRoberts_dns idasSlCrank_dns idasAkzoNob_dns idasHeat2D_bnd idasHessian_ASA_FSA idasRoberts_ASAi_dns idasRoberts_FSA_dns idasSlCrank_FSA_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/parallel ; then
     AC_CONFIG_COMMANDS([idas_par_ex],
     [
     IN_FILE="examples/idas/parallel/Makefile_ex"
     SOLVER="IDAS"
     SOLVER_LIB="sundials_idas"
     SOLVER_FLIB=""
     EXAMPLES="idasBruss_ASAp_kry_bbd_p idasBruss_kry_bbd_p idasFoodWeb_kry_p idasHeat2D_kry_bbd_p idasBruss_FSA_kry_bbd_p idasFoodWeb_kry_bbd_p idasHeat2D_FSA_kry_bbd_p idasHeat2D_kry_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

fi


# KINSOL module
if test "X${KINSOL_ENABLED}" = "Xyes"; then

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([kinsol_ser_ex_bl],
       [
       IN_FILE="examples/kinsol/serial/Makefile_ex"
       SOLVER="KINSOL"
       SOLVER_LIB="sundials_kinsol"
       SOLVER_FLIB=""
       EXAMPLES="kinFerTron_dns kinFoodWeb_kry kinKrylovDemo_ls kinLaplace_bnd kinRoboKin_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([kinsol_ser_ex],
       [
       IN_FILE="examples/kinsol/serial/Makefile_ex"
       SOLVER="KINSOL"
       SOLVER_LIB="sundials_kinsol"
       SOLVER_FLIB=""
       EXAMPLES="kinFerTron_dns kinFoodWeb_kry kinKrylovDemo_ls kinLaplace_bnd kinRoboKin_dns"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([kinsol_fser_ex_bl],
       [
       IN_FILE="examples/kinsol/fcmix_serial/Makefile_ex"
       SOLVER="KINSOL"
       SOLVER_LIB="sundials_kinsol"
       SOLVER_FLIB="sundials_fkinsol"
       EXAMPLES="fkinDiagon_kry"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([kinsol_fser_ex],
       [
       IN_FILE="examples/kinsol/fcmix_serial/Makefile_ex"
       SOLVER="KINSOL"
       SOLVER_LIB="sundials_kinsol"
       SOLVER_FLIB="sundials_fkinsol"
       EXAMPLES="fkinDiagon_kry"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/parallel ; then
     AC_CONFIG_COMMANDS([kinsol_par_ex],
     [
     IN_FILE="examples/kinsol/parallel/Makefile_ex"
     SOLVER="KINSOL"
     SOLVER_LIB="sundials_kinsol"
     SOLVER_FLIB=""
     EXAMPLES="kinFoodWeb_kry_bbd_p kinFoodWeb_kry_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_parallel ; then
     AC_CONFIG_COMMANDS([kinsol_fpar_ex],
     [
     IN_FILE="examples/kinsol/fcmix_parallel/Makefile_ex"
     SOLVER="KINSOL"
     SOLVER_LIB="sundials_kinsol"
     SOLVER_FLIB="sundials_fkinsol"
     EXAMPLES="fkinDiagon_kry_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

fi


# CPODES module
if test "X${CPODES_ENABLED}" = "Xyes"; then

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cpodes/serial ; then
     if test "X${LAPACK_ENABLED}" = "Xyes"; then
       AC_CONFIG_COMMANDS([cpodes_ser_ex_bl],
       [
       IN_FILE="examples/cpodes/serial/Makefile_ex"
       SOLVER="CPODES"
       SOLVER_LIB="sundials_cpodes"
       SOLVER_FLIB=""
       EXAMPLES="cpsAdvDiff_bnd cpsAdvDiff_non cpsNewtCrd_dns cpsPend_dns cpsRoberts_dns cpsVanDPol_non"
       EXAMPLES_BL="cpsAdvDiff_bndL cpsPend_dnsL cpsRoberts_dnsL"
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     else
       AC_CONFIG_COMMANDS([cpodes_ser_ex],
       [
       IN_FILE="examples/cpodes/serial/Makefile_ex"
       SOLVER="CPODES"
       SOLVER_LIB="sundials_cpodes"
       SOLVER_FLIB=""
       EXAMPLES="cpsAdvDiff_bnd cpsAdvDiff_non cpsNewtCrd_dns cpsPend_dns cpsRoberts_dns cpsVanDPol_non"
       EXAMPLES_BL=""
       ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
       ])
     fi
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cpodes/parallel ; then
     AC_CONFIG_COMMANDS([cpodes_par_ex],
     [
     IN_FILE="examples/cpodes/parallel/Makefile_ex"
     SOLVER="CPODES"
     SOLVER_LIB="sundials_cpodes"
     SOLVER_FLIB=""
     EXAMPLES="cpsHeat2D_kry_bbd_p"
     EXAMPLES_BL=""
     ${SHELL} bin/makefile-update.sh "${IN_FILE}" "${SOLVER}" "${EXAMPLES}" "${EXAMPLES_BL}" "${SOLVER_LIB}" "${SOLVER_FLIB}"
     ])
  fi

fi

]) dnl END SUNDIALS_POST_PROCESSING

#------------------------------------------------------------------
# PRINT STATUS REPORT
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_REPORT],
[

if test "X${SUNDIALS_WARN_FLAG}" = "Xyes"; then
echo "
***************
*   WARNING   *
***************

At least one warning was issued. Some features were disabled.

Review the configure output and/or the contents of config.log 
before proceeding with the build.
"
fi

echo "
------------------------------
SUNDIALS Configuration Summary
------------------------------"

echo "
Configuration
-------------

  Host System:               ${host}
  Build System:              ${build}

  C Preprocessor:            ${CPP} 
  C Preprocessor Flags:      ${CPPFLAGS}
  C Compiler:	             ${CC}
  C Compiler Flags           ${CFLAGS}
  C Linker:                  ${CC}
  Linker Flags:              ${LDFLAGS}
  Libraries:                 ${LIBS}"

if test "X${F77_OK}" = "Xyes"; then
echo "
  Fortran Compiler:          ${F77}
  Fortran Compiler Flags:    ${FFLAGS}
  Fortran Linker:            ${F77_LNKR}
  Extra Fortran Libraries:   ${FLIBS}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_C_COMP_OK}" = "Xyes"; then
echo "
  MPI Root Directory:        ${MPI_ROOT_DIR}
  MPI Include Directory:     ${MPI_INC_DIR}
  MPI Library Directory:     ${MPI_LIB_DIR}
  MPI Flags:                 ${MPI_FLAGS}
  Extra MPI Libraries:       ${MPI_LIBS}

  Using MPI-C script?        ${USE_MPICC_SCRIPT}
  MPI-C:                     ${MPICC}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${F77_EXAMPLES_ENABLED}" = "Xyes" && test "X${MPI_F77_COMP_OK}" = "Xyes"; then
echo "
  Using MPI-Fortran script?  ${USE_MPIF77_SCRIPT}
  MPI-Fortran:               ${MPIF77}
  MPI-Fortran Linker:        ${MPIF77_LNKR}"
fi

# Determine SOURCE, BUILD, and EXEC_PREFIX directories
cv_srcdir=`( cd ${srcdir} ; pwd )`
cv_builddir=`pwd`
if test "X${exec_prefix}" = "XNONE"; then
  cv_exec_prefix="${prefix}"
else
  cv_exec_prefix="${exec_prefix}"
fi

echo "
  srcdir:                    ${cv_srcdir}
  builddir:                  ${cv_builddir}
  prefix:                    ${prefix}
  exec_prefix:               ${cv_exec_prefix}
  includedir:                ${includedir}
  libdir:                    ${libdir}"

if test "X${EXAMPLES_ENABLED}" = "Xyes"; then
echo "  examples installed in:     ${EXS_INSTDIR}"
fi

echo "
Modules
-------
"

if test "X${CVODE_ENABLED}" = "Xyes"; then
  THIS_LINE="CVODE"
  if test "X${FCMIX_ENABLED}" = "Xyes"; then
    THIS_LINE="${THIS_LINE}  FCVODE"
  fi
  echo "  ${THIS_LINE}"
fi

if test "X${CVODES_ENABLED}" = "Xyes"; then
  THIS_LINE="CVODES"
  echo "  ${THIS_LINE}"
fi

if test "X${IDA_ENABLED}" = "Xyes"; then
  THIS_LINE="IDA"
  if test "X${FCMIX_ENABLED}" = "Xyes"; then
    THIS_LINE="${THIS_LINE}    FIDA"
  fi
  echo "  ${THIS_LINE}"
fi

if test "X${IDAS_ENABLED}" = "Xyes"; then
  THIS_LINE="IDAS"
  echo "  ${THIS_LINE}"
fi

if test "X${KINSOL_ENABLED}" = "Xyes"; then
  THIS_LINE="KINSOL"
  if test "X${FCMIX_ENABLED}" = "Xyes"; then
    THIS_LINE="${THIS_LINE} FKINSOL"
  fi
  echo "  ${THIS_LINE}"
fi

if test "X${CPODES_ENABLED}" = "Xyes"; then
  THIS_LINE="CPODES"
  echo "  ${THIS_LINE}"
fi

if test "X${EXAMPLES_ENABLED}" = "Xyes"; then
echo "
Examples
--------
"

echo "  Serial C examples:         ${SERIAL_C_EXAMPLES}"
echo "  Parallel C examples:       ${PARALLEL_C_EXAMPLES}"
echo "  Serial Fortran examples:   ${SERIAL_F77_EXAMPLES}"
echo "  Parallel Fortran examples: ${PARALLEL_F77_EXAMPLES}"

fi


echo "  
  Type 'make' and then 'make install' to build and install ${PACKAGE_STRING}."



echo "
----------------------------------
Finished SUNDIALS Configure Script
----------------------------------
"

]) dnl END SUNDIALS_REPORT
