# -----------------------------------------------------------------
# $Revision: 1.48 $
# $Date: 2006-11-08 00:48:19 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# SUNDIALS autoconf macros
# -----------------------------------------------------------------

#------------------------------------------------------------------
# GREETING
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_GREETING],
[

# Say Hi!
if test "X${SUNDIALS_DEV}" = "Xyes"; then
echo "
---------------------------------------------
Running SUNDIALS Development Configure Script
---------------------------------------------
"
else
echo "
---------------------------------
Running SUNDIALS Configure Script
---------------------------------
"
fi

]) dnl END SUNDIALS_GREETING

#------------------------------------------------------------------
# PERFORM INITIALIZATIONS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_INITIALIZE],
[

# Specify directory containing auxillary build tools and M4 files
AC_CONFIG_AUX_DIR([config])

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

# Test if various environment variables exist.
# If not, we may set some default values for them
if test "X${CPPFLAGS}" = "X"; then
  CPPFLAGS_PROVIDED="no"
else
  CPPFLAGS_PROVIDED="yes"
fi
if test "X${CFLAGS}" = "X"; then
  CFLAGS_PROVIDED="no"
else
  CFLAGS_PROVIDED="yes"
fi
if test "X${LDFLAGS}" = "X"; then
  LDFLAGS_PROVIDED="no"
else
  LDFLAGS_PROVIDED="yes"
fi
if test "X${FFLAGS}" = "X"; then
  FFLAGS_PROVIDED="no"
else
  FFLAGS_PROVIDED="yes"
fi

# Set defaults for config/sundials_config.in file
F77_MANGLE_MACRO1=""
F77_MANGLE_MACRO2=""
F77_CASE=""
F77_UNDERSCORES=""
PRECISION_LEVEL=""
GENERIC_MATH_LIB=""
F77_MPI_COMM_F2C=""

# Initialize enable status of various modules, options, and features
# to their default values

CVODE_ENABLED="yes"
CVODES_ENABLED="yes"
IDA_ENABLED="yes"
#IDAS_ENABLED="yes"
IDAS_ENABLED="no"
KINSOL_ENABLED="yes"
#CPODES_ENABLED="yes"
CPODES_ENABLED="no"
FCMIX_ENABLED="yes"
MPI_ENABLED="yes"
STB_ENABLED="no"
EXAMPLES_ENABLED="no"
F77_EXAMPLES_ENABLED="no"
DEV_EXAMPLES_ENABLED="no"

# Initialize variables that may NOT necessarily be initialized
# during normal execution. Should NOT use uninitialized variables
F77_OK="no"
MPI_C_COMP_OK="no"
MPI_F77_COMP_OK="no"
STB_PARALLEL_OK="no"

# This variable is set to "yes" if an AC_MSG_WARN statement
# was executed
SUNDIALS_WARN_FLAG="no"

]) dnl END SUNDIALS_INITIALIZE

#------------------------------------------------------------------
# FIND ARCHIVER (not used)
#------------------------------------------------------------------

# Check for archiver (ar)
# Define AR="ar rc" if found
# AC_SUBST is called automatically for AR
AC_DEFUN([SUNDIALS_CHECK_AR],
[

AC_CHECK_TOOL([AR],[ar],["none"])

if test "X${AR}" = "Xnone"; then
  AC_MSG_ERROR([cannot find archiver])
else
  AR="${AR} rc"
fi

]) dnl END SUNDIALS_CHECK_AR

#------------------------------------------------------------------
# FIND ARCHIVE INDEX UTILITY (not used)
#------------------------------------------------------------------

# Check for ranlib utility
# Defines RANLIB="ranlib" if found
# AC_SUBST is called automatically for RANLIB
AC_DEFUN([SUNDIALS_CHECK_RANLIB],
[

AC_CHECK_TOOL([RANLIB],[ranlib],["none"])

if test "X${RANLIB}" = "Xnone"; then
  AC_MSG_WARN([cannot find ranlib])
  SUNDIALS_WARN_FLAG="yes"
fi

]) dnl END SUNDIALS_CHECK_RANLIB

#------------------------------------------------------------------
# TEST ENABLES
#
# The following variables may be changed here (default value in []):
#
#   CVODE_ENABLED    - enable CVODE module [yes]
#   CVODES_ENABLED   - enable CVODES module [yes]
#   IDA_ENABLED      - enable IDA module [yes]
#   KINSOL_ENABLED   - enable KINSOL module [yes]
#   FCMIX_ENABLED    - enable Fortran-C inerfaces [yes]
#   MPI_ENABLED      - enable parallel support [yes]
#   EXAMPLES_ENABLED - enable example programs [no]
#   STB_ENABLED      - enable sundialsTB Matlab interfaces [no]
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
#AC_ARG_ENABLE(ida,
#[AC_HELP_STRING([--disable-idas],[disable configuration of IDAS])],
#[
#if test "X${enableval}" = "Xno"; then
#  IDAS_ENABLED="no"
#fi
#],
#[
#if test -d ${srcdir}/src/idas  ; then
#  IDAS_ENABLED="yes"
#else
#  IDAS_ENABLED="no"
#fi
#])

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
#AC_ARG_ENABLE(cpodes,
#[AC_HELP_STRING([--disable-cpodes],[disable configuration of CPODES])],
#[
#if test "X${enableval}" = "Xno"; then
#  CPODES_ENABLED="no"
#fi
#],
#[
#if test -d ${srcdir}/src/cpodes ; then
#  CPODES_ENABLED="yes"
#else
#  CPODES_ENABLED="no"
#fi
#])

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

# Check if user wants to disable support for MPI.
# If not, then make certain source directory actually exists
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

# Check if user wants to enable all examples
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

# Check if user wants to enable support for MEX compilation
# sundialsTB is NOT built by default
# If yes, then make certain source directory actually exists
AC_ARG_ENABLE([sundialsTB],
[AC_HELP_STRING([--enable-sundialsTB],[enable configuration of sundialsTB])],
[
if test "X${enableval}" = "Xno"; then
  STB_ENABLED="no"
else
  if test -d ${srcdir}/sundialsTB ; then
    STB_ENABLED="yes"
  else
    STB_ENABLED="no"
  fi
fi
])

# Fortran examples are enabled only if both FCMIX and EXAMPLES are enabled
if test "X${FCMIX_ENABLED}" = "Xyes" && test "X${EXAMPLES_ENABLED}" = "Xyes"; then
  F77_EXAMPLES_ENABLED="yes"
fi

]) dnl END SUNDIALS_ENABLES

#------------------------------------------------------------------
# TEST DEVELOPMENT ENABLES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEV_ENABLES],
[

IDAS_ENABLED="yes"
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
if test -d ${srcdir}/src/idas ; then
  IDAS_ENABLED="yes"
else
  IDAS_ENABLED="no"
fi
])

CPODES_ENABLED="yes"
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


]) dnl END SUNDIALS_DEV_ENABLES

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

  SUNDIALS_CHECK_CC

fi

]) dnl END SUNDIALS_SET_CC
 

AC_DEFUN([SUNDIALS_CHECK_CC],
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

# If CFLAGS is not provided, set defaults
if test "X${CFLAGS_PROVIDED}" = "Xno"; then
  SUNDIALS_DEFAULT_CFLAGS
fi

# Add any additional CFLAGS
AC_MSG_CHECKING([for extra C compiler flags])
AC_ARG_WITH(cflags,
[AC_HELP_STRING([--with-cflags=ARG],[add extra C compiler flags])],
[
AC_MSG_RESULT([${withval}])
CFLAGS="${CFLAGS} ${withval}"
],
[
AC_MSG_RESULT([none])
])

# Set CPP to command that runs C preprocessor
AC_PROG_CPP

# If CPPFLAGS is not provided, set defaults
if test "X${CPPFLAGS_PROVIDED}" = "Xno"; then
  SUNDIALS_DEFAULT_CPPFLAGS
fi

# Add any additional CPPFLAGS
AC_MSG_CHECKING([for extra C/C++ preprocessor flags])
AC_ARG_WITH(cppflags,
[AC_HELP_STRING([--with-cppflags=ARG],[add extra C/C++ preprocessor flags])],
[
AC_MSG_RESULT([${withval}])
CPPFLAGS="${CPPFLAGS} ${withval}"
],
[
AC_MSG_RESULT([none])
])

# If LDFLAGS is not provided, set defaults
if test "X${LDFLAGS_PROVIDED}" = "Xno"; then
  SUNDIALS_DEFAULT_LDFLAGS
fi

# Add any additional linker flags 
AC_MSG_CHECKING([for extra linker flags])
AC_ARG_WITH(ldflags,
[AC_HELP_STRING([--with-ldflags=ARG],[add extra linker flags])],
[
AC_MSG_RESULT([${withval}])
LDFLAGS="${LDFLAGS} ${withval}"
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

#------------------------------------------------------------------
# DEFAULT CFLAGS
#
# Set default CFLAGS (called only if CFLAGS is not provided)
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_CFLAGS],
[

# Extract executable name of C compiler, for optional use below
CC_EXEC=`basename "${CC}"`

case $host in 

  # IA-32 system running Linux
  i*-pc-linux-*)
    #if test "X${GCC}" = "Xyes"; then
    #  if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
    #    CFLAGS="-ffloat-store"
    #  fi
    #fi
    ;;

  *)
    CFLAGS=""
    ;;

esac

]) dnl END SUNDIALS_DEFAULT_CFLAGS

#------------------------------------------------------------------
# DEFAULT CPPFLAGS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_CPPFLAGS],
[

# Set default CPPFLAGS (called only if CPPFLAGS is not provided)

CPPFLAGS=""

]) dnl END SUNDIALS_DEFAULT_CFLAGS

#------------------------------------------------------------------
# DEFAULT LDFLAGS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_LDFLAGS],
[

# Set default LDFLAGS (called only if LDFLAGS is not provided)

case $build_os in

  *cygwin*)
     # Under cygwin, if building shared libraries, add -no-undefined
     # and explicitly specify library extension
     if test "X${enable_shared}" = "Xyes"; then
       LDFLAGS="-no-undefined"
     fi
     ;;

  *mingw*)
     # Under mingw, if building shared libraries, add -no-undefined
     # and explicitly specify library extension
     if test "X${enable_shared}" = "Xyes"; then
       LDFLAGS="-no-undefined"
     fi
     ;;


  *)
     LDFLAGS=""
     ;;

esac

]) dnl END SUNDIALS_DEFAULT_CFLAGS

#------------------------------------------------------------------
# CHECK FORTRAN COMPILER
#
# If a working compiler cannot be found, set F77_OK="no" and issue
# a warning that Fortran support will be disabled.
#
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_F77],
[

if test "X${F77}" = "X"; then

  AC_MSG_WARN([cannot find a Fortran compiler])
  echo ""
  echo "   Unable to find a working Fortran compiler"
  echo ""
  echo "   Try using F77 to explicitly specify a C compiler"
  echo ""
  echo "   Disabling F77 support..."
  echo ""

  F77_OK="no"
  F77_EXAMPLES_ENABLED="no"
  SUNDIALS_WARN_FLAG="yes"

else

  F77_OK="yes"
  SUNDIALS_CHECK_F77

fi

]) dnl END SUNDIALS_SET_F77


AC_DEFUN([SUNDIALS_CHECK_F77],
[

AC_LANG_PUSH([Fortran 77])

# If FFLAGS is not provided, set defaults
if test "X${FFLAGS_PROVIDED}" = "Xno"; then
  SUNDIALS_DEFAULT_FFLAGS
fi        

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

# Determine Fortran name mangling scheme
# Default is lower case with one underscore
AC_MSG_CHECKING([Fortran name-mangling scheme])
RUN_F77_WRAPPERS="yes"
F77_WRAPPER_CHECK_OK="yes"
F77_FUNC_CASE="lower"
F77_FUNC_UNDERSCORES="one"

# If user provides the number of underscores, overwrite default value
AC_ARG_WITH(f77underscore, 
[AC_HELP_STRING([--with-f77underscore=ARG],[specify number of underscores to append to function names (none/one/two) [AUTO]],[])],
[
if test "X${withval}" = "Xnone" || test "X${withval}" = "Xone" || test "X${withval}" = "Xtwo"; then
  F77_FUNC_UNDERSCORES="${withval}"
else
  AC_MSG_RESULT([failed])
  AC_MSG_ERROR([invalid input for --with-f77underscore])
fi
RUN_F77_WRAPPERS="no"
])

# If user provides the case, overwrite default value
AC_ARG_WITH(f77case, 
[AC_HELP_STRING([--with-f77case=ARG   ],[specify case of function names (lower/upper) [AUTO]],[])],
[
if test "X${withval}" = "Xupper" || test "X${withval}" = "Xlower"; then
  F77_FUNC_CASE="${withval}"
else
  AC_MSG_RESULT([failed])
  AC_MSG_ERROR([invalid input for --with-f77case])
fi
RUN_F77_WRAPPERS="no"
])

# Determine how to properly mangle function names so Fortran subroutines can
# call C functions included in SUNDIALS libraries
if test "X${RUN_F77_WRAPPERS}" = "Xyes"; then
  SUNDIALS_F77_WRAPPERS
fi

# Based on F77_FUNC_CASE and F77_FUNC_UNDERSCORES, define C preprocessor macros
# F77_FUNC and F77_FUNC_
if test "X${F77_WRAPPER_CHECK_OK}" = "Xyes"; then
  AC_MSG_RESULT([case: ${F77_FUNC_CASE}, underscores: ${F77_FUNC_UNDERSCORES}])
  if test "X${F77_FUNC_CASE}" = "Xlower"; then
    if test "X${F77_FUNC_UNDERSCORES}" = "Xnone"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[name],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) name"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) name"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xone"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[name ## _],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) name ## _"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) name ## _"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xtwo"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[name ## __],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) name ## __"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) name ## __"
    fi
  elif test "X${F77_FUNC_CASE}" = "Xupper"; then
    if test "X${F77_FUNC_UNDERSCORES}" = "Xnone"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[NAME],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) NAME"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) NAME"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xone"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[NAME ## _],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) NAME ## _"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) NAME ## _"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xtwo"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[NAME ## __],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) NAME ## __"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) NAME ## __"
    fi
  fi
fi

# Check if we must use a Fortran compiler to link the Fortran examples
# (default is to use CC)
AC_MSG_CHECKING([which linker to use])
SUNDIALS_F77_LNKR_CHECK
AC_MSG_RESULT([${F77_LNKR_OUT}])

AC_LANG_POP([Fortran 77])

]) dnl END SUNDIALS_SET_F77

#------------------------------------------------------------------
# DETERMINE F77 LINKER
#
# Set F77_LNKR variable
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_LNKR_CHECK],
[

# Perform the next test only if using a C++ compiler to build SUNDIALS.
# Note: SUNDIALS_CPLUSPLUS_CHECK macro was called by SUNDIALS_SET_CC,
#       so we need only check the value of USING_CPLUSPLUS_COMP.
if test "X${USING_CPLUSPLUS_COMP}" = "Xyes"; then

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

  # Note: F77_LNKR variable is used by serial Fortran example Makefile's
  if test "X${F77_LNKR_CHECK_OK}" = "Xyes"; then
    F77_LNKR="\$(CC)"
    F77_LNKR_OUT="${CC}"
  else
    F77_LNKR="\$(F77)"
    F77_LNKR_OUT="${F77}"
  fi

  ],
  [

  # If compilation fails, then just use F77
  F77_LNKR="\$(F77)"
  F77_LNKR_OUT="${F77}"

  ])

else

  # Using a C compiler so use F77 to link serial Fortran examples
  F77_LNKR="\$(F77)"
  F77_LNKR_OUT="${F77}"

fi

]) dnl SUNDIALS_F77_LNKR_CHECK

#------------------------------------------------------------------
# DETERMINE F77 NAME-MANGLING SCHEME
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_WRAPPERS],
[

# Replacement macro for AC_F77_WRAPPERS which has too many problems
# (based on implementation of AC_F77_WRAPPERS minus the "fluff")

# Provide variable description templates for config.hin and config.h files
# Required by autoheader utility
AH_TEMPLATE(F77[_FUNC],
            [FCMIX: Define name-mangling macro])

# Remaining test pertains to Fortran programming language
AC_LANG_PUSH([Fortran 77])

# Compile a dummy Fortran subroutine
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

F77_WRAPPER_CHECK_OK="no"
for i in "sundials" "SUNDIALS"
do
  for j in "" "_" "__"
  do
    F77_MANGLED_NAME="${i}${j}"
    AC_LINK_IFELSE([AC_LANG_CALL([],[${F77_MANGLED_NAME}])],[F77_WRAPPER_CHECK_OK="yes" ; break 2])
  done
done

AC_LANG_POP([C])

# If test succeeded, then set F77_FUNC_CASE and F77_FUNC_UNDERSCORES variables
if test "X${F77_WRAPPER_CHECK_OK}" = "Xyes"; then
  # Determine case (lower, upper)
  if test "X${i}" = "Xsundials"; then
    F77_FUNC_CASE="lower"
  else
    F77_FUNC_CASE="upper"
  fi
  # Determine number of underscores to append (none, one, two)
  if test "X${j}" = "X"; then
    F77_FUNC_UNDERSCORES="none"
  elif test "X${j}" = "X_"; then
    F77_FUNC_UNDERSCORES="one"
  else
    F77_FUNC_UNDERSCORES="two"
  fi
# If test failed, then tell user to use '--with-f77case' and '--with-f77underscore'
# and disable Fortran support
else
  AC_MSG_RESULT([UNKNOWN])
  AC_MSG_WARN([cannot determine Fortran name-mangling scheme])
  echo ""
  echo "   Unable to determine name-mangling scheme required by Fortran"
  echo "   compiler."
  echo ""
  echo "   Try using --with-f77case and --with-f77underscore to explicitly"
  echo "   specify the appropriate name-mangling scheme."
  echo ""
  echo "   Disabling Fortran support..."
  echo ""
  F77_OK="no"
  F77_EXAMPLES_ENABLED="no"
  SUNDIALS_WARN_FLAG="yes"
fi

# Set LIBS environment variable back to original value
LIBS="${SAVED_LIBS}"

],
[

# If a compilation error occurred, then disable Fortran support
AC_MSG_RESULT([ERROR])
echo ""
echo "   Unable to compile test program using given Fortran compiler."
echo ""
echo "   Disabling Fortran support..."
echo ""
F77_OK="no"
F77_EXAMPLES_ENABLED="no"
SUNDIALS_WARN_FLAG="yes"

])

# Reset language (remove 'Fortran 77' from stack)
AC_LANG_POP([Fortran 77])

# Remove temporary file
rm -f f77_wrapper_check.${ac_objext}

]) dnl END SUNDIALS_SET_F77

#------------------------------------------------------------------
# DEFAULT FFLAGS
#
# Set default FFLAGS (called only if FFLAGS is not provided)
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_FFLAGS],
[

# Extract executable name of Fortran compiler, for optional use below
F77_EXEC=`basename "${F77}"`

case $host in

  # IA-32 system running Linux
  i*-pc-linux-*)
    #if test "X${G77}" = "Xyes"; then
    #  if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
    #    FFLAGS="-ffloat-store"
    #  fi
    #fi
    ;;

  # SGI/IRIX
  mips-sgi-irix*) 
    FFLAGS="-64"
    ;;

  # Compaq/Tru64
  *-dec-osf*)
    if test "X${F77_EXEC}" = "Xf77"; then
      FFLAGS="-O1"
    fi
    ;;

  *)
    FFLAGS=""
    ;;

esac

]) dnl END SUNDIALS_DEFAULT_FFLAGS

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
[AC_HELP_STRING([--with-mpicc[[[[=ARG]]]]],[specify MPI-C compiler to use @<:@mpicc@:>@],
                [                                ])],
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
  SUNDIALS_CHECK_CC_WITH_MPI
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

AC_DEFUN([SUNDIALS_CHECK_CC_WITH_MPI],
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

]) dnl END SUNDIALS_CHECK_CC_WITH_MPI

#------------------------------------------------------------------
# CHECK MPI-F77 COMPILER
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
[AC_HELP_STRING([--with-mpif77[[[[=ARG]]]]],[specify MPI-Fortran compiler to use @<:@mpif77@:>@],
                [                                ])],
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
  SUNDIALS_CHECK_F77_WITH_MPI
fi

]) dnl END SUNDIALS_SET_MPIF77

#------------------------------------------------------------------
# DETERMINE MPI-FORTRAN LINKER
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

  # Note: MPIF77_LNKR variable is used by parallel Fortran example Makefile's
  if test "X${MPIF77_LNKR_CHECK_OK}" = "Xyes"; then
    MPIF77_LNKR="\$(MPICC)"
    MPIF77_LNKR_OUT="${MPICC}"
  else
    MPIF77_LNKR="\$(MPIF77)"
    MPIF77_LNKR_OUT="${MPIF77}"
  fi

  ],
  [

  # If compilation fails, then just use MPIF77
  MPIF77_LNKR="\$(MPIF77)"
  MPIF77_LNKR_OUT="${MPIF77}"

  ])

else

  # Using a C compiler so use MPIF77 to link parallel Fortran examples
  MPIF77_LNKR="\$(MPIF77)"
  MPIF77_LNKR_OUT="${MPIF77}"

fi
AC_MSG_RESULT([${MPIF77_LNKR_OUT}])


]) dnl SUNDIALS_MPIF77_LNKR_CHECK

#------------------------------------------------------------------
# TEST MPI-FORTRAN COMPILER
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
# TEST FORTRAN COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_F77_WITH_MPI],
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
# SUNDIALS_CHECK_CC_WITH_MPI has been executed
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
    if test "X${F77_LNKR}" = "X\$(CC)"; then
      MPIF77_LNKR="\$(MPICC)"
      MPIF77_LNKR_OUT="${MPICC}"
    elif test "X${F77_LNKR}" = "X\$(F77)"; then
      MPIF77_LNKR="\$(MPIF77)"
      MPIF77_LNKR_OUT="${MPIF77}"
    fi
    AC_MSG_RESULT([${MPIF77_LNKR_OUT}])
  fi

else
  MPI_F77_COMP_OK="no"
fi

# Restore FFLAGS, LDFLAGS and LIBS
FFLAGS="${SAVED_FFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

AC_LANG_POP([Fortran 77])

]) dnl END SUNDIALS_CHECK_F77_WITH_MPI

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
# Test for MPI_Comm_spawn() funciton:
#   (1) NO  : sundialsTB will NOT provide parallel support
#   (2) YES : sundialsTB will be confoigured with parallel support
 
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

  # Check if MPI implementation supports MPI_Comm_spawn() from
  # MPI-2 specification
  if test "X${STB_ENABLED}" = "Xyes"; then
    AC_MSG_CHECKING([for MPI_Comm_spawn() from MPI-2 specification])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include "mpi.h"]],
    [[
        int c;
        char **v;
        MPI_Info info;
        MPI_Comm comm;
        MPI_Init(&c, &v);
        c = MPI_Comm_spawn(*v, v, c, info, c, comm, &comm, &c);
        MPI_Finalize();
    ]])],
    [STB_PARALLEL_OK="yes"],
    [STB_PARALLEL_OK="no"])
     AC_MSG_RESULT([${STB_PARALLEL_OK}])
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

#------------------------------------------------------------------
# CHECK MATLAB
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_MATLAB],
[

AC_REQUIRE([AC_CANONICAL_HOST])

# Set the STB_OS variable, depending on the build OS
case $build_os in
  *cygwin*)
    STB_OS="cygwin"
    ;;
  *mingw*)
    STB_OS="mingw"
    ;;
  *) 
    STB_OS="other"
    ;;
esac

# Under cygwin, check if CFLGAGS contains -mno-cygwin 
if test "X${STB_OS}" = "Xcygwin"; then
  AC_MSG_CHECKING([if CFLAGS contains -mno-cygwin])
  cv_no_cygwin_exists="no"
  for i in ${CFLAGS}
  do
    if test "X${i}" = "X-mno-cygwin"; then
      cv_no_cygwin_exists="yes"
    fi
  done
  AC_MSG_RESULT([${cv_no_cygwin_exists}])
  if test "X${cv_no_cygwin_exists}" = "Xno"; then
    AC_MSG_WARN([compilation of sundialsTB mex files may fail])
    echo ""
    echo "   Building under CYGWIN without -mno-cygwin"
    echo ""
    echo "   Beware that compilation of the sundialsTB mex files may fail."
    echo ""
    SUNDIALS_WARN_FLAG="yes"
  fi 
fi

AC_ARG_WITH([],[      ],[])

# Find Matlab program. If --with-matlab=<matlab> was passed to the configure 
# script, then we only check if it exists and is an executable file. We'll try
# to run it later... Otherwise, use AC_PATH_PROG to search for matlab.
MATLAB_CMD="none"
AC_ARG_WITH([matlab], 
[AC_HELP_STRING([--with-matlab=MATLAB], [specify Matlab executable])],
[
AC_MSG_CHECKING([for matlab])
if test -f ${withval} ; then
  AC_MSG_RESULT([${withval}])
  MATLAB_CMD="${withval}"
else
  AC_MSG_RESULT([none])
  AC_MSG_WARN([invalid value '${withval}' for --with-matlab])
  echo ""
  echo "   ${withval} does not exist"
  echo ""
  echo "   Disabling compilation of sundialsTB mex files..."
  echo ""
  SUNDIALS_WARN_FLAG="yes"
  STB_ENABLED="no"
fi
],
[
AC_PATH_PROG(MATLAB_CMD, matlab, "none")
if test "X${MATLAB_CMD}" = "Xnone"; then
  STB_ENABLED="no"
fi
])

MATLAB_CMD_FLAGS="-nojvm -nosplash"

# Set MEXOPTS (MEX options file)
if test "X${STB_ENABLED}" = "Xyes"; then

  MEXOPTS=""
  AC_ARG_WITH([mexopts], 
  AC_HELP_STRING([--with-mexopts=ARG], [use MEX options file ARG [[standard]]]),
  [
  cv_mexopts_file="${withval}"
  if test -f ${cv_mexopts_file} ; then
    cv_mexopts_dir=`AS_DIRNAME(["${cv_mexopts_file}"])`
    # MEX options file is located under the current working directory
    if test "X${cv_mexopts_dir}" = "X${cv_mexopts_file}"; then
      cv_mexopts_dir="."
      cv_mexopts_name="${cv_mexopts_file}"
      cv_mexopts_file="${cv_mexopts_dir}/${cv_mexopts_name}"
    fi
    MEXOPTS="-f ${cv_mexopts_file}"
  else
    AC_MSG_WARN([invalid value '${cv_mexopts_file}' for --with-mexopts])
    echo ""
    echo "   ${cv_mexopts_file} does not exist"
    echo ""
    echo "   Disabling compilation of sundialsTB mex files..."
    echo ""
    SUNDIALS_WARN_FLAG="yes"
    STB_ENABLED="no"
  fi
  ])

fi

# Set MEXFLAGS and MEXLDADD
if test "X${STB_ENABLED}" = "Xyes"; then

  AC_MSG_CHECKING([for MEX compiler compiler flags])
  AC_ARG_WITH([mexflags], 
  [AC_HELP_STRING([--with-mexflags=ARG], [specify MEX compiler flags])],
  [
  AC_MSG_RESULT([${withval}])
  MEXFLAGS="${withval}"
  ],
  [
  # If MEXFLAGS is not defined, set it to the default -O
  if test "X${MEXFLAGS}" = "X"; then
    MEXFLAGS="-O"
  fi
  AC_MSG_RESULT([${MEXFLAGS}])
  ])

  AC_MSG_CHECKING([for additional MEX linker flags])
  AC_ARG_WITH([mexldadd], 
  [AC_HELP_STRING([--with-mexldadd=ARG], [specify additional MEX linker flags])],
  [
  AC_MSG_RESULT([${withval}])
  MEXLDADD="${withval}"
  ],
  [
  # If MEXLDADD is not defined, none are used
  if test "X${MEXLDADD}" = "X"; then
    MEXLDADD=""
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${MEXLDADD}])
  fi
  ])

fi

# Run matlab and try the MEX compiler. Set extension for MEX files (set MEXEXT)
if test "X${STB_ENABLED}" = "Xyes"; then

  AC_MSG_CHECKING([if the Matlab MEX compiler works])

  # Create test file cvmextest.c (temporary file)
  cat > cvmextest.c <<_END_MEX_C
#include "mex.h"
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{}
_END_MEX_C

  # Create test file cvmextest_script.m (temporary file)
  cat > cvmextest_script.m <<_END_MEX_M
mex ${MEXOPTS} ${MEXFLAGS} -output cvmextest cvmextest.c ${MEXLDADD}
exit
_END_MEX_M

  # Run matlab in batch mode
  # Warning: File descriptors and I/O redirection can be problematic
  ( eval ${MATLAB_CMD} ${MATLAB_CMD_FLAGS} -r cvmextest_script ) 2>/dev/null 1>&2

  # Get exit status of previous command (meaning eval statement)
  cv_status=$?

  # MEX test succeeded
  if test "${cv_status}" = "0"; then
     AC_MSG_RESULT([yes])
     AC_MSG_CHECKING([for MEX file extension])
     if test -f cvmextest.dll ; then
        MEXEXT="dll"
     elif test -f cvmextest.mex ; then
        MEXEXT="mex"
     elif test -f cvmextest.mexaxp ; then
        MEXEXT="mexaxp"
     elif test -f cvmextest.mexglx ; then
        MEXEXT="mexglx"
     elif test -f cvmextest.mexhp7 ; then
        MEXEXT="mexhp7"
     elif test -f cvmextest.mexhpux ; then
        MEXEXT="mexhpux"
     elif test -f cvmextest.mexrs6 ; then
        MEXEXT="mexrs6"
     elif test -f cvmextest.mexsg ; then
        MEXEXT="mexsg"
     elif test -f cvmextest.mexsol ; then
        MEXEXT="mexsol"
     else
	MEXEXT=""
     fi
     AC_MSG_RESULT([${MEXEXT}])
  # MEX test failed
  else 
     AC_MSG_RESULT([no])
     STB_ENABLED="no"
  fi

  # Remove temporary files
  rm -f cvmextest*

fi

# Determine where to install sundialsTB
if test "X${STB_ENABLED}" = "Xyes"; then

  AC_ARG_WITH([],[           ],[])
  AC_MSG_CHECKING([where to install sundialsTB])
  AC_ARG_WITH([sundialsTB-instdir],
  [AC_HELP_STRING([--with-sundialsTB-instdir=STBINSTDIR], [install sundialsTB in STBINSTDIR @<:@MATLAB/toolbox@:>@])],
  [
    STB_INSTDIR="${withval}"
  ],
  [
    if test "X${MATLAB}" = "X"; then
      STB_INSTDIR="no"
    else
      STB_INSTDIR="${MATLAB}/toolbox/sundialsTB"
    fi
  ])

  # Set STB_PATH (usually same as STB_INSTDIR)
  if test "X${STB_OS}" = "Xcygwin"; then
    STB_PATH=`cygpath -a -m "${STB_INSTDIR}"`
  elif test "X${STB_OS}" = "Xmingw"; then
    STB_PATH=`cd "${STB_INSTDIR}" > /dev/null 2>&1 && pwd -W`
  else
    STB_PATH="${STB_INSTDIR}"
  fi

  AC_MSG_RESULT([${STB_PATH}])

fi

AC_ARG_WITH([],[        ],[])

]) dnl END SUNDIALS_SET_MEX

#------------------------------------------------------------------
# ADD SOME MORE STUFF TO configure --help
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_MORE_HELP],
[

AC_ARG_WITH([],[          ],[])
AC_ARG_WITH([],[NOTES],[])
AC_ARG_WITH([],[  Both --with-examples-instdir and --with-sundialsTB-instdir can be set to "no",],[])
AC_ARG_WITH([],[  in which case the examples and sundialsTB, respectively, are built but not installed.],[])
AC_ARG_WITH([],[  Enabling the compilation of the examples (--enable-examples) but disabling their],[])
AC_ARG_WITH([],[  installation (--with-examples-instdir=no) can be used to test the SUNDIALS libraries.],[])

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

# Check if the Fortran update script (config/fortran_update.in) is needed
if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" || test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
  BUILD_F77_UPDATE_SCRIPT="yes";
else
  BUILD_F77_UPDATE_SCRIPT="no"
fi

# Where should we install the examples?
AC_MSG_CHECKING([where to install the SUNDIALS examples])
AC_ARG_WITH([],[   ],[])
AC_ARG_WITH([examples-instdir],
[AC_HELP_STRING([--with-examples-instdir=EXINSTDIR], [install SUNDIALS examples in EXINSTDIR @<:@EPREFIX/examples@:>@])],
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

# Development examples are disabled
SERIAL_DEV_C_EXAMPLES="disabled"
PARALLEL_DEV_C_EXAMPLES="disabled"
SERIAL_DEV_F77_EXAMPLES="disabled"
PARALLEL_DEV_F77_EXAMPLES="disabled"

]) dnl END SUNDIALS_SET_EXAMPLES

#------------------------------------------------------------------
# SET DEVELOPER EXAMPLES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_DEV_EXAMPLES],
[

DEV_EXAMPLES_ENABLED="${EXAMPLES_ENABLED}"
SERIAL_DEV_C_EXAMPLES="${SERIAL_C_EXAMPLES}"
PARALLEL_DEV_C_EXAMPLES="${PARALLEL_C_EXAMPLES}"
SERIAL_DEV_F77_EXAMPLES="${SERIAL_F77_EXAMPLES}"
PARALLEL_DEV_F77_EXAMPLES="${PARALEL_F77_EXAMPLES}"

]) dnl END SUNDIALS_SET_DEV_EXAMPLES

#------------------------------------------------------------------
# BUILD MODULES LIST
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_BUILD_MODULES_LIST],
[

# Initialize the list of Makefiles to be created
SUNDIALS_MAKEFILES="Makefile"

# Initialize list of additional configure files to be created
SUNDIALS_CONFIGFILES="src/sundials/sundials_config.h:src/sundials/sundials_config.in"
SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} sundials-config:config/sundials-config.in"

# Initialize lists of solver modules, example modules, and sundialsTB modules
SLV_MODULES="src/sundials"
if test "X${SUNDIALS_DEV}" = "Xyes" && test -f ${srcdir}/src/sundials/Makefile.dev.in ; then
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/sundials/Makefile:src/sundials/Makefile.dev.in"
else
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/sundials/Makefile"
fi

EXS_MODULES=""
STB_MODULES=""

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
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/fcmix_serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/fcmix_serial/Makefile"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvode/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvode/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/cvode/serial/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/parallel/Makefile"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/fcmix_parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvode/fcmix_parallel/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvode/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvode/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/cvode/parallel/Makefile"
  fi

fi

# CVODES module
if test "X${CVODES_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/cvodes"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/cvodes/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvodes/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvodes/serial/Makefile"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvodes/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvodes/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/cvodes/serial/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvodes/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cvodes/parallel/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvodes/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/cvodes/parallel/Makefile"
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
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/fcmix_serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/fcmix_serial/Makefile"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/ida/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/ida/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/ida/serial/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/parallel/Makefile"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/fcmix_parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/ida/fcmix_parallel/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/ida/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/ida/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/ida/parallel/Makefile"
  fi

fi

# IDAS module
if test "X${IDAS_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/idas"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/idas/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/idas/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/idas/serial/Makefile"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/idas/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/idas/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/idas/serial/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/idas/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/idas/parallel/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/idas/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/idas/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/idas/parallel/Makefile"
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
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/fcmix_serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/fcmix_serial/Makefile"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/kinsol/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/kinsol/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/kinsol/serial/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/parallel/Makefile"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/fcmix_parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/kinsol/fcmix_parallel/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/kinsol/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/kinsol/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} "
  fi

fi

# CPODES module
if test "X${CPODES_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/cpodes"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} src/cpodes/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cpodes/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cpodes/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cpodes/serial/Makefile"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cpodes/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cpodes/serial"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/cpodes/serial/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cpodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cpodes/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/cpodes/parallel/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cpodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cpodes/parallel"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} test_examples/cpodes/parallel/Makefile"
  fi

fi

# Add Fortran update script to the list of additional files to be generated
if test "X${BUILD_F77_UPDATE_SCRIPT}" = "Xyes"; then
  SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} examples/fortran_update.sh:config/fortran_update.in"
fi

# sundialsTB modules
if test "X${STB_ENABLED}" = "Xyes"; then

  if test "X${CVODES_ENABLED}" = "Xyes" || test "X${IDA_ENABLED}" = "Xyes" || test "X${KINSOL_ENABLED}" = "Xyes"; then
    SLV_MODULES="${SLV_MODULES} sundialsTB"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} sundialsTB/Makefile:sundialsTB/Makefile.in"
    SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} sundialsTB/startup_STB.m:sundialsTB/startup_STB.in"
  fi

  if test "X${CVODES_ENABLED}" = "Xyes"; then
    STB_MODULES="${STB_MODULES} cvodes/cvm/src"
    SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} sundialsTB/cvodes/cvm/src/setup.m:sundialsTB/cvodes/cvm/src/setup.m.in"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} sundialsTB/cvodes/cvm/src/Makefile:sundialsTB/cvodes/cvm/src/Makefile.in"
  fi

  if test "X${IDA_ENABLED}" = "Xyes"; then
    STB_MODULES="${STB_MODULES} idas/idm/src"
    SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} sundialsTB/idas/idm/src/setup.m:sundialsTB/idas/idm/src/setup.m.in"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} sundialsTB/idas/idm/src/Makefile:sundialsTB/idas/idm/src/Makefile.in"
  fi

  if test "X${KINSOL_ENABLED}" = "Xyes"; then
    STB_MODULES="${STB_MODULES} kinsol/kim/src"
    SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} sundialsTB/kinsol/kim/src/setup.m:sundialsTB/kinsol/kim/src/setup.m.in"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} sundialsTB/kinsol/kim/src/Makefile:sundialsTB/kinsol/kim/src/Makefile.in"
  fi

fi

]) dnl END SUNDIALS_BUILD_MODULES_LIST

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
  C Preporcessor Flags:      ${CPPFLAGS}
  C Compiler:	             ${CC}
  C Compiler Flags           ${CFLAGS}
  C Linker:                  ${CC}
  Linker Flags:              ${LDFLAGS}
  Libraries:                 ${LIBS}"

if test "X${F77_EXAMPLES_ENABLED}" = "Xyes"; then
echo "
  Fortran Compiler:          ${F77}
  Fortran Compiler Flags:    ${FFLAGS}
  Fortran Linker:            ${F77_LNKR_OUT}
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
  MPI-Fortran Linker:        ${MPIF77_LNKR_OUT}"
fi

if test "X${STB_ENABLED}" = "Xyes"; then
echo "
  Matlab executable:         ${MATLAB_CMD}
  Extension of MEX files:    ${MEXEXT}
  Mex options file:          ${MEXOPTS}
  Mex compiler options:      ${MEXFLAGS}
  Mex linker flags:          ${MEXLDADD}"
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
if test "X${STB_ENABLED}" = "Xyes"; then
echo "  sundialsTB installed in:   ${STB_PATH}"
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

if test "X${STB_ENABLED}" = "Xyes"; then
  THIS_LINE="sundialsTB:"
  if test "X${CVODES_ENABLED}" = "Xyes" && test -d ${srcdir}/sundialsTB/cvodes ; then
    THIS_LINE="${THIS_LINE} cvodes"
  fi
  if test "X${IDA_ENABLED}" = "Xyes" && test -d ${srcdir}/sundialsTB/idas ; then
    THIS_LINE="${THIS_LINE} idas"
  fi
  if test "X${KINSOL_ENABLED}" = "Xyes"; then
     THIS_LINE="${THIS_LINE} kinsol"
  fi
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

if test "X${DEV_EXAMPLES_ENABLED}" = "Xyes"; then
echo "  Serial C examples (dev):   ${SERIAL_DEV_C_EXAMPLES}"
echo "  Parallel C examples (dev): ${PARALLEL_DEV_C_EXAMPLES}"
fi
fi


echo "  
  Type 'make' and then 'make install' to build and install ${PACKAGE_STRING}."



echo "
----------------------------------
Finished SUNDIALS Configure Script
----------------------------------
"

]) dnl END SUNDIALS_REPORT
