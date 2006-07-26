# -----------------------------------------------------------------
# $Revision: 1.41 $
# $Date: 2006-07-26 22:58:54 $
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

# Set CC to default value
if test "X${CC}" = "X"; then
  CC="cc"
fi

# Unfortunately, autoconf gets confused if F77 is undefined (it tries
# to test the compiler before it actually searches for the compiler)
if test "X${F77}" = "X"; then
  F77="f77"
fi

# Allow user to override default flags if corresponding environment
# variable is already defined
if test "X${CFLAGS}" = "X"; then
  CFLAGS_USER_OVERRIDE="no"
else
  CFLAGS_USER_OVERRIDE="yes"
fi
if test "X${CXXFLAGS}" = "X"; then
  CXXFLAGS_USER_OVERRIDE="no"
else
  CXXFLAGS_USER_OVERRIDE="yes"
fi
if test "X${FFLAGS}" = "X"; then
  FFLAGS_USER_OVERRIDE="no"
else
  FFLAGS_USER_OVERRIDE="yes"
fi
if test "X${FLIBS}" = "X"; then
  FLIBS_USER_OVERRIDE="no"
else
  FLIBS_USER_OVERRIDE="yes"
fi

# Set defaults for config/sundials_config.in file
F77_MANGLE_MACRO1=""
F77_MANGLE_MACRO2=""
F77_CASE=""
F77_UNDERSCORES=""
PRECISION_LEVEL=""
GENERIC_MATH_LIB=""
F77_MPI_COMM_F2C=""

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

# Check if user wants to disable Fortran support (FCMIX components)
F77_ENABLED="yes"
AC_ARG_ENABLE([f77],
[AC_HELP_STRING([--disable-f77], [disable Fortran support])],
[
if test "X${enableval}" = "Xno"; then
  F77_ENABLED="no"
  FCMIX_ENABLED="no"
fi
],
[
if test "X${CVODE_ENABLED}" = "Xno" && test "X${KINSOL_ENABLED}" = "Xno" && test "X${IDA_ENABLED}" = "Xno"; then
  F77_ENABLED="no"
  FCMIX_ENABLED="no"
else
  F77_ENABLED="yes"
  FCMIX_ENABLED="yes"
fi
])

# Check if user wants to disable support for MPI
# If not, then make certain source directory actually exists
MPI_ENABLED="yes"
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
EXAMPLES_ENABLED="no"
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
STB_ENABLED="no"
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

# Initialize variables that may NOT necessarily be initialized
# during normal execution
# Should NOT use uninitialized variables
CXX_ENABLED="no"
IDAS_ENABLED="no"
CXX_OK="no"
F77_OK="no"
MPI_CXX_COMP_OK="no"
MPI_F77_COMP_OK="no"
STB_PARALLEL_OK="no"

]) dnl END SUNDIALS_ENABLES

#------------------------------------------------------------------
# TEST DEVELOPMENT ENABLES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEV_ENABLES],
[

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
        
]) dnl END SUNDIALS_DEV_ENABLES

#------------------------------------------------------------------
# CHECK C COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_CC],
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

# Determine absolute pathname for specified C compiler
CC_TEMP1="${CC}"
CC_TEMP2=`basename "${CC}"`
CC_EXEC="${CC_TEMP2}"
# If only the executable name was given, then determine absolute pathname
if test "X${CC_TEMP1}" = "X${CC_TEMP2}"; then
  AC_PATH_PROG([CC_COMP],[${CC}],[none])
  # We already know the C compiler exists, so this test is provided only
  # for the sake of completeness
  if test "X${CC_COMP}" = "Xnone"; then
    CC="${CC_TEMP1}"
  else
    CC="${CC_COMP}"
  fi
fi

AC_MSG_CHECKING([for extra C compiler flags])
AC_ARG_WITH(cflags,
[AC_HELP_STRING([--with-cflags=ARG],[add extra C compiler flags])],
[
AC_MSG_RESULT([${withval}])
USER_CFLAGS="${withval}"
],
[
AC_MSG_RESULT([none])
USER_CFLAGS=""
])

# Set CPP to command that runs C preprocessor
# Defines CPP
AC_PROG_CPP

AC_MSG_CHECKING([for extra C/C++ preprocessor flags])
AC_ARG_WITH(cppflags,
[AC_HELP_STRING([--with-cppflags=ARG],[add extra C/C++ preprocessor flags])],
[
AC_MSG_RESULT([${withval}])
if test "X${CPPFLAGS}" = "X"; then
  CPPFLAGS="${withval}"
else
  CPPFLAGS="${CPPFLAGS} ${withval}"
fi
],
[
AC_MSG_RESULT([none])
])

# Add default flag(s) to CFLAGS only if not already defined
if test "X${CFLAGS_USER_OVERRIDE}" = "Xno"; then
  SUNDIALS_DEFAULT_CFLAGS
fi

if test "X${USER_CFLAGS}" = "X"; then
  CFLAGS="${CFLAGS}"
else
  CFLAGS="${CFLAGS} ${USER_CFLAGS}"
fi

AC_MSG_CHECKING([for extra linker flags])
AC_ARG_WITH(ldflags,
[AC_HELP_STRING([--with-ldflags=ARG],[add extra linker flags])],
[
AC_MSG_RESULT([${withval}])
if test "X${LDFLAGS}" = "X"; then
  LDFLAGS="${withval}"
else
  LDFLAGS="${LDFLAGS} ${withval}"
fi
],
[
AC_MSG_RESULT([none])
])

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
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_CFLAGS],
[

# Note: Although NOT "used", we will keep this particular test

# Let user-specified CFLAGS variable pass through
case $host in 

  # IA-32 system running Linux
  i*-pc-linux-*)

    if test "X${GCC}" = "Xyes"; then
      if test "X${CFLAGS}" = "X"; then
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CFLAGS=""
        fi
      else
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CFLAGS="${CFLAGS}"
        fi
      fi
    fi

  ;;

esac

]) dnl END SUNDIALS_DEFAULT_CFLAGS

#------------------------------------------------------------------
# CHECK FORTRAN COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_F77],
[

AC_LANG_PUSH([Fortran 77])

F77_OK="yes"

# Note: This check may no longer be needed
if test "X${F77}" = "X"; then

  # If F77="" then issue warning (means did NOT find a valid Fortran compiler)
  # Do NOT abort since Fortran compiler is NOT required
  AC_MSG_WARN([cannot find Fortran compiler])
  echo ""
  echo "   Unable to find a functional Fortran compiler."
  echo ""
  echo "   Disabling Fortran support..."
  echo ""
  F77_OK="no"
  SUNDIALS_WARN_FLAG="yes"

else

  # Provide variable description templates for config.hin and config.h files
  # Required by autoheader utility
  AH_TEMPLATE([SUNDIALS_UNDERSCORE_NONE],
              [FCMIX: Do NOT append any underscores to functions names])
  AH_TEMPLATE([SUNDIALS_UNDERSCORE_ONE],
              [FCMIX: Append ONE underscore to function names])
  AH_TEMPLATE([SUNDIALS_UNDERSCORE_TWO],
              [FCMIX: Append TWO underscores to function names])

  # Provided in case SUNDIALS_F77_WRAPPERS cannot determine name-mangling scheme
  AC_ARG_WITH(f77underscore, 
  [AC_HELP_STRING([--with-f77underscore=ARG],[specify number of underscores to append to function names (none/one/two) [AUTO]],[])],
  [
  UNDERSCORE_ARG_GIVEN="yes"
  if test "X${withval}" = "Xnone"; then
    AC_DEFINE([SUNDIALS_UNDERSCORE_NONE],[1],[])
    F77_UNDERSCORES="#define SUNDIALS_UNDERSCORE_NONE 1"
  elif test "X${withval}" = "Xone"; then
    AC_DEFINE([SUNDIALS_UNDERSCORE_ONE],[1],[])
    F77_UNDERSCORES="#define SUNDIALS_UNDERSCORE_ONE 1"
  elif test "X${withval}" = "Xtwo"; then
    AC_DEFINE([SUNDIALS_UNDERSCORE_TWO],[1],[])
    F77_UNDERSCORES="#define SUNDIALS_UNDERSCORE_TWO 1"
  else
    AC_MSG_ERROR([invalid input])
  fi
  ],
  [
  UNDERSCORE_ARG_GIVEN="no"
  ])

  # Provide variable description templates for config.hin and config.h files
  # Required by autoheader utility
  AH_TEMPLATE([SUNDIALS_CASE_UPPER],
              [FCMIX: Make function names uppercase])
  AH_TEMPLATE([SUNDIALS_CASE_LOWER],
              [FCMIX: Make function names lowercase])

  # Provided in case SUNDIALS_F77_WRAPPERS cannot determine name-mangling scheme
  AC_ARG_WITH(f77case, 
  [AC_HELP_STRING([--with-f77case=ARG   ],[specify case of function names (lower/upper) [AUTO]],[])],
  [
  CASE_ARG_GIVEN="yes"
  if test "X${withval}" = "Xupper"; then
    AC_DEFINE([SUNDIALS_CASE_UPPER],[1],[])
    F77_CASE="#define SUNDIALS_CASE_UPPER 1"
  elif test "X${withval}" = "Xlower"; then
    AC_DEFINE([SUNDIALS_CASE_LOWER],[1],[])
    F77_CASE="#define SUNDIALS_CASE_LOWER 1"
  else
    AC_MSG_ERROR([invalid input])
  fi
  ],
  [
  CASE_ARG_GIVEN="no"
  ])

  # Check if user used "--with-f77underscore" or "--with-f77case"
  # If both flags were used, then just continue
  if test "X${UNDERSCORE_ARG_GIVEN}" = "Xyes" || test "X${CASE_ARG_GIVEN}" = "Xyes"; then
    # Only used "--with-f77underscore" flag, so set default case
    if test "X${UNDERSCORE_ARG_GIVEN}" = "Xyes" && test "X${CASE_ARG_GIVEN}" = "Xno"; then
      AC_DEFINE([SUNDIALS_CASE_LOWER],[1],[])
      F77_CASE="#define SUNDIALS_CASE_LOWER 1"
    # Only used "--with-f77case" flag, so set default number of underscores
    elif test "X${UNDERSCORE_ARG_GIVEN}" = "Xno" && test "X${CASE_ARG_GIVEN}" = "Xyes"; then
      AC_DEFINE([SUNDIALS_UNDERSCORE_ONE],[1],[])
      F77_UNDERSCORES="#define SUNDIALS_UNDERSCORE_ONE 1"
    fi
    RUN_F77_WRAPPERS="no"
  # Only call SUNDIALS_F77_WRAPPERS if user did NOT use "--with-f77underscore" or "--with-f77case" flags
  else
    RUN_F77_WRAPPERS="yes"
  fi

  # Determine absolute pathname for specified Fortran compiler
  F77_TEMP1="${F77}"
  F77_TEMP2=`basename "${F77}"`
  # SUNDIALS_DEFAULT_FFLAGS needs just the executable name
  F77_EXEC="${F77_TEMP2}"
  # If only the executable name was given, then determine absolute pathname
  if test "X${F77_TEMP1}" = "X${F77_TEMP2}"; then
    AC_PATH_PROG([F77_COMP],[${F77}],[none])
    # We already know the Fortran compiler exists, so this test is provided only
    # for the sake of completeness
    if test "X${F77_COMP}" = "Xnone"; then
      F77="${F77_TEMP1}"
    else
      F77="${F77_COMP}"
    fi
  fi

  AC_MSG_CHECKING([for extra Fortran compiler flags])
  AC_ARG_WITH(fflags,
  [AC_HELP_STRING([--with-fflags=ARG],[add extra Fortran compiler flags])],
  [
  AC_MSG_RESULT([${withval}])
  USER_FFLAGS="${withval}"
  ],
  [
  AC_MSG_RESULT([none])
  USER_FFLAGS=""
  ])
  
  # Add default flag(s) to FFLAGS only if not already defined
  if test "X${FFLAGS_USER_OVERRIDE}" = "Xno"; then
    SUNDIALS_DEFAULT_FFLAGS
  fi

  if test "X${USER_FFLAGS}" = "X"; then
    FFLAGS="${FFLAGS}"
  else
    FFLAGS="${FFLAGS} ${USER_FFLAGS}"
  fi

  # Add any required linker flags to FLIBS
  # Note: if FLIBS is defined, it is left unchanged
  AC_F77_LIBRARY_LDFLAGS

fi

# Determine how to properly mangle function names so Fortran subroutines can
# call C functions included in SUNDIALS libraries
# Defines C preprocessor macros F77_FUNC and F77_FUNC_
if test "X${RUN_F77_WRAPPERS}" = "Xyes"; then
  if test "X${F77_OK}" = "Xyes"; then
    SUNDIALS_F77_WRAPPERS
  fi
fi

# Check if we must use a Fortran compiler to link the Fortran examples
# (default is to use CC)
if test "X${F77_OK}" = "Xyes"; then
  SUNDIALS_F77_LNKR_CHECK
fi

AC_LANG_POP([Fortran 77])

]) dnl END SUNDIALS_SET_F77

#------------------------------------------------------------------
# DETERMINE F77 LINKER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_F77_LNKR_CHECK],
[

# Note: SUNDIALS_CPLUSPLUS_CHECK macro was called by SUNDIALS_SET_CC,
# so we need only check the value of USING_CPLUSPLUS_COMP
if test "X${USING_CPLUSPLUS_COMP}" = "Xyes"; then
  RUN_F77_LNKR_CHECK="yes"
else
  RUN_F77_LNKR_CHECK="no"
fi

AC_MSG_CHECKING([which linker to use])
# Perform the next test only if using a C++ compiler to build SUNDIALS
if test "X${RUN_F77_LNKR_CHECK}" = "Xyes"; then

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
AC_MSG_RESULT([${F77_LNKR_OUT}])

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

AC_MSG_CHECKING([Fortran name-mangling scheme])

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
    F77_FUNC_CASE="lowercase"
  else
    F77_FUNC_CASE="uppercase"
  fi
  # Determine number of underscores to append (none, one, two)
  if test "X${j}" = "X"; then
    F77_FUNC_UNDERSCORES="no underscores"
  elif test "X${j}" = "X_"; then
    F77_FUNC_UNDERSCORES="one underscore"
  else
    F77_FUNC_UNDERSCORES="two underscores"
  fi
  AC_MSG_RESULT([${F77_FUNC_CASE} with ${F77_FUNC_UNDERSCORES}])
  # Set exported macro definition (F77_FUNC)
  if test "X${F77_FUNC_CASE}" = "Xlowercase"; then
    if test "X${F77_FUNC_UNDERSCORES}" = "Xno underscores"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[name],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) name"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) name"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xone underscore"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[name ## _],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) name ## _"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) name ## _"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xtwo underscores"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[name ## __],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) name ## __"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) name ## __"
    fi
  elif test "X${F77_FUNC_CASE}" = "Xuppercase"; then
    if test "X${F77_FUNC_UNDERSCORES}" = "Xno underscores"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[NAME],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) NAME"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) NAME"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xone underscore"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[NAME ## _],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) NAME ## _"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) NAME ## _"
    elif test "X${F77_FUNC_UNDERSCORES}" = "Xtwo underscores"; then
      AC_DEFINE(F77[_FUNC(name,NAME)],[NAME ## __],[])
      F77_MANGLE_MACRO1="#define F77_FUNC(name,NAME) NAME ## __"
      F77_MANGLE_MACRO2="#define F77_FUNC_(name,NAME) NAME ## __"
    fi
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

])

# Reset language (remove 'Fortran 77' from stack)
AC_LANG_POP([Fortran 77])

# Remove temporary file
rm -f f77_wrapper_check.${ac_objext}

]) dnl END SUNDIALS_SET_F77

#------------------------------------------------------------------
# DEFAULT FFLAGS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_FFLAGS],
[

# FIXME: Should IRIX and Tru64 options overwrite FFLAGS?
case $host in

  # IA-32 system running Linux
  # Let user-specified FFLAGS variable pass through
  i*-pc-linux-*)

    if test "X${G77}" = "Xyes"; then
      if test "X${FFLAGS}" = "X"; then
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          FFLAGS=""
        fi
      else
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          FFLAGS="${FFLAGS}"
        fi
      fi
    fi

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

esac

]) dnl END SUNDIALS_DEFAULT_FFLAGS

#------------------------------------------------------------------
# CHECK C++ COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_CXX],
[

AC_LANG_PUSH([C++])

CXX_OK="yes"

# CXX and CCC are common so check both
if test "X${CXX}" = "X"; then
  if test "X${CCC}" = "X"; then
    CXX=""
  else
    CXX="${CCC}"
  fi
fi

# Sets GXX="yes" if CXX="g++" (GNU C++ compiler)
# Search for C++ compiler given by user first
AC_PROG_CXX(${CXX} g++ CC gcc c++ cxx)

# If CXX="" then issue warning (means did NOT find a valid C++ compiler)
# Do NOT abort since C++ compiler is NOT required
if test "X${CXX}" = "X"; then

  AC_MSG_WARN([cannot find C++ compiler])
  echo ""
  echo "   Unable to find a functional C++ compiler."
  echo ""
  echo "   Disabling dependent modules..."
  echo ""
  CXX_OK="no"
  SUNDIALS_WARN_FLAG="yes"

else

  # Determine absolute pathname for specified C++ compiler
  CXX_TEMP1="${CXX}"
  CXX_TEMP2=`basename "${CXX}"`
  CXX_EXEC="${CXX_TEMP2}"
  # If only the executable name was given, then determine absolute pathname
  if test "X${CXX_TEMP1}" = "X${CXX_TEMP2}"; then
    AC_PATH_PROG([CXX_COMP],[${CXX}],[none])
    # We already know the C++ compiler exists, so this test is provided only
    # for the sake of completeness
    if test "X${CXX_COMP}" = "Xnone"; then
      CXX="${CXX_TEMP1}"
    else
      CXX="${CXX_COMP}"
    fi
  fi

  AC_MSG_CHECKING([for extra C++ compiler flags])
  AC_ARG_WITH(cxxflags,
  [AC_HELP_STRING([--with-cxxflags=ARG],[add extra C++ compiler flags])],
  [
  AC_MSG_RESULT([${withval}])
  USER_CXXFLAGS="${withval}"
  ],
  [
  AC_MSG_RESULT([none])
  USER_CXXFLAGS=""
  ])

  # Set CXXCPP to command that runs C++ preprocessor
  AC_PROG_CXXCPP

  # Add default flag(s) to CXXFLAGS only if not already defined
  if test "X${CXXFLAGS_USER_OVERRIDE}" = "Xno"; then
    SUNDIALS_DEFAULT_CXXFLAGS
  fi

  if test "X${USER_CXXFLAGS}" = "X"; then
    CXXFLAGS="${CXXFLAGS}"
  else
    CXXFLAGS="${CXXFLAGS} ${USER_CXXFLAGS}"
  fi

  # Check for complex.h header file (required)
  AC_CHECK_HEADER([complex.h])

fi

AC_LANG_POP([C++])

]) dnl END SUNDIALS_SET_CXX

#------------------------------------------------------------------
# DEFAULT CXXFLAGS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_CXXFLAGS],
[

# Note: Although NOT "used", we will keep this particular test

# Let user-specified CXXFLAGS variable pass through
case $host in 

  # IA-32 system running Linux
  i*-pc-linux-*)

    if test "X${GXX}" = "Xyes"; then
      if test "X${CXXFLAGS}" = "X"; then
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CXXFLAGS=""
        fi
      else
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CXXFLAGS="${CXXFLAGS}"
        fi
      fi
    fi

  ;;

esac

]) dnl END SUNDIALS_DEFAULT_CXXFLAGS

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

# Do NOT even check for MPI support for Fortran if serial Fortran compiler does NOT work
# Need serial Fortran compiler to determine name-mangline scheme, and FCMIX libraries will
# only be built if a serial Fortran compiler is found
if test "X${F77_OK}" = "Xyes"; then
  # Check MPI-Fortran compiler (either MPI compiler script or regular Fortran compiler)
  if test "X${USE_MPIF77_SCRIPT}" = "Xyes"; then
    SUNDIALS_CHECK_MPIF77
  else
    MPIF77_COMP="${F77}"
    MPIF77="${F77}"
    SUNDIALS_CHECK_F77_WITH_MPI
  fi

else

  AC_MSG_WARN([serial F77 does not work so we do not even bother with MPI-F77])

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

# Only continue if actually found a valid Fortran compiler
if test "X${F77_OK}" = "Xyes"; then

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

else

  AC_MSG_WARN([cannot find Fortran compiler])
  MPI_F77_COMP_OK="no"
  SUNDIALS_WARN_FLAG="yes"

fi

# Restore FFLAGS, LDFLAGS and LIBS
FFLAGS="${SAVED_FFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

AC_LANG_POP([Fortran 77])

]) dnl END SUNDIALS_CHECK_F77_WITH_MPI

#------------------------------------------------------------------
# CHECK MPI-C++ COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_MPICXX],
[

AC_MSG_CHECKING([if using MPI-C++ script])
AC_ARG_WITH(mpicxx,
[AC_HELP_STRING([--with-mpicxx[[[[=ARG]]]]],[specify MPI-C++ compiler to use @<:@mpiCC@:>@],
                [                                ])],
[
if test "X${withval}" = "Xno"; then
  USE_MPICXX_SCRIPT="no"
else
  USE_MPICXX_SCRIPT="yes"
  MPICXX_COMP="${withval}"
fi
],
[
USE_MPICXX_SCRIPT="yes"
MPICXX_COMP="mpiCC"
])
AC_MSG_RESULT([${USE_MPICXX_SCRIPT}])

# Check MPI-C++ compiler (either MPI compiler script or regular C++ compiler)
if test "X${USE_MPICXX_SCRIPT}" = "Xyes"; then
  SUNDIALS_CHECK_MPICXX
else
  MPICXX_COMP="${CXX}"
  MPICXX="${CXX}"
  SUNDIALS_CHECK_CXX_WITH_MPI
fi

AC_ARG_WITH([],[     ],[])

]) dnl END SUNDIALS_SET_MPICXX

#------------------------------------------------------------------
# TEST MPI-C++ COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_MPICXX],
[

# Test the MPI-C++ compiler (meaning test MPICXX_COMP)
# Check if MPI-C++ compiler can be found

AC_MSG_CHECKING([if absolute path to ${MPICXX_COMP} was given])

# CASE 1: MPICXX_COMP was found (cannot check if executable because the
# "-x" flag is NOT portable)
if test -f ${MPICXX_COMP} ; then

  AC_MSG_RESULT([yes])
  MPICXX_COMP_EXISTS="yes"
  # Determine MPI_INC_DIR and MPI_LIB_DIR for use by Makefile
  MPI_BASE_DIR=`AS_DIRNAME(["${MPICXX_COMP}"])`

# CASE 2: MPICXX_COMP could NOT be found and MPI_ROOT_DIR was NOT specified,
# so search in PATH
else

  AC_MSG_RESULT([no])

  if test "X${MPI_ROOT_DIR}" = "X"; then

    # Try to find location of executable (perhaps directory was entered
    # incorrectly)
    TEMP_MPICXX_COMP=`basename "${MPICXX_COMP}"`
    AC_PATH_PROG([MPICXX_COMP],[${TEMP_MPICXX_COMP}],[none])
    # Cannot find executable in PATH
    if test "X$S{MPICXX_COMP}" = "Xnone"; then
      MPICXX_COMP_EXISTS="no"
      MPICXX_COMP=""
    # Found executable and set MPICXX_COMP to absolute pathname
    else
      MPICXX_COMP_EXISTS="yes"
      MPI_BASE_DIR=`AS_DIRNAME(["${MPICXX_COMP}"])`
    fi

  # CASE 3: MPICXX_COMP could NOT be found, but MPI_ROOT_DIR was specified
  else

    AC_MSG_CHECKING([if ${MPICXX_COMP} exists in ${MPI_ROOT_DIR}/bin])
    # MPICXX_COMP should really only contain an executable name
    # Found location of MPICXX_COMP
    if test -f ${MPI_ROOT_DIR}/bin/${MPICXX_COMP} ; then
      AC_MSG_RESULT([yes])
      MPICXX_COMP_EXISTS="yes"
      MPICXX_COMP="${MPI_ROOT_DIR}/bin/${MPICXX_COMP}"
    # Could NOT find MPICXX_COMP anywhere
    else
      AC_MSG_RESULT([no])
      MPICXX_COMP_EXISTS="no"
      MPICXX_COMP=""
    fi

  fi

fi

# Issue warning message if MPICXX_COMP does NOT exist, else set MPICXX
if test "X${MPICXX_COMP_EXISTS}" = "Xyes"; then

  MPICXX="${MPICXX_COMP}"
  MPI_CXX_COMP_OK="yes"

  # Note that we do not have to worry about empty MPI_INC_DIR and MPI_LIB_DIR
  # here as they were set in SUNDIALS_CHECK_MPICC

else

  AC_MSG_WARN([cannot find MPI-C++ compiler])
  echo ""
  echo "   Unable to find a functional MPI-C++ compiler."
  echo ""
  echo "   Try using --with-mpicxx to specify a MPI-C++ compiler script,"
  echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
  echo "   to specify the locations of all relevant MPI files, or"
  echo "   --with-mpi-root to specify the base installation directory"
  echo "   of the MPI implementation to be used."
  echo ""
  echo "   Disabling dependent modules...."
  echo ""
  MPICXX=""
  MPI_CXX_COMP_OK="no"
  SUNDIALS_WARN_FLAG="yes"

fi

]) dnl END SUNDIALS_CHECK_MPICXX

#------------------------------------------------------------------
# TEST C++ COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_CXX_WITH_MPI],
[

# Test if we can compile MPI programs using the C++ compiler
# and current MPI settings

AC_MSG_NOTICE([Testing C++ with MPI settings])

AC_LANG_PUSH([C++])

# Save copies of CPPFLAGS, LDFLAGS and LIBS
# Temporarily overwritten so we can test MPI implementation
SAVED_CPPFLAGS="${CPPFLAGS}"
SAVED_LDFLAGS="${LDFLAGS}"
SAVED_LIBS="${LIBS}"

# Only continue if actually found a valid C++ compiler
if test "X${CXX_OK}" = "Xyes"; then

  # This may seem redundant, but we are not guaranteed that
  # either SUNDIALS_CHECK_CC_WITH_MPI or SUNDIALS_CHECK_F77_WITH_MPI
  # has been executed
  # Determine location of MPI header files (find MPI include directory)
  MPI_EXISTS="yes"

  AC_MSG_CHECKING([for location of MPI implementation])

  # MPI include directory was NOT explicitly specified so check if MPI root
  # directory was given by user
  if test "X${MPI_INC_DIR}" = "X"; then
    # MPI root directory was NOT given so issue a warning message
    if test "X${MPI_ROOT_DIR}" = "X"; then
      AC_MSG_RESULT([not found])
      MPI_EXISTS="no"
      AC_MSG_WARN([cannot find MPI implementation files])
      echo ""
      echo "   Unable to find MPI implementation files."
      echo ""
      echo "   Try using --with-mpicxx to specify a MPI-C++ compiler script,"
      echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
      echo "   to specify the locations of all relevant MPI files, or"
      echo "   --with-mpi-root to specify the base installation directory"
      echo "   of the MPI implementation to be used."
      echo ""
      echo "   Disabling dependent modules..."
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
        CPPLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
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
      echo "   Unable to compile MPI program using C++ compiler because"
      echo "   MPI libraries were not specified."
      echo ""
      echo "   Try using --with-mpi-libdir and --with-mpi-libs to"
      echo "   specify the location and names of the MPI libraries."
      echo ""
      echo "   Disabling dependent modules..."
      echo ""
      MPI_CXX_COMP_OK="no"
    # MPI libraries were specified so update LIBS
    else
      AC_MSG_RESULT([${MPI_LIBS}])
      if test "X${LIBS}" = "X"; then
        LIBS="${MPI_LIBS}"
      else
        LIBS="${LIBS} ${MPI_LIBS}"
      fi
      # Set the MPI_CXX_COMP_OK variable to NULL so we can conditionally execute
      # the next test
      MPI_CXX_COMP_OK=""
    fi

    if test "X${MPI_CXX_COMP_OK}" = "X"; then
      AC_MSG_CHECKING([if C++ compiler can compile MPI programs])
      AC_LINK_IFELSE(
      [AC_LANG_PROGRAM([[#include "mpi.h"]],[[int c; car **v; MPI_Init(&c,&v);]])],
      [AC_MSG_RESULT([yes])
       MPI_CXX_COMP_OK="yes"],
      [AC_MSG_RESULT([no])
       AC_MSG_WARN([C++ compiler cannot compile MPI programs])
       echo ""
       echo "   Unable to compile MPI program using C++ compiler."
       echo ""
       echo "   Try using --with-mpicxx to specify a MPI-C++ compiler script,"
       echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
       echo "   to specify the locations of all relevant MPI files, or"
       echo "   --with-mpi-root to specify the base installation directory"
       echo "   of the MPI implementation to be used."
       echo ""
       echo "   Disabling dependent modules..."
       echo ""
       MPI_CXX_COMP_OK="no"
       SUNDIALS_WARN_FLAG="yes"])
    fi

  else
    MPI_CXX_COMP_OK="no"
  fi

else

  AC_MSG_WARN([cannot find C++ compiler])
  MPI_CXX_COMP_OK="no"
  SUNDIALS_WARN_FLAG="yes"

fi

# Restore CPPFLAGS, LDFLAGS and LIBS
CPPFLAGS="${SAVED_CPPFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

AC_LANG_POP([C++])

]) dnl END SUNDIALS_CHECK_CXX_WITH_MPI

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
  if test "X${F77_ENABLED}" = "Xyes"; then
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

AC_ARG_WITH([],[      ],[])

# Find Matlab program. If --with-matlab=<matlab> was passed to the configure 
# script, then we only check if it exists and is an executable file. We'll try
# to run it later... Otherwise, use AC_PATH_PROG to search for matlab.
MATLAB_CMD="none"
AC_ARG_WITH([matlab], 
[AC_HELP_STRING([--with-matlab=MATLAB], [specify Matlab executable])],
[
if test -x ${withval} ; then
  MATLAB_CMD="${withval}"
else
  AC_MSG_WARN([invalid value '${withval}' for --with-matlab])
  STB_ENABLED="no"
fi
],
[
AC_PATH_PROG(MATLAB_CMD, matlab, "none")
if test "X${MATLAB_CMD}" = "Xnone"; then
  STB_ENABLED="no"
fi
])

# Under cygwin, standardize MATLAB_CMD
if test "X${STB_ENABLED}" = "Xyes"; then
  case $build_os in
    *cygwin*)
      MATLAB_CMD=`cygpath -a -m "${MATLAB_CMD}"`
      ;;
  esac
fi

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

# Decide whether to enable parallel support in sundialsTB
if test "X${STB_ENABLED}" = "Xyes"; then
  if test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_C_COMP_OK}" = "Xyes"; then
    STB_PARALLEL="yes"
  else
    STB_PARALLEL="no"
  fi
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
  ( eval ${MATLAB_CMD} -nojvm -nosplash -r cvmextest_script ) 2>/dev/null 1>&2

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
    STB_INSTDIR="${MATLAB}/toolbox/sundialsTB"
  ])
  AC_MSG_RESULT([${STB_INSTDIR}])

fi

# Display message if Matlab support had to be disabled
if test "X${STB_ENABLED}" = "Xno"; then

  AC_MSG_WARN([Matlab support was disabled])
  echo ""
  echo "   Matlab configuration for sundialsTB was disabled."
  echo ""

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
# ENABLE EXAMPLES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_ENABLE_EXAMPLES],
[

# Check libtool settings to determine which compiler and linker commands
# should be used
# Must export LIBTOOL_CMD, COMPILER_PREFIX and LINKER_PREFIX via AC_SUBST
# If building shared libraries, then use libtool
if test "X${enable_shared}" = "Xyes"; then

  LIBTOOL_CMD="LIBTOOL = ${LIBTOOL}"
  COMPILER_PREFIX="\$(LIBTOOL) --mode=compile"
  LINKER_PREFIX="\$(LIBTOOL) --mode=link"
  OBJ_EXT="lo"

# If building static libraries, then use regular C compiler
else

  LIBTOOL_CMD=""
  COMPILER_PREFIX=""
  LINKER_PREFIX=""
  OBJ_EXT="o"

fi

# Check if we need to build the Fortran examples
# Recall F77_ENABLED is set based upon CVODE_ENABLED and KINSOL_ENABLED,
# so the check here is rather simple:
# F77_EXAMPLES = F77_ENABLED
F77_EXAMPLES="no";
if test "X${F77_ENABLED}" = "Xyes"; then
  F77_EXAMPLES="yes"
fi

# Check if we need to build the C++ examples
CXX_EXAMPLES="no"
if test "X${CXX_ENABLED}" = "Xyes"; then
  if test "X${CVODES_ENABLED}" = "Xyes" || test "X${IDAS_ENABLED}" = "Xyes"; then
    CXX_EXAMPLES="yes"
  fi
fi

# Check if examples can actually be built
SERIAL_C_EXAMPLES="no"
SERIAL_CXX_EXAMPLES="no"
SERIAL_F77_EXAMPLES="no"

PARALLEL_C_EXAMPLES="no"
PARALLEL_CXX_EXAMPLES="no"
PARALLEL_F77_EXAMPLES="no"

# Check C examples

AC_MSG_CHECKING([if we can build serial C examples])
SERIAL_C_EXAMPLES="yes"
AC_MSG_RESULT([${SERIAL_C_EXAMPLES}])

if test "X${MPI_ENABLED}" = "Xyes"; then

  AC_MSG_CHECKING([if we can build parallel C examples])
  if test "X${MPI_C_COMP_OK}" = "Xyes"; then
    PARALLEL_C_EXAMPLES="yes"
  fi
  AC_MSG_RESULT([${PARALLEL_C_EXAMPLES}])

fi

# Check Fortran examples
if test "X${F77_EXAMPLES}" = "Xyes"; then

  AC_MSG_CHECKING([if we can build serial Fortran examples])
  if test "X${F77_OK}" = "Xyes"; then
    SERIAL_F77_EXAMPLES="yes"
  fi
  AC_MSG_RESULT([${SERIAL_F77_EXAMPLES}])

  if test "X${MPI_ENABLED}" = "Xyes"; then

    AC_MSG_CHECKING([if we can build parallel Fortran examples])
    if test "X${MPI_F77_COMP_OK}" = "Xyes"; then
      PARALLEL_F77_EXAMPLES="yes"
    fi
    AC_MSG_RESULT([${PARALLEL_F77_EXAMPLES}])

  fi

fi

# Check if the Fortran update script (config/fortran_update.in) should
# be initialized
BUILD_F77_UPDATE_SCRIPT="no"
if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" || test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
  BUILD_F77_UPDATE_SCRIPT="yes";
fi

# Check C++ examples
if test "X${CXX_EXAMPLES}" = "Xyes"; then

  AC_MSG_CHECKING([if we can build serial C++ examples])
  if test "X${CXX_OK}" = "Xyes"; then
    SERIAL_CXX_EXAMPLES="yes"
  fi
  AC_MSG_RESULT([${SERIAL_CXX_EXAMPLES}])

  if test "X${MPI_ENABLED}" = "Xyes"; then

    AC_MSG_CHECKING([if we can build parallel C++ examples])
    if test "X${MPI_CXX_COMP_OK}" = "Xyes"; then
      PARALLEL_CXX_EXAMPLES="yes"
    fi
    AC_MSG_RESULT([${PARALLEL_CXX_EXAMPLES}])

  fi

fi

# Where should we install the examples?

AC_ARG_WITH([],[   ],[])

AC_MSG_CHECKING([where to install the SUNDIALS examples])
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


# Disable developer examples be default so output logic still works
SERIAL_DEV_C_EXAMPLES="no"
PARALLEL_DEV_C_EXAMPLES="no"

RAN_ENABLE_DEV_EXAMPLES="no"

]) dnl END SUNDIALS_ENABLE_EXAMPLES

#------------------------------------------------------------------
# ENABLE DEVELOPER EXAMPLES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_ENABLE_DEV_EXAMPLES],
[

RAN_ENABLE_DEV_EXAMPLES="yes"

# Check if developer examples can actually be built
SERIAL_DEV_C_EXAMPLES="no"
PARALLEL_DEV_C_EXAMPLES="no"

# Check developer examples

AC_MSG_CHECKING([if we can build serial C developer examples])
SERIAL_DEV_C_EXAMPLES="yes"
AC_MSG_RESULT([${SERIAL_DEV_C_EXAMPLES}])

if test "X${MPI_ENABLED}" = "Xyes"; then
  AC_MSG_CHECKING([if we can build parallel C developer examples])
  if test "X${MPI_C_COMP_OK}" = "Xyes"; then
    PARALLEL_DEV_C_EXAMPLES="yes"
  fi
  AC_MSG_RESULT([${PARALLEL_DEV_C_EXAMPLES}])
fi

]) dnl END SUNDIALS_ENABLE_DEV_EXAMPLES

#------------------------------------------------------------------
# BUILD MODULES LIST
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_BUILD_MODULES_LIST],
[

# Initialize the list of Makefiles to be created
SUNDIALS_MAKEFILES="Makefile"

# Initialize list of additional configure files to be created
SUNDIALS_CONFIGFILES="src/sundials/sundials_config.h:src/sundials/sundials_config.in"
SUNDIALS_CONFIGFILES="${SUNDIALS_CONFIGFILES} sundials-config:sundials-config.in"

# Initialize lists of solver modules, exampel modules, and sundialsTB modules
SLV_MODULES="src/sundials"
EXS_MODULES=""
STB_MODULES=""

# NVECTOR modules
if test -d ${srcdir}/src/nvec_ser ; then
  SLV_MODULES="${SLV_MODULES} src/nvec_ser"
fi

if test -d ${srcdir}/src/nvec_par && test "X${MPI_C_COMP_OK}" = "Xyes"; then
  SLV_MODULES="${SLV_MODULES} src/nvec_par"
fi

if test -d ${srcdir}/src/nvec_spcpar && test "X${MPI_C_COMP_OK}" = "Xyes"; then
  SLV_MODULES="${SLV_MODULES} src/nvec_spcpar"
fi

# CVODE module
if test "X${CVODE_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/cvode"

  if test "X${FCMIX_ENABLED}" = "Xyes" && test -d ${srcdir}/src/cvode/fcmix ; then
    SLV_MODULES="${SLV_MODULES} src/cvode/fcmix"
  fi

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/serial"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/fcmix_serial"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvode/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvode/serial"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/parallel"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvode/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvode/fcmix_parallel"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvode/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvode/parallel"
  fi

fi

# CVODES module
if test "X${CVODES_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/cvodes"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/cvodes/serial"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvodes/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvodes/serial"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/cvodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/cvodes/parallel"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/cvodes/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/cvodes/parallel"
  fi

fi

# IDA module
if test "X${IDA_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/ida"

  if test "X${FCMIX_ENABLED}" = "Xyes" && test -d ${srcdir}/src/ida/fcmix ; then
    SLV_MODULES="${SLV_MODULES} src/ida/fcmix"
  fi

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/serial"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/fcmix_serial"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/ida/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/ida/serial"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/parallel"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/ida/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/ida/fcmix_parallel"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/ida/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/ida/parallel"
  fi

fi

# IDAS module
if test "X${IDAS_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/idas"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/idas/serial"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/idas/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/idas/serial"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/idas/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/idas/parallel"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/idas/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/idas/parallel"
  fi

fi

# KINSOL module
if test "X${KINSOL_ENABLED}" = "Xyes"; then

  SLV_MODULES="${SLV_MODULES} src/kinsol"

  if test "X${FCMIX_ENABLED}" = "Xyes" && test -d ${srcdir}/src/kinsol/fcmix ; then
    SLV_MODULES="${SLV_MODULES} src/kinsol/fcmix"
  fi

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/serial ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/serial"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_serial ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/fcmix_serial"
  fi

  if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/kinsol/serial ; then
    EXS_MODULES="${EXS_MODULES} test_examples/kinsol/serial"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/parallel"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/examples/kinsol/fcmix_parallel ; then
    EXS_MODULES="${EXS_MODULES} examples/kinsol/fcmix_parallel"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/test_examples/kinsol/parallel ; then
    EXS_MODULES="${EXS_MODULES} test_examples/kinsol/parallel"
  fi

fi

# Update the list of Makefiles to be created for all solver and example moduls
for i in ${SLV_MODULES}
do
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ${i}/Makefile"
done
for i in ${EXS_MODULES}
do
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ${i}/Makefile"
done

# Add Fortran update script to the list of additional files to be generated
if test "X${BUILD_F77_UPDATE_SCRIPT}" = "Xyes"; then
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} examples/fortran_update.sh:config/fortran_update.in"
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

if test "X${CXX_ENABLED}" = "Xyes" && test "X${CXX_OK}" = "Xyes"; then
echo "
  C++ Preprocessor:          ${CPPCXX}
  C++ Preprocessor Flags:    ${CPPFLAGS}
  C++ Compiler:              ${CXX}
  C++ Compiler Flags:        ${CXXFLAGS}"
fi

if test "X${F77_ENABLED}" = "Xyes" && test "X${F77_OK}" = "Xyes"; then
echo "
  Fortran Compiler:          ${F77}
  Fortran Compiler Flags:    ${FFLAGS}
  Fortran Linker:            ${F77_LNKR_OUT}
  Extra Fortran Libraries:   ${FLIBS}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && (test "X${MPI_C_COMP_OK}" = "Xyes" || test "X${MPI_CXX_COMP_OK}" = "Xyes" || test "X${MPI_F77_COMP_OK}" = "Xyes"); then
echo "
  MPI Root Directory:        ${MPI_ROOT_DIR}
  MPI Include Directory:     ${MPI_INC_DIR}
  MPI Library Directory:     ${MPI_LIB_DIR}
  MPI Flags:                 ${MPI_FLAGS}
  Extra MPI Libraries:       ${MPI_LIBS}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_C_COMP_OK}" = "Xyes"; then
echo "  
  Using MPI-C script?        ${USE_MPICC_SCRIPT}
  MPI-C:                     ${MPICC}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${CXX_ENABLED}" = "Xyes" && test "X${MPI_CXX_COMP_OK}" = "Xyes"; then
echo "  
  Using MPI-C++ script?      ${USE_MPICXX_SCRIPT}
  MPI-C++:                   ${MPICXX}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${F77_ENABLED}" = "Xyes" && test "X${MPI_F77_COMP_OK}" = "Xyes"; then
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
echo "  sundialsTB installed in:   ${STB_INSTDIR}"
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

if test "X${SERIAL_C_EXAMPLES}" = "Xyes"; then
  C_SERIAL_BOX="YES"
else
  C_SERIAL_BOX="NO "
fi
C_PARALLEL_BOX="N/A"
if test "X${PARALLEL_C_EXAMPLES}" = "Xyes"; then
  C_PARALLEL_BOX="YES"
elif test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_C_COMP_OK}" = "Xno"; then
  C_PARALLEL_BOX="NO "
fi

C_DEV_SERIAL_BOX="N/A"
if test "X${SERIAL_DEV_C_EXAMPLES}" = "Xyes"; then
  C_DEV_SERIAL_BOX="YES"
elif test "X${CVODES_ENABLED}" = "Xyes"; then
  C_DEV_SERIAL_BOX="NO "
fi
C_DEV_PARALLEL_BOX="N/A"
if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes"; then
  C_DEV_PARALLEL_BOX="YES"
elif test "X${CVODES_ENABLED}" = "Xyes" && test "X${MPI_C_COMP_OK}" = "Xno"; then
  C_DEV_PARALLEL_BOX="NO "
fi

CXX_SERIAL_BOX="N/A"
if test "X${SERIAL_CXX_EXAMPLES}" = "Xyes"; then
  CXX_SERIAL_BOX="YES"
elif test "X${CXX_ENABLED}" = "Xyes" && test "X${CXX_OK}" = "Xno"; then
  CXX_SERIAL_BOX="NO "
fi
CXX_PARALLEL_BOX="N/A"
if test "X${PARALLEL_CXX_EXAMPLES}" = "Xyes"; then
  CXX_PARALLEL_BOX="YES"
elif test "X${CXX_ENABLED}" = "Xyes" && test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_CXX_COMP_OK}" = "Xno"; then
  CXX_PARALLEL_BOX="NO "
fi
CXX_DEV_SERIAL_BOX="N/A"
CXX_DEV_PARALLEL_BOX="N/A"

F77_SERIAL_BOX="N/A"
if test "X${SERIAL_F77_EXAMPLES}" = "Xyes"; then
  F77_SERIAL_BOX="YES"
elif test "X${F77_ENABLED}" = "Xyes" && test "X${F77_OK}" = "Xno"; then
  F77_SERIAL_BOX="NO "
fi
F77_DEV_SERIAL_BOX="N/A"
F77_DEV_PARALLEL_BOX="N/A"
if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes"; then
  F77_PARALLEL_BOX="YES"
elif test "X${F77_ENABLED}" = "Xyes" && test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_F77_COMP_OK}" = "Xno"; then
  F77_PARALLEL_BOX="NO "
fi

if test "X${MPI_ENABLED}" = "Xyes"; then

  if test "X${RAN_ENABLE_DEV_EXAMPLES}" = "Xno"; then

    echo "            SERIAL    PARALLEL"
    echo "          ---------------------"
    echo " C       |    ${C_SERIAL_BOX}   |    ${C_PARALLEL_BOX}   |"
    echo "          ---------------------"

    if test "X${CXX_ENABLED}" = "Xyes"; then
      echo " C++     |    ${CXX_SERIAL_BOX}   |    ${CXX_PARALLEL_BOX}   |"
      echo "          ---------------------"
    fi

    if test "X${F77_ENABLED}" = "Xyes"; then
      echo " Fortran |    ${F77_SERIAL_BOX}   |    ${F77_PARALLEL_BOX}   |"
      echo "          ---------------------"
    fi

  elif test "X${RAN_ENABLE_DEV_EXAMPLES}" = "Xyes"; then

    echo "            SERIAL    PARALLEL    DEV_SERIAL    DEV_PARALLEL"
    echo "          ---------------------------------------------------"
    echo " C       |    ${C_SERIAL_BOX}   |    ${C_PARALLEL_BOX}   |      ${C_DEV_SERIAL_BOX}     |      ${C_DEV_PARALLEL_BOX}     |"
    echo "          ---------------------------------------------------"

    if test "X${CXX_ENABLED}" = "Xyes"; then
      echo " C++     |    ${CXX_SERIAL_BOX}   |    ${CXX_PARALLEL_BOX}   |      ${CXX_DEV_SERIAL_BOX}     |      ${CXX_DEV_PARALLEL_BOX}     |"
      echo "          ---------------------------------------------------"
    fi

    if test "X${F77_ENABLED}" = "Xyes"; then
      echo " Fortran |    ${F77_SERIAL_BOX}   |    ${F77_PARALLEL_BOX}   |      ${F77_DEV_SERIAL_BOX}     |      ${F77_DEV_PARALLEL_BOX}     |"
      echo "          ---------------------------------------------------"
    fi

  fi

else

  if test "X${RAN_ENABLE_DEV_EXAMPLES}" = "Xno"; then

    echo "           SERIAL"
    echo "          --------"
    echo " C       |   ${C_SERIAL_BOX}  |"
    echo "          --------"

    if test "X${CXX_ENABLED}" = "Xyes"; then
      echo " C++     |   ${CXX_SERIAL_BOX}  |"
      echo "          --------"
    fi

    if test "X${F77_ENABLED}" = "Xyes"; then
      echo " Fortran |   ${F77_SERIAL_BOX}  |"
      echo "          --------"
   fi

  elif test "X${RAN_ENABLE_DEV_EXAMPLES}" = "Xyes"; then

    echo "           SERIAL   DEV_SERIAL"
    echo "          ---------------------"
    echo " C       |   ${C_SERIAL_BOX}  |     ${C_DEV_SERIAL_BOX}    |"
    echo "          ---------------------"

    if test "X${CXX_ENABLED}" = "Xyes"; then
      echo " C++     |   ${CXX_SERIAL_BOX}  |     ${CXX_DEV_SERIAL_BOX}    |"
      echo "          ---------------------"
    fi

    if test "X${F77_ENABLED}" = "Xyes"; then
      echo " Fortran |   ${F77_SERIAL_BOX}  |     ${F77_DEV_SERIAL_BOX}    |"
      echo "          ---------------------"
    fi

  fi

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
