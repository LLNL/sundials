# aclocal.m4 generated automatically by aclocal 1.6.3 -*- Autoconf -*-

# Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
# Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

# -----------------------------------------------------------------
# $Revision: 1.18 $
# $Date: 2004-10-21 20:18:22 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see sundials/shared/LICENSE.
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

])

#------------------------------------------------------------------
# PERFORM INITIALIZATIONS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_INITIALIZE],
[

# Specify directory containing auxillary build tools and M4 files
AC_CONFIG_AUX_DIR([config])

# Make input filename DOS compatible (change config.h.in to config.hin)
AC_CONFIG_HEADERS([config.h:config.hin])

# Make user aware of copyright notice (input COPYRIGHT information)
AC_COPYRIGHT(
[
Copyright (c) 2002, The Regents of the University of California.
Produced at the Lawrence Livermore National Laboratory.
All rights reserved.
For details, see sundials/shared/LICENSE.
])

# Specify root of source tree
# Given file is guaranteed to exist in all SUNDIALS packages
AC_CONFIG_SRCDIR([shared/source/nvector.c])

# Get host information
# AC_CANONICAL_BUILD defines the following variables: build, build_cpu,
# build_vendor, and build_os
AC_CANONICAL_BUILD
# AC_CANONICAL_HOST defines the following variables: host, host_cpu,
# host_vendor, and host_os
AC_CANONICAL_HOST

# Overwrite default installation path (/usr/local)
# DEFAULT_PREFIX is defined for later use
DEFAULT_PREFIX=`pwd`
AC_PREFIX_DEFAULT(`pwd`)

# Set MAKE if necessary
# Must include @SET_MAKE@ in each Makefile.in file
# AC_SUBST is called automatically for SET_MAKE
AC_PROG_MAKE_SET

# Defines INSTALL (sets to path of "install" program)
# Also sets INSTALL_PROGRAM and INSTALL_SCRIPT
AC_PROG_INSTALL

# Set CC and F77 to default values
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

])

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

])

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
fi

])

#------------------------------------------------------------------
# CHECK C COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_CC],
[

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
elif test "X${withval}" = "Xdouble"; then
  AC_MSG_RESULT([double])
  AC_DEFINE([SUNDIALS_DOUBLE_PRECISION],[1],[])
  FLOAT_TYPE="double"
elif test "X${withval}" = "Xextended"; then
  AC_MSG_RESULT([long double])
  AC_DEFINE([SUNDIALS_EXTENDED_PRECISION],[1],[])
  FLOAT_TYPE="extended"
else
  AC_MSG_ERROR([invalid input])
fi

],
[

# Use 'double' by default
AC_MSG_RESULT([double])
AC_DEFINE([SUNDIALS_DOUBLE_PRECISION],[1],[])
FLOAT_TYPE="double"

])

AC_ARG_WITH([],[ ],[])

USER_CC="${CC}"
# unset is NOT portable so just undefine/unset CC by setting to "" (NULL)
CC=""

# Defines CC and sets GCC="yes" if CC="gcc"
# Search for C compiler given by user first
AC_PROG_CC(${USER_CC} cc gcc)

# If CC="" then abort (means did NOT find a valid C compiler)
if test "X${CC}" = "X"; then
  CC_OK="no"
  AC_MSG_ERROR([cannot find C compiler])
else

  CC_OK="yes"

  # Determine absolute pathname for specified C compiler
  CC_TEMP1="${CC}"
  CC_TEMP2=`basename "${CC}"`
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
  [AC_HELP_STRING([--with-cflags],[add extra C compiler flags])],
  [AC_MSG_RESULT([${withval}])
   USER_CFLAGS="${withval}"],
  [AC_MSG_RESULT([none])
   USER_CFLAGS=""])

  # Set CPP to command that runs C preprocessor
  # Defines CPP
  AC_PROG_CPP

  AC_MSG_CHECKING([for extra C/C++ preprocessor flags])
  AC_ARG_WITH(cppflags,
  [AC_HELP_STRING([--with-cppflags],[add extra C/C++ preprocessor flags])],
  [AC_MSG_RESULT([${withval}])
   if test "X${CPPFLAGS}" = "X"; then
     CPPFLAGS="${withval}"
   else
     CPPFLAGS="${CPPFLAGS} ${withval}"
   fi],
  [AC_MSG_RESULT([none])])

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
  [AC_HELP_STRING([--with-ldflags],[add extra linker flags])],
  [AC_MSG_RESULT([${withval}])
   if test "X${LDFLAGS}" = "X"; then
     LDFLAGS="${withval}"
   else
     LDFLAGS="${LDFLAGS} ${withval}"
   fi],
  [AC_MSG_RESULT([none])])

  AC_MSG_CHECKING([for extra libraries])
  AC_ARG_WITH(libs,
  [AC_HELP_STRING([--with-libs],[add extra libraries])],
  [AC_MSG_RESULT([${withval}])
   if test "X${LIBS}" = "X"; then
     LIBS="${withval}"
   else
     LIBS="${LIBS} ${withval}"
   fi],
  [AC_MSG_RESULT([none])])

  # Remaining tests pertain to C programming language
  AC_LANG([C])

  # Defines STDC_HEADERS if the following header files are found: stdlib.h,
  # stdarg.h, string.h, and float.h
  # We really only need stdlib.h and float.h
  AC_HEADER_STDC

  AC_CHECK_HEADERS([stdlib.h float.h math.h])

  TEMP_LIBS="${LIBS}"
  LIBS=""

  # Check if math library contains pow() and sqrt() functions (required)
  # May update LIBS (meaning add additional library, namely libm)
  AC_CHECK_LIB([m],pow,[],[AC_MSG_ERROR([cannot find pow function])])
  LIBS_TEMP="${LIBS}"
  LIBS=""
  AC_CHECK_LIB([m],sqrt,[],[AC_MSG_ERROR([cannt find sqrt function])])
  # For the sake of clarity, do NOT include the same library twice
  if test "X${LIBS}" = "X"; then
    if test "X${LIBS_TEMP}" = "X"; then
      LIBS=""
    else
      LIBS="${LIBS_TEMP}"
    fi
  else
    if test "X${LIBS_TEMP}" = "X"; then
      LIBS="${LIBS}"
    else
      if test "X${LIBS}" = "X${LIBS_TEMP}"; then
        LIBS="${LIBS}"
      else
        LIBS="${LIBS} ${LIBS_TEMP}"
      fi
    fi
  fi

  AC_MSG_CHECKING([for additional required C libraries])
  if test "X${LIBS}" = "X"; then
    if test "X${TEMP_LIBS}" = "X"; then
      LIBS=""
    else
      LIBS="${TEMP_LIBS}"
    fi
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${LIBS}])
    if test "X${TEMP_LIBS}" = "X"; then
      LIBS="${LIBS}"
    else
      LIBS="${LIBS} ${TEMP_LIBS}"
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

  # F77_FLOAT_TYPE exported via AC_SUBST (see configure.ac) - used by FCMIX Makefile's
  F77_FLOAT_TYPE="${FLOAT_TYPE}"

  # Defines EGREP and exports via AC_SUBST - used by FCMIX Makefile's
  AC_PROG_EGREP

  # Defines FGREP and exports via AC_SUBST - used by FCMIX Makefile's
  AC_PROG_FGREP

fi

])

#------------------------------------------------------------------
# DEFAULT CFLAGS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_CFLAGS],
[

# Specify "-ffloat-store" flag if using gcc on an IA-32 system (recommended)
case $host in 

  # IA-32 system running Linux
  i*-pc-linux-*)

    if test "X${GCC}" = "Xyes"; then
      if test "X${CFLAGS}" = "X"; then
        # If user wants extra precision (long double), then let program store
        # floating-point values in registers
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CFLAGS="-ffloat-store"
        fi
      else
        # If user wants extra precision (long double), then let program store
        # floating-point values in registers
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CFLAGS="${CFLAGS} -ffloat-store"
        fi
      fi
    fi

  ;;

esac

])

#------------------------------------------------------------------
# CHECK FORTRAN COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_F77],
[

AC_ARG_ENABLE([],[],[])

# Fortran compiler is only required if examples are enabled and building KINSOL
# and/or CVODE
# For now just check if examples were enabled (default)
EXAMPLES_ENABLED="yes"
AC_ARG_ENABLE(examples,
[AC_HELP_STRING([--disable-examples],[disable configuration of examples])],
[
if test "X${enableval}" = "Xno"; then
  EXAMPLES_ENABLED="no"
fi
])

AC_ARG_ENABLE([],[ ],[])

# Must have kinsol directory in order to build KINSOL
KINSOL_ENABLED="no"
if test -d ${srcdir}/kinsol ; then
  KINSOL_ENABLED="yes"
fi

# Check if user wants to skip KINSOL
AC_ARG_ENABLE(kinsol,
[AC_HELP_STRING([--disable-kinsol],[disable configuration of KINSOL])],
[
if test "X${enableval}" = "Xno"; then
  KINSOL_ENABLED="no"
fi
])

# Must have cvode directory in order to build CVODE
CVODE_ENABLED="no"
if test -d ${srcdir}/cvode ; then
  CVODE_ENABLED="yes"
fi

# Check if user wants to skip CVODE
AC_ARG_ENABLE(cvode,
[AC_HELP_STRING([--disable-cvode],[disable configuration of CVODE])],
[
if test "X${enableval}" = "Xno"; then
  CVODE_ENABLED="no"
fi
])

# Check if we need Fortran
if test "X${CVODE_ENABLED}" = "Xyes" || test "X${KINSOL_ENABLED}" = "Xyes"; then
  if test "X${EXAMPLES_ENABLED}" = "Xno"; then
    F77_ENABLED="no"
  else
    F77_ENABLED="yes"
  fi
else
  F77_ENABLED="no"
fi

AC_ARG_WITH([],[  ],[])

# Check if user wants to skip Fortran
AC_ARG_WITH([f77],
[AC_HELP_STRING([--without-f77], [disable Fortran examples])],
[
  if test "X${withval}" = "Xno"; then
    F77_ENABLED="no"
  fi
])

# Provide variable description templates for config.hin and config.h files
# Required by autoheader utility
AH_TEMPLATE([SUNDIALS_UNDERSCORE_NONE],
            [FCMIX: Do NOT append any underscores to functions names])
AH_TEMPLATE([SUNDIALS_UNDERSCORE_ONE],
            [FCMIX: Append ONE underscore to function names])
AH_TEMPLATE([SUNDIALS_UNDERSCORE_TWO],
            [FCMIX: Append TWO underscores to function names])

# Provided in case AC_F77_WRAPPERS cannot determine name mangling scheme
AC_ARG_WITH(f77underscore, 
[AC_HELP_STRING([--with-f77underscore=ARG],[specify number of underscores to append to function names (none/one/two) [AUTO]],[])],
[

UNDERSCORE_ARG_GIVEN="yes"
if test "X${withval}" = "Xnone"; then
  AC_DEFINE([SUNDIALS_UNDERSCORE_NONE],[1],[])
elif test "X${withval}" = "Xone"; then
  AC_DEFINE([SUNDIALS_UNDERSCORE_ONE],[1],[])
elif test "X${withval}" = "Xtwo"; then
  AC_DEFINE([SUNDIALS_UNDERSCORE_TWO],[1],[])
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

# Provided in case AC_F77_WRAPPERS cannot determine name mangling scheme
AC_ARG_WITH(f77case, 
[AC_HELP_STRING([--with-f77case=ARG   ],[specify case of function names (lower/upper) [AUTO]],[])],
[

CASE_ARG_GIVEN="yes"
if test "X${withval}" = "Xupper"; then
  AC_DEFINE([SUNDIALS_CASE_UPPER],[1],[])
elif test "X${withval}" = "Xlower"; then
  AC_DEFINE([SUNDIALS_CASE_LOWER],[1],[])
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
  # Only used "--with-f77case" flag, so set default number of underscores
  elif test "X${UNDERSCORE_ARG_GIVEN}" = "Xno" && test "X${CASE_ARG_GIVEN}" = "Xyes"; then
    AC_DEFINE([SUNDIALS_UNDERSCORE_ONE],[1],[])
  fi
  RUN_F77_WRAPPERS="no"
# Only call AC_F77_WRAPPERS if user did NOT use "--with-f77underscore" or "--with-f77case" flags
else
  RUN_F77_WRAPPERS="yes"
fi

# Only continue checks if Fortran support is "required"
if test "X${F77_ENABLED}" = "Xyes"; then

  F77_OK="yes"

  USER_F77="${F77}"
  # unset is NOT portable so just undefine/unset F77 by setting to "" (NULL)
  F77=""

  # Defines F77 and sets G77="yes" if F77="g77"
  # Search for Fortran compiler given by user first
  AC_PROG_F77(${USER_F77} f77 g77)

  # If F77="" then issue warning (means did NOT find a valid Fortran compiler)
  # Do NOT abort since Fortran compiler is NOT required
  if test "X${F77}" = "X"; then
    AC_MSG_WARN([cannot find Fortran compiler])
    echo ""
    echo "   Unable to find a functional Fortran compiler."
    echo ""
    echo "   Disabling serial Fortran examples..."
    echo ""
    F77_OK="no"
  else

    # Determine absolute pathname for specified Fortran compiler
    F77_TEMP1="${F77}"
    F77_TEMP2=`basename "${F77}"`
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

    # Determine how to properly mangle function names so Fortran subroutines can
    # call C functions included in SUNDIALS libraries
    # Defines C preprocessor macros F77_FUNC and F77_FUNC_
    if test "X${RUN_F77_WRAPPERS}" = "Xyes"; then
      AC_F77_WRAPPERS
    fi

    AC_MSG_CHECKING([for extra Fortran compiler flags])
    AC_ARG_WITH(fflags,
    [AC_HELP_STRING([--with-fflags],[add extra Fortran compiler flags])],
    [AC_MSG_RESULT([${withval}])
     USER_FFLAGS="${withval}"],
    [AC_MSG_RESULT([none])
     USER_FFLAGS=""])

    # Add default flag(s) to FFLAGS only if not already defined
    if test "X${FFLAGS_USER_OVERRIDE}" = "Xno"; then
      SUNDIALS_DEFAULT_FFLAGS
    fi

    if test "X${USER_FFLAGS}" = "X"; then
      FFLAGS="${FFLAGS}"
    else
      FFLAGS="${FFLAGS} ${USER_FFLAGS}"
    fi

    TEMP_FLIBS="${FLIBS}"
    FLIBS=""

    # Add any required linker flags to FLIBS
    AC_F77_LIBRARY_LDFLAGS

    AC_MSG_CHECKING([for additional required Fortran linker flags])
    if test "X${FLIBS}" = "X"; then
      if test "X${TEMP_FLIBS}" = "X"; then
        FLIBS=""
      else
        FLIBS="${TEMP_FLIBS}"
      fi
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${FLIBS}])
      if test "X${TEMP_FLIBS}" = "X"; then
        FLIBS="${FLIBS}"
      else
        if test "X${FLIBS}" = "X${TEMP_FLIBS}"; then
          FLIBS="${FLIBS}"
        else
          FLIBS="${FLIBS} ${TEMP_FLIBS}"
        fi
      fi
    fi
  fi
fi

])

#------------------------------------------------------------------
# DEFAULT FFLAGS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_DEFAULT_FFLAGS],
[

# FIXME: Should IRIX and Tru64 options overwrite FFLAGS?
case $host in

  # IA-32 system running Linux
  i*-pc-linux-*)

    if test "X${G77}" = "Xyes"; then
      if test "X${FFLAGS}" = "X"; then
        # If user wants extra precision (long double), then let program store
        # floating-point values in registers
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          FFLAGS="-ffloat-store"
        fi
      else
        # If user wants extra precision (long double), then let program store
        # floating-point values in registers
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          FFLAGS="${FFLAGS} -ffloat-store"
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

    if test "X${F77}" = "Xf77"; then
      FFLAGS="-O1"
    fi

  ;;

esac

])

#------------------------------------------------------------------
# CHECK C++ COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_CXX],
[

AC_ARG_WITH([],[   ],[])

# C++ compiler is only required if xs4c directory exists
if test -d ${srcdir}/xs4c ; then
  CXX_ENABLED="yes"
else
  CXX_ENABLED="no"
fi

# Check if user wants to skip C++
AC_ARG_WITH([cxx],
[AC_HELP_STRING([--without-cxx], [skip C++ support])],
[
if test "X${withval}" = "Xno"; then
  CXX_ENABLED="no"
fi
])

# Only continue checks if C++ support is required
if test "X${CXX_ENABLED}" = "Xyes"; then

  CXX_OK="yes"

  # CXX and CCC are common so check both
  if "X${CXX}" = "X"; then
    if "X${CCC}" = "X"; then
      USER_CXX=""
    else
      USER_CXX="${CCC}"
    fi
  else
    USER_CXX="${CXX}"
  fi
  # unset is NOT portable so just undefine/unset CXX by setting to "" (NULL)
  CXX=""

  # Defines CXX and sets GXX="yes" if CXX="g++" (GNU C++ compiler)
  AC_PROG_CXX(${USER_CXX} g++ CC gcc c++ cxx)

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
  else

    # Determine absolute pathname for specified C++ compiler
    CXX_TEMP1="${CXX}"
    CXX_TEMP2=`basename "${CXX}"`
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
    [AC_HELP_STRING([--with-cxxflags],[add extra C++ compiler flags])],
    [AC_MSG_RESULT([${withval}])
     USER_CXXFLAGS="${withval}"],
    [AC_MSG_RESULT([none])
     USER_CXXFLAGS=""])

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

    # Remaining test pertains to C++ programming language
    AC_LANG([C++])

    # Check for complex.h header file (required)
    AC_CHECK_HEADER([complex.h])

  fi
fi

])

#------------------------------------------------------------------
# DEFAULT CXXFLAGS
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_PROG_CXXFLAGS],
[

# Specify "-ffloat-store" flag if using g++ on an IA-32 system (recommended)
case $host in 

  # IA-32 system running Linux
  i*-pc-linux-*)

    if test "X${GXX}" = "Xyes"; then
      if test "X${CXXFLAGS}" = "X"; then
        # If user wants extra precision (long double), then let program store
        # floating-point values in registers
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CXXFLAGS="-ffloat-store"
        fi
      else
        # If user wants extra precision (long double), then let program store
        # floating-point values in registers
        if test "X${FLOAT_TYPE}" = "Xsingle" || test "X${FLOAT_TYPE}" = "Xdouble"; then
          CXXFLAGS="${CXXFLAGS} -ffloat-store"
        fi
      fi
    fi

  ;;

esac

])

#------------------------------------------------------------------
# CHECK MPI-C COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_MPICC],
[

AC_ARG_WITH([],[    ],[])

# MPI is only required if nvec_par directory exists (parallel NVECTOR
# implementation) and MPI support is enabled
MPI_ENABLED="no"
if test -d ${srcdir}/nvec_par ; then
  MPI_ENABLED="yes"
fi

# Check if user wants to skip MPI
AC_ARG_WITH([mpi],
[AC_HELP_STRING([--without-mpi],[skip MPI support])],
[
if test "X${withval}" = "Xno"; then
  MPI_ENABLED="no"
fi
])

# Only continue checks if MPI support is enabled
if test "X${MPI_ENABLED}" = "Xyes"; then

  # MPI root directory
  AC_ARG_WITH(mpi-root,
  [AC_HELP_STRING([--with-mpi-root=MPIROOT],[use MPI root directory])],
  [MPI_ROOT_DIR="${withval}"],
  [MPI_ROOT_DIR=""])

  # MPI include directory
  AC_ARG_WITH(mpi-incdir,
  [AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@])],
  [MPI_INC_DIR="${withval}"],
  [MPI_INC_DIR=""])

  # MPI library directory
  AC_ARG_WITH(mpi-libdir,
  [AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@])],
  [MPI_LIB_DIR="${withval}"],
  [MPI_LIB_DIR=""])

  # MPI libraries
  AC_ARG_WITH(mpi-libs,
  [AC_HELP_STRING([--with-mpi-libs=ARG],[MPI libraries @<:@-lmpi@:>@])],
  [MPI_LIBS="${withval}"],
  [MPI_LIBS=""])

  # MPI-C compiler
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
  [USE_MPICC_SCRIPT="yes"
   MPICC_COMP="mpicc"])

  # Check MPI-C compiler (either MPI compiler script or regular C compiler)
  if test "X${USE_MPICC_SCRIPT}" = "Xyes"; then
    SUNDIALS_CHECK_MPICC
  else
    MPICC_COMP="${CC}"
    MPICC="${CC}"
    SUNDIALS_CHECK_CC_WITH_MPI
  fi

fi

])

#------------------------------------------------------------------
# TEST MPI-C COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_MPICC],
[

# Test MPI-C compiler (meaning test MPICC_COMP)
# Check if MPI-C compiler can be found

# CASE 1: MPICC_COMP was found (cannot check if executable because the
# "-x" flag is NOT portable)
if test -f ${MPICC_COMP} ; then
  MPICC_COMP_EXISTS="yes"
  # Determine MPI_INC_DIR and MPI_LIB_DIR for use by Makefile
  MPI_BASE_DIR=`AS_DIRNAME(["${MPICC_COMP}"])`
  MPI_INC_DIR="${MPI_BASE_DIR}/../include"
  MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
# CASE 2: MPICC_COMP could NOT be found and MPI_ROOT_DIR was NOT specified,
# so search in PATH
else
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
      MPI_INC_DIR="${MPI_BASE_DIR}/../include"
      MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
    fi
  # CASE 3: MPICC_COMP could NOT be found, but MPI_ROOT_DIR was specified
  else
    # MPICC_COMP should really only contain an executable name
    # Found location of MPICC_COMP
    if test -f ${MPI_ROOT_DIR}/bin/${MPICC_COMP} ; then
      MPICC_COMP_EXISTS="yes"
      MPICC_COMP="${MPI_ROOT_DIR}/bin/${MPICC_COMP}"
      MPI_INC_DIR="${MPI_ROOT_DIR}/include"
      MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    # Could NOT find MPICC_COMP anywhere
    else
      MPICC_COMP_EXISTS="no"
      MPICC_COMP=""
    fi
  fi
fi

# Issue warning message if MPICC_COMP does NOT exist, else set MPICC
if test "X${MPICC_COMP_EXISTS}" = "Xyes"; then
  MPICC="${MPICC_COMP}"
  MPI_C_COMP_OK="yes"
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
fi

])

#------------------------------------------------------------------
# TEST C COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_CC_WITH_MPI],
[

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
    echo "   Disabling the parallel NVECTOR module and all parallel examples..."
    echo ""
  # MPI root directory was given so set MPI_INC_DIR accordingly
  # Update CPPFLAGS
  else
    MPI_INC_DIR="${MPI_ROOT_DIR}/include"
    if test "X${CPPFLAGS}" = "X"; then
      CPPFLAGS="-I${MPI_INC_DIR}"
    else
      CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
    fi
  fi
# MPI include directory was specified so update CPPFLAGS
else
  if test "X${CPPFLAGS}" = "X"; then
    CPPFLAGS="-I${MPI_INC_DIR}"
  else
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
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

  # Check if user specified which MPI libraries must be included
  # If no libraries are given, then assume libmpi.[a/so] exists
  if test "X${MPI_LIBS}" = "X"; then
    MPI_LIBS="-lmpi"
    if test "X${LIBS}" = "X"; then
      LIBS="${MPI_LIBS}"
    else
      LIBS="${LIBS} ${MPI_LIBS}"
    fi
  # MPI libraries were specified so update LIBS
  else
    if test "X${LIBS}" = "X"; then
      LIBS="${MPI_LIBS}"
    else
      LIBS="${LIBS} ${MPI_LIBS}"
    fi
  fi

  # Testing C language program with MPI
  AC_LANG([C])

  AC_MSG_CHECKING([if C compiler can compile MPI programs])
  AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[#include "mpi.h"]],[[int c; char **v; MPI_Init(&c,&v);]])],
  [AC_MSG_RESULT(yes)
   MPI_C_COMP_OK="yes"],
  [AC_MSG_RESULT(no)
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
   MPI_C_COMP_OK="no"])
else
  MPI_C_COMP_OK="no"
fi

# Restore CPPFLAGS, LDFLAGS and LIBS
CPPFLAGS="${SAVED_CPPFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

])

#------------------------------------------------------------------
# CHECK MPI-F77 COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_MPIF77],
[

# Check MPI-F77 compiler only if MPI support is enabled and F77_ENABLED="yes"
# Do NOT base check on F77_OK flag because even if a valid serial Fortran
# compiler could NOT be found the MPIF77 compiler script may still work
# F77_ENABLED flag also indicates if "--without-f77" flag was given by user
if test "X${MPI_ENABLED}" = "Xyes" && test "X${F77_ENABLED}" = "Xyes" ; then

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
  [USE_MPIF77_SCRIPT="yes"
   MPIF77_COMP="mpif77"])

  # Check MPI-Fortran compiler (either MPI compiler script or regular Fortran compiler)
  if test "X${USE_MPIF77_SCRIPT}" = "Xyes"; then
    SUNDIALS_CHECK_MPIF77
  else
    MPIF77_COMP="${F77}"
    MPIF77="${F77}"
    SUNDIALS_CHECK_F77_WITH_MPI
  fi

fi

])

#------------------------------------------------------------------
# TEST MPI-FORTRAN COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_MPIF77],
[

# Test the MPI-Fortran compiler (meaning test MPIF77_COMP)
# Check if MPI-Fortran compiler can be found

# CASE 1: MPIF77_COMP was found (cannot check if executable because the
# "-x" flag is NOT portable)
if test -f ${MPIF77_COMP} ; then
  MPIF77_COMP_EXISTS="yes"
  # Determine MPI_INC_DIR and MPI_LIB_DIR for use by Makefile
  MPI_BASE_DIR=`AS_DIRNAME(["${MPIF77_COMP}"])`
  MPI_INC_DIR="${MPI_BASE_DIR}/../include"
  MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
# CASE 2: MPIF77_COMP could NOT be found and MPI_ROOT_DIR was NOT specified,
# so search in PATH
else
  if test "X${MPI_ROOT_DIR}" = "X"; then
    # Try to find location of executable (perhaps directory was entered
    # incorrectly)
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
      MPI_INC_DIR="${MPI_BASE_DIR}/../include"
      MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
    fi
  # CASE 3: MPIF77_COMP could NOT be found, but MPI_ROOT_DIR was specified
  else
    # MPIF77_COMP should really only contain an executable name
    # Found location of MPIF77_COMP
    if test -f ${MPI_ROOT_DIR}/bin/${MPIF77_COMP} ; then
      MPIF77_COMP_EXISTS="yes"
      MPIF77_COMP="${MPI_ROOT_DIR}/bin/${MPIF77_COMP}"
      MPI_INC_DIR="${MPI_ROOT_DIR}/include"
      MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    # Could NOT find MPIF77_COMP anywhere
    else
      MPIF77_COMP_EXISTS="no"
      MPIF77_COMP=""
    fi
  fi
fi

# Issue warning message if MPIF77_COMP does NOT exist, else set MPIF77
if test "X${MPIF77_COMP_EXISTS}" = "Xyes"; then
  MPIF77="${MPIF77_COMP}"
  MPI_F77_COMP_OK="yes"
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
fi

])

#------------------------------------------------------------------
# TEST FORTRAN COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_F77_WITH_MPI],
[

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
      echo "   Try using --with-mpif77 to specify a MPI-Fortran compiler script,"
      echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
      echo "   to specify the locations of all relevant MPI files, or"
      echo "   --with-mpi-root to specify the base installation directory"
      echo "   of the MPI implementation to be used."
      echo ""
      echo "   Disabling all parallel Fortran examples..."
      echo ""
    # MPI root directory was given so set MPI_INC_DIR accordingly
    # Update FFLAGS
    else
      MPI_INC_DIR="${MPI_ROOT_DIR}/include"
      if test "X${FFLAGS}" = "X"; then
        FFLAGS="-I${MPI_INC_DIR}"
      else
        FFLAGS="${FFLAGS} -I${MPI_INC_DIR}"
      fi
    fi
  # MPI include directory was specified so update FFLAGS
  else
    if test "X${FFLAGS}" = "X"; then
      FFLAGS="-I${MPI_INC_DIR}"
    else
      FFLAGS="${FFLAGS} -I${MPI_INC_DIR}"
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

    # Check if user specified which MPI libraries must be included
    # If no libraries are given, then assume libmpi.[a/so] exists
    if test "X${MPI_LIBS}" = "X"; then
      MPI_LIBS="-lmpi"
      if test "X${LIBS}" = "X"; then
        LIBS="${MPI_LIBS}"
      else
        LIBS="${LIBS} ${MPI_LIBS}"
      fi
    # MPI libraries were specified so update LIBS
    else
      if test "X${LIBS}" = "X"; then
        LIBS="${MPI_LIBS}"
      else
        LIBS="${LIBS} ${MPI_LIBS}"
      fi
    fi

    # Testing Fortran language program with MPI
    AC_LANG([Fortran 77])

    AC_MSG_CHECKING([if Fortran compiler can compile MPI programs])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([],
    [      
          INCLUDE "mpif.h"
          CALL MPI_INIT(IER)
    ])],
    [AC_MSG_RESULT(yes)
     MPI_F77_COMP_OK="yes"],
    [AC_MSG_RESULT(no)
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
     MPI_F77_COMP_OK="no"])
  else
    MPI_F77_COMP_OK="no"
  fi
else
  AC_MSG_WARN([cannot find Fortran compiler])
  MPI_F77_COMP_OK="no"
fi

# Restore FFLAGS, LDFLAGS and LIBS
FFLAGS="${SAVED_FFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

])

#------------------------------------------------------------------
# CHECK MPI-C++ COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_SET_MPICXX],
[

# Check MPI-C++ compiler only if MPI support is enabled and CXX_ENABLED="yes"
# Do NOT base check on CXX_OK flag because even if a valid serial C++ compiler
# could NOT be found the MPIC++ compiler script may still work
# CXX_ENABLED flag also indicates if "--without-cxx" flag was given by user
if test "X${MPI_ENABLED}" = "Xyes" && test "X${CXX_ENABLED}" = "Xyes"; then

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

  # Check MPI-C++ compiler (either MPI compiler script or regular C++ compiler)
  if test "X${USE_MPICXX_SCRIPT}" = "Xyes"; then
    SUNDIALS_CHECK_MPICXX
  else
    MPICXX_COMP="${CXX}"
    MPICXX="${CXX}"
    SUNDIALS_CHECK_CXX_WITH_MPI
  fi

fi

AC_ARG_WITH([],[     ],[])

])

#------------------------------------------------------------------
# TEST MPI-C++ COMPILER
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_MPICXX],
[

# Test the MPI-C++ compiler (meaning test MPICXX_COMP)
# Check if MPI-C++ compiler can be found

# CASE 1: MPICXX_COMP was found (cannot check if executable because the
# "-x" flag is NOT portable)
if test -f ${MPICXX_COMP} ; then
  MPICXX_COMP_EXISTS="yes"
  # Determine MPI_INC_DIR and MPI_LIB_DIR for use by Makefile
  MPI_BASE_DIR=`AS_DIRNAME(["${MPICXX_COMP}"])`
  MPI_INC_DIR="${MPI_BASE_DIR}/../include"
  MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
# CASE 2: MPICXX_COMP could NOT be found and MPI_ROOT_DIR was NOT specified,
# so search in PATH
else
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
      MPI_INC_DIR="${MPI_BASE_DIR}/../include"
      MPI_LIB_DIR="${MPI_BASE_DIR}/../lib"
    fi
  # CASE 3: MPICXX_COMP could NOT be found, but MPI_ROOT_DIR was specified
  else
    # MPICXX_COMP should really only contain an executable name
    # Found location of MPICXX_COMP
    if test -f ${MPI_ROOT_DIR}/bin/${MPICXX_COMP} ; then
      MPICXX_COMP_EXISTS="yes"
      MPICXX_COMP="${MPI_ROOT_DIR}/bin/${MPICXX_COMP}"
      MPI_INC_DIR="${MPI_ROOT_DIR}/include"
      MPI_LIB_DIR="${MPI_ROOT_DIR}/lib"
    # Could NOT find MPICXX_COMP anywhere
    else
      MPICXX_COMP_EXISTS="no"
      MPICXX_COMP=""
    fi
  fi
fi

# Issue warning message if MPICXX_COMP does NOT exist, else set MPICXX
if test "X${MPICXX_COMP_EXISTS}" = "Xyes"; then
  MPICXX="${MPICXX_COMP}"
  MPI_CXX_COMP_OK="yes"
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
fi

])

#------------------------------------------------------------------
# TEST C++ COMPILER WITH MPI
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_CHECK_CXX_WITH_MPI],
[

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
      echo "   Try using --with-mpicxx to specify a MPI-C++ compiler script,"
      echo "   --with-mpi-incdir, --with-mpi-libdir and --with-mpi-libs"
      echo "   to specify the locations of all relevant MPI files, or"
      echo "   --with-mpi-root to specify the base installation directory"
      echo "   of the MPI implementation to be used."
      echo ""
      echo "   Disabling dependent modules..."
      echo ""
    # MPI root directory was given so set MPI_INC_DIR accordingly
    # Update CPPFLAGS
    else
      MPI_INC_DIR="${MPI_ROOT_DIR}/include"
      if test "X${CPPFLAGS}" = "X"; then
        CPPFLAGS="-I${MPI_INC_DIR}"
      else
        CPPLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
      fi
    fi
  # MPI include directory was specified so update CPPFLAGS
  else
    if test "X${CPPFLAGS}" = "X"; then
      CPPFLAGS="-I${MPI_INC_DIR}"
    else
      CPPFLAGS="${CPPFLAGS} -I${MPI_INC_DIR}"
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

    # Check if user specified which MPI libraries must be included
    # If no libraries are given, then assume libmpi.[a/so] exists
    if test "X${MPI_LIBS}" = "X"; then
      MPI_LIBS="-lmpi"
      if test "X${LIBS}" = "X"; then
        LIBS="${MPI_LIBS}"
      else
        LIBS="${LIBS} ${MPI_LIBS}"
      fi
    # MPI libraries were specified so update LIBS
    else
      if test "X${LIBS}" = "X"; then
        LIBS="${MPI_LIBS}"
      else
        LIBS="${LIBS} ${MPI_LIBS}"
      fi
    fi

    # Testing C++ language program with MPI
    AC_LANG([C++])

    AC_MSG_CHECKING([if C++ compiler can compile MPI programs])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include "mpi.h"]],[[int c; car **v; MPI_Init(&c,&v);]])],
    [AC_MSG_RESULT(yes)
     MPI_CXX_COMP_OK="yes"],
    [AC_MSG_RESULT(no)
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
     MPI_CXX_COMP_OK="no"])
  else
    MPI_CXX_COMP_OK="no"
  fi
else
  AC_MSG_WARN([cannot find C++ compiler])
  MPI_CXX_COMP_OK="no"
fi

# Restore CPPFLAGS, LDFLAGS and LIBS
CPPFLAGS="${SAVED_CPPFLAGS}"
LDFLAGS="${SAVED_LDFLAGS}"
LIBS="${SAVED_LIBS}"

])

#------------------------------------------------------------------
# ENABLE MODULES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_ENABLE_MODULES],
[

# Must have cvode directory in order to build CVODE
CVODE_ENABLED="no"
if test -d ${srcdir}/cvode ; then
  CVODE_ENABLED="yes"
fi

# Check if user wants to skip CVODE
AC_ARG_ENABLE(cvode,
[AC_HELP_STRING([--disable-cvode],[disable configuration of CVODE])],
[
if test "X${enableval}" = "Xno"; then
  CVODE_ENABLED="no"
fi
])

# Must have cvodes directory in order to build CVODES
CVODES_ENABLED="no"
if test -d ${srcdir}/cvodes ; then
  CVODES_ENABLED="yes"
fi

# Check if user wants to skip CVODES
AC_ARG_ENABLE(cvodes,
[AC_HELP_STRING([--disable-cvodes],[disable configuration of CVODES])],
[
if test "X${enableval}" = "Xno"; then
  CVODES_ENABLED="no"
fi
])

# Must have ida directory in order to build IDA
IDA_ENABLED="no"
if test -d ${srcdir}/ida  ; then
  IDA_ENABLED="yes"
fi

# Check if user wants to skip IDA
AC_ARG_ENABLE(ida,
[AC_HELP_STRING([--disable-ida],[disable configuration of IDA])],
[
if test "X${enableval}" = "Xno"; then
  IDA_ENABLED="no"
fi
])

# Must have kinsol directory in order to build KINSOL
KINSOL_ENABLED="no"
if test -d ${srcdir}/kinsol ; then
  KINSOL_ENABLED="yes"
fi

# Check if user wants to skip KINSOL
AC_ARG_ENABLE(kinsol,
[AC_HELP_STRING([--disable-kinsol],[disable configuration of KINSOL])],
[
if test "X${enableval}" = "Xno"; then
  KINSOL_ENABLED="no"
fi
])

IDAS_ENABLED="no"

AC_ARG_ENABLE([],[  ],[])

])

#------------------------------------------------------------------
# ENABLE DEVELOPER MODULES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_ENABLE_DEV_MODULES],
[

# Must have idas directory in order to build IDAS
IDAS_ENABLED="no"
if test -d ${srcdir}/idas ; then
  IDAS_ENABLED="yes"
fi

# Check if user wants to skip IDAS
AC_ARG_ENABLE(idas,
[AC_HELP_STRING([--disable-idas],[disable configuration of IDAS])],
[
if test "X${enableval}" = "Xno"; then
  IDAS_ENABLED="no"
fi
])

AC_ARG_ENABLE([],[   ],[])

])

#------------------------------------------------------------------
# ENABLE EXAMPLES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_ENABLE_EXAMPLES],
[

# Examples are built by default, but user can override
EXAMPLES_ENABLED="yes"
AC_ARG_ENABLE(examples,
[AC_HELP_STRING([--disable-examples],[disable configuration of examples])],
[
if test "X${enableval}" = "Xno"; then
  EXAMPLES_ENABLED="no"
fi
])

# Check if we need to build the Fortran examples
# Recall F77_ENABLED is set in SUNDIALS_SET_F77 based upon CVODE_ENABLED,
# KINSOL_ENABLED and EXAMPLES_ENABLED, so the check here is rather simple
F77_EXAMPLES="no";
if test "X${F77_ENABLED}" = "Xyes"; then
  F77_EXAMPLES="yes"
fi

# Check if we need to build the C++ examples
CXX_EXAMPLES="no"
if test "X${EXAMPLES_ENABLED}" = "Xyes" && test "X${CXX_ENABLED}" = "Xyes"; then
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
if test "X${EXAMPLES_ENABLED}" = "Xyes"; then

  AC_MSG_CHECKING([if we can build serial C examples])
  if test "X${CC_OK}" = "Xyes"; then
    SERIAL_C_EXAMPLES="yes"
  fi
  AC_MSG_RESULT([${SERIAL_C_EXAMPLES}])

  if test "X${MPI_ENABLED}" = "Xyes"; then

    AC_MSG_CHECKING([if we can build parallel C examples])
    if test "X${MPI_C_COMP_OK}" = "Xyes"; then
      PARALLEL_C_EXAMPLES="yes"
    fi
    AC_MSG_RESULT([${PARALLEL_C_EXAMPLES}])

  fi

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

# Disable developer examples be default so output logic still works
PARALLEL_DEV_C_EXAMPLES="no"

RAN_ENABLE_DEV_EXAMPLES="no"

])

#------------------------------------------------------------------
# ENABLE DEVELOPER EXAMPLES
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_ENABLE_DEV_EXAMPLES],
[

RAN_ENABLE_DEV_EXAMPLES="yes"

# Check if developer examples can actually be built
PARALLEL_DEV_C_EXAMPLES="no"

# Check C developer examples
if test "X${EXAMPLES_ENABLED}" = "Xyes" && test "X${CVODES_ENABLED}" = "Xyes"; then

  if test "X${MPI_ENABLED}" = "Xyes"; then
    AC_MSG_CHECKING([if we can build parallel C developer examples])
    if test "X${MPI_C_COMP_OK}" = "Xyes"; then
      PARALLEL_DEV_C_EXAMPLES="yes"
    fi
    AC_MSG_RESULT([${PARALLEL_DEV_C_EXAMPLES}])
  fi

fi

])

#------------------------------------------------------------------
# BUILD MODULE LIST
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_BUILD_MODULES_LIST],
[

EXAMPLE_MODULES=""
SUNDIALS_MAKEFILES="Makefile"

# SHARED module (always required)
MODULES="shared/source"
SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} shared/source/Makefile"

# NVECTOR modules (serial and parallel)
NVEC_MODULES=""
if test -d ${srcdir}/nvec_ser ; then
  NVEC_MODULES="nvec_ser"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} nvec_ser/Makefile"
fi

if test -d ${srcdir}/nvec_par && test "X${MPI_C_COMP_OK}" = "Xyes"; then
  NVEC_MODULES="${NVEC_MODULES} nvec_par"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} nvec_par/Makefile"
fi

# CVODE module
if test "X${CVODE_ENABLED}" = "Xyes"; then

  MODULES="${MODULES} cvode/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/source/Makefile"

  MODULES="${MODULES} cvode/fcmix"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvode/examples_ser ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} cvode/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvode/examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} cvode/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/examples_par/Makefile"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvode/fcmix/examples_ser ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} cvode/fcmix/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvode/fcmix/examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} cvode/fcmix/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvode/fcmix/examples_par/Makefile"
  fi

fi

# CVODES module
if test "X${CVODES_ENABLED}" = "Xyes"; then

  MODULES="${MODULES} cvodes/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvodes/examples_ser ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} cvodes/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvodes/examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} cvodes/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/examples_par/Makefile"
  fi

  if test "X${PARALLEL_DEV_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/cvodes/test_examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} cvodes/test_examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} cvodes/test_examples_par/Makefile"
  fi

fi

# IDA module
if test "X${IDA_ENABLED}" = "Xyes"; then

  MODULES="${MODULES} ida/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/ida/examples_ser ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} ida/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/ida/examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} ida/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} ida/examples_par/Makefile"
  fi

fi

# IDAS module
if test "X${IDAS_ENABLED}" = "Xyes"; then

  MODULES="${MODULES} idas/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/source/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/idas/examples_ser ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} idas/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/idas/examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} idas/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} idas/examples_par/Makefile"
  fi

fi

# KINSOL module
if test "X${KINSOL_ENABLED}" = "Xyes"; then

  MODULES="${MODULES} kinsol/source"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/source/Makefile"

  MODULES="${MODULES} kinsol/fcmix"
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/Makefile"

  if test "X${SERIAL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/kinsol/examples_ser ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} kinsol/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_C_EXAMPLES}" = "Xyes" && test -d ${srcdir}/kinsol/examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} kinsol/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/examples_par/Makefile"
  fi

  if test "X${SERIAL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/kinsol/fcmix/examples_ser ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} kinsol/fcmix/examples_ser"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/examples_ser/Makefile"
  fi

  if test "X${PARALLEL_F77_EXAMPLES}" = "Xyes" && test -d ${srcdir}/kinsol/fcmix/examples_par ; then
    EXAMPLE_MODULES="${EXAMPLE_MODULES} kinsol/fcmix/examples_par"
    SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} kinsol/fcmix/examples_par/Makefile"
  fi

fi

# Fortran update script
if test "X${BUILD_F77_UPDATE_SCRIPT}" = "Xyes"; then
  SUNDIALS_MAKEFILES="${SUNDIALS_MAKEFILES} config/fortran_update"
fi

])

#------------------------------------------------------------------
# DETERMINE INSTALLATION PATH
#------------------------------------------------------------------

AC_DEFUN([SUNDIALS_INSTALL_PATH],
[

AC_MSG_CHECKING([for 'include' directory])
if test "X${prefix}" = "XNONE"; then
  SUNDIALS_INC_DIR="${DEFAULT_PREFIX}/include"
else
  SUNDIALS_INC_DIR="${prefix}/include"
fi
AC_MSG_RESULT([${SUNDIALS_INC_DIR}])

AC_MSG_CHECKING([for 'lib' directory])
if test "X${exec_prefix}" = "XNONE"; then
  if test "X${prefix}" = "XNONE"; then
    SUNDIALS_LIB_DIR="${DEFAULT_PREFIX}/lib"
  else
    SUNDIALS_LIB_DIR="${prefix}/lib"
  fi
else
  SUNDIALS_LIB_DIR="${exec_prefix}/lib"
fi
AC_MSG_RESULT([${SUNDIALS_LIB_DIR}])

if test -d ${SUNDIALS_INC_DIR} ; then
  :
else
  AS_MKDIR_P(${SUNDIALS_INC_DIR})
fi

if test -d ${SUNDIALS_LIB_DIR} ; then
  :
else
  AS_MKDIR_P(${SUNDIALS_LIB_DIR})
fi

])

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
Configuration:
--------------

  Host System:             ${host}
  Build System:            ${build}
  Source Code Location:    ${srcdir}
  Install Path (include):  ${SUNDIALS_INC_DIR}
  Install Path (lib):      ${SUNDIALS_LIB_DIR}

  C Preprocessor:          ${CPP} ${CPPFLAGS}
  C Compiler:	           ${CC} ${CFLAGS}
  C Linker:                ${CC} ${LDFLAGS} ${LIBS}"

if test "X${CXX_ENABLED}" = "Xyes" && test "X${CXX_OK}" = "Xyes"; then
echo "
  C++ Preprocessor:        ${CPPCXX} ${CPPFLAGS}
  C++ Compiler:            ${CXX} ${CXXFLAGS}"
fi

if test "X${F77_ENABLED}" = "Xyes" && test "X${F77_OK}" = "Xyes"; then
echo "
  F77 Compiler:            ${F77} ${FFLAGS}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_C_COMP_OK}" = "Xyes"; then
echo "  
  MPI-C:                   ${MPICC}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${CXX_ENABLED}" = "Xyes" && test "X${MPI_CXX_COMP_OK}" = "Xyes"; then
echo "  
  MPI-C++:                 ${MPICXX}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${F77_ENABLED}" = "Xyes" && test "X${MPI_F77_COMP_OK}" = "Xyes"; then
echo "  
  MPI-F77:                 ${MPIF77}"
fi

if test "X${MPI_ENABLED}" = "Xyes" && test "X${MPI_C_COMP_OK}" = "Xyes" && test "X${USE_MPICC_SCRIPT}" = "Xno"; then
echo "
  MPI Include Directory:   ${MPI_INC_DIR}
  MPI Library Directory:   ${MPI_LIB_DIR}
  MPI Libraries:           ${MPI_LIBS}"
fi

echo "  
  Type 'make' and then 'make install' to build and install ${PACKAGE_STRING}"

echo "
Modules:
--------
"

if test "X${KINSOL_ENABLED}" = "Xyes"; then
echo "  KINSOL"
fi
if test "X${CVODE_ENABLED}" = "Xyes"; then
echo "  CVODE"
fi
if test "X${CVODES_ENABLED}" = "Xyes"; then
echo "  CVODES"
fi
if test "X${IDA_ENABLED}" = "Xyes"; then
echo "  IDA"
fi
if test "X${IDAS_ENABLED}" = "Xyes"; then
echo "  IDAS"
fi

if test "X${EXAMPLES_ENABLED}" = "Xyes"; then

echo "
Examples:
---------

  Type 'make examples' to build the following examples:
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
----------------------------------
Finished SUNDIALS Configure Script
----------------------------------
"

])

# libtool.m4 - Configure libtool for the host system. -*-Shell-script-*-

# serial 46 AC_PROG_LIBTOOL

AC_DEFUN([AC_PROG_LIBTOOL],
[AC_REQUIRE([AC_LIBTOOL_SETUP])dnl

# This can be used to rebuild libtool when needed
LIBTOOL_DEPS="$ac_aux_dir/ltmain.sh"

# Always use our own libtool.
LIBTOOL='$(SHELL) $(top_builddir)/libtool'
AC_SUBST(LIBTOOL)dnl

# Prevent multiple expansion
define([AC_PROG_LIBTOOL], [])
])

AC_DEFUN([AC_LIBTOOL_SETUP],
[AC_PREREQ(2.13)dnl
AC_REQUIRE([AC_ENABLE_SHARED])dnl
AC_REQUIRE([AC_ENABLE_STATIC])dnl
AC_REQUIRE([AC_ENABLE_FAST_INSTALL])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_LD])dnl
AC_REQUIRE([AC_PROG_LD_RELOAD_FLAG])dnl
AC_REQUIRE([AC_PROG_NM])dnl
AC_REQUIRE([LT_AC_PROG_SED])dnl

AC_REQUIRE([AC_PROG_LN_S])dnl
AC_REQUIRE([AC_DEPLIBS_CHECK_METHOD])dnl
AC_REQUIRE([AC_OBJEXT])dnl
AC_REQUIRE([AC_EXEEXT])dnl
dnl

_LT_AC_PROG_ECHO_BACKSLASH
# Only perform the check for file, if the check method requires it
case $deplibs_check_method in
file_magic*)
  if test "$file_magic_cmd" = '$MAGIC_CMD'; then
    AC_PATH_MAGIC
  fi
  ;;
esac

AC_CHECK_TOOL(RANLIB, ranlib, :)
AC_CHECK_TOOL(STRIP, strip, :)

ifdef([AC_PROVIDE_AC_LIBTOOL_DLOPEN], enable_dlopen=yes, enable_dlopen=no)
ifdef([AC_PROVIDE_AC_LIBTOOL_WIN32_DLL],
enable_win32_dll=yes, enable_win32_dll=no)

AC_ARG_ENABLE(libtool-lock,
  [  --disable-libtool-lock  avoid locking (might break parallel builds)])
test "x$enable_libtool_lock" != xno && enable_libtool_lock=yes

# Some flags need to be propagated to the compiler or linker for good
# libtool support.
case $host in
*-*-irix6*)
  # Find out which ABI we are using.
  echo '[#]line __oline__ "configure"' > conftest.$ac_ext
  if AC_TRY_EVAL(ac_compile); then
    case `/usr/bin/file conftest.$ac_objext` in
    *32-bit*)
      LD="${LD-ld} -32"
      ;;
    *N32*)
      LD="${LD-ld} -n32"
      ;;
    *64-bit*)
      LD="${LD-ld} -64"
      ;;
    esac
  fi
  rm -rf conftest*
  ;;

*-*-sco3.2v5*)
  # On SCO OpenServer 5, we need -belf to get full-featured binaries.
  SAVE_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -belf"
  AC_CACHE_CHECK([whether the C compiler needs -belf], lt_cv_cc_needs_belf,
    [AC_LANG_SAVE
     AC_LANG_C
     AC_TRY_LINK([],[],[lt_cv_cc_needs_belf=yes],[lt_cv_cc_needs_belf=no])
     AC_LANG_RESTORE])
  if test x"$lt_cv_cc_needs_belf" != x"yes"; then
    # this is probably gcc 2.8.0, egcs 1.0 or newer; no need for -belf
    CFLAGS="$SAVE_CFLAGS"
  fi
  ;;

ifdef([AC_PROVIDE_AC_LIBTOOL_WIN32_DLL],
[*-*-cygwin* | *-*-mingw* | *-*-pw32*)
  AC_CHECK_TOOL(DLLTOOL, dlltool, false)
  AC_CHECK_TOOL(AS, as, false)
  AC_CHECK_TOOL(OBJDUMP, objdump, false)

  # recent cygwin and mingw systems supply a stub DllMain which the user
  # can override, but on older systems we have to supply one
  AC_CACHE_CHECK([if libtool should supply DllMain function], lt_cv_need_dllmain,
    [AC_TRY_LINK([],
      [extern int __attribute__((__stdcall__)) DllMain(void*, int, void*);
      DllMain (0, 0, 0);],
      [lt_cv_need_dllmain=no],[lt_cv_need_dllmain=yes])])

  case $host/$CC in
  *-*-cygwin*/gcc*-mno-cygwin*|*-*-mingw*)
    # old mingw systems require "-dll" to link a DLL, while more recent ones
    # require "-mdll"
    SAVE_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS -mdll"
    AC_CACHE_CHECK([how to link DLLs], lt_cv_cc_dll_switch,
      [AC_TRY_LINK([], [], [lt_cv_cc_dll_switch=-mdll],[lt_cv_cc_dll_switch=-dll])])
    CFLAGS="$SAVE_CFLAGS" ;;
  *-*-cygwin* | *-*-pw32*)
    # cygwin systems need to pass --dll to the linker, and not link
    # crt.o which will require a WinMain@16 definition.
    lt_cv_cc_dll_switch="-Wl,--dll -nostartfiles" ;;
  esac
  ;;
  ])
esac

_LT_AC_LTCONFIG_HACK

])

# AC_LIBTOOL_HEADER_ASSERT
# ------------------------
AC_DEFUN([AC_LIBTOOL_HEADER_ASSERT],
[AC_CACHE_CHECK([whether $CC supports assert without backlinking],
    [lt_cv_func_assert_works],
    [case $host in
    *-*-solaris*)
      if test "$GCC" = yes && test "$with_gnu_ld" != yes; then
        case `$CC --version 2>/dev/null` in
        [[12]].*) lt_cv_func_assert_works=no ;;
        *)        lt_cv_func_assert_works=yes ;;
        esac
      fi
      ;;
    esac])

if test "x$lt_cv_func_assert_works" = xyes; then
  AC_CHECK_HEADERS(assert.h)
fi
])# AC_LIBTOOL_HEADER_ASSERT

# _LT_AC_CHECK_DLFCN
# --------------------
AC_DEFUN([_LT_AC_CHECK_DLFCN],
[AC_CHECK_HEADERS(dlfcn.h)
])# _LT_AC_CHECK_DLFCN

# AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE
# ---------------------------------
AC_DEFUN([AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE],
[AC_REQUIRE([AC_CANONICAL_HOST])
AC_REQUIRE([AC_PROG_NM])
AC_REQUIRE([AC_OBJEXT])
# Check for command to grab the raw symbol name followed by C symbol from nm.
AC_MSG_CHECKING([command to parse $NM output])
AC_CACHE_VAL([lt_cv_sys_global_symbol_pipe], [dnl

# These are sane defaults that work on at least a few old systems.
# [They come from Ultrix.  What could be older than Ultrix?!! ;)]

# Character class describing NM global symbol codes.
symcode='[[BCDEGRST]]'

# Regexp to match symbols that can be accessed directly from C.
sympat='\([[_A-Za-z]][[_A-Za-z0-9]]*\)'

# Transform the above into a raw symbol and a C symbol.
symxfrm='\1 \2\3 \3'

# Transform an extracted symbol line into a proper C declaration
lt_cv_global_symbol_to_cdecl="sed -n -e 's/^. .* \(.*\)$/extern char \1;/p'"

# Transform an extracted symbol line into symbol name and symbol address
lt_cv_global_symbol_to_c_name_address="sed -n -e 's/^: \([[^ ]]*\) $/  {\\\"\1\\\", (lt_ptr) 0},/p' -e 's/^$symcode \([[^ ]]*\) \([[^ ]]*\)$/  {\"\2\", (lt_ptr) \&\2},/p'"

# Define system-specific variables.
case $host_os in
aix*)
  symcode='[[BCDT]]'
  ;;
cygwin* | mingw* | pw32*)
  symcode='[[ABCDGISTW]]'
  ;;
hpux*) # Its linker distinguishes data from code symbols
  lt_cv_global_symbol_to_cdecl="sed -n -e 's/^T .* \(.*\)$/extern char \1();/p' -e 's/^$symcode* .* \(.*\)$/extern char \1;/p'"
  lt_cv_global_symbol_to_c_name_address="sed -n -e 's/^: \([[^ ]]*\) $/  {\\\"\1\\\", (lt_ptr) 0},/p' -e 's/^$symcode* \([[^ ]]*\) \([[^ ]]*\)$/  {\"\2\", (lt_ptr) \&\2},/p'"
  ;;
irix* | nonstopux*)
  symcode='[[BCDEGRST]]'
  ;;
osf*)
  symcode='[[BCDEGQRST]]'
  ;;
solaris* | sysv5*)
  symcode='[[BDT]]'
  ;;
sysv4)
  symcode='[[DFNSTU]]'
  ;;
esac

# Handle CRLF in mingw tool chain
opt_cr=
case $host_os in
mingw*)
  opt_cr=`echo 'x\{0,1\}' | tr x '\015'` # option cr in regexp
  ;;
esac

# If we're using GNU nm, then use its standard symbol codes.
if $NM -V 2>&1 | egrep '(GNU|with BFD)' > /dev/null; then
  symcode='[[ABCDGISTW]]'
fi

# Try without a prefix undercore, then with it.
for ac_symprfx in "" "_"; do

  # Write the raw and C identifiers.
lt_cv_sys_global_symbol_pipe="sed -n -e 's/^.*[[ 	]]\($symcode$symcode*\)[[ 	]][[ 	]]*\($ac_symprfx\)$sympat$opt_cr$/$symxfrm/p'"

  # Check to see that the pipe works correctly.
  pipe_works=no
  rm -f conftest*
  cat > conftest.$ac_ext <<EOF
#ifdef __cplusplus
extern "C" {
#endif
char nm_test_var;
void nm_test_func(){}
#ifdef __cplusplus
}
#endif
int main(){nm_test_var='a';nm_test_func();return(0);}
EOF

  if AC_TRY_EVAL(ac_compile); then
    # Now try to grab the symbols.
    nlist=conftest.nm
    if AC_TRY_EVAL(NM conftest.$ac_objext \| $lt_cv_sys_global_symbol_pipe \> $nlist) && test -s "$nlist"; then
      # Try sorting and uniquifying the output.
      if sort "$nlist" | uniq > "$nlist"T; then
	mv -f "$nlist"T "$nlist"
      else
	rm -f "$nlist"T
      fi

      # Make sure that we snagged all the symbols we need.
      if egrep ' nm_test_var$' "$nlist" >/dev/null; then
	if egrep ' nm_test_func$' "$nlist" >/dev/null; then
	  cat <<EOF > conftest.$ac_ext
#ifdef __cplusplus
extern "C" {
#endif

EOF
	  # Now generate the symbol file.
	  eval "$lt_cv_global_symbol_to_cdecl"' < "$nlist" >> conftest.$ac_ext'

	  cat <<EOF >> conftest.$ac_ext
#if defined (__STDC__) && __STDC__
# define lt_ptr void *
#else
# define lt_ptr char *
# define const
#endif

/* The mapping between symbol names and symbols. */
const struct {
  const char *name;
  lt_ptr address;
}
lt_preloaded_symbols[[]] =
{
EOF
	  sed "s/^$symcode$symcode* \(.*\) \(.*\)$/  {\"\2\", (lt_ptr) \&\2},/" < "$nlist" >> conftest.$ac_ext
	  cat <<\EOF >> conftest.$ac_ext
  {0, (lt_ptr) 0}
};

#ifdef __cplusplus
}
#endif
EOF
	  # Now try linking the two files.
	  mv conftest.$ac_objext conftstm.$ac_objext
	  save_LIBS="$LIBS"
	  save_CFLAGS="$CFLAGS"
	  LIBS="conftstm.$ac_objext"
	  CFLAGS="$CFLAGS$no_builtin_flag"
	  if AC_TRY_EVAL(ac_link) && test -s conftest$ac_exeext; then
	    pipe_works=yes
	  fi
	  LIBS="$save_LIBS"
	  CFLAGS="$save_CFLAGS"
	else
	  echo "cannot find nm_test_func in $nlist" >&AC_FD_CC
	fi
      else
	echo "cannot find nm_test_var in $nlist" >&AC_FD_CC
      fi
    else
      echo "cannot run $lt_cv_sys_global_symbol_pipe" >&AC_FD_CC
    fi
  else
    echo "$progname: failed program was:" >&AC_FD_CC
    cat conftest.$ac_ext >&5
  fi
  rm -f conftest* conftst*

  # Do not use the global_symbol_pipe unless it works.
  if test "$pipe_works" = yes; then
    break
  else
    lt_cv_sys_global_symbol_pipe=
  fi
done
])
global_symbol_pipe="$lt_cv_sys_global_symbol_pipe"
if test -z "$lt_cv_sys_global_symbol_pipe"; then
  global_symbol_to_cdecl=
  global_symbol_to_c_name_address=
else
  global_symbol_to_cdecl="$lt_cv_global_symbol_to_cdecl"
  global_symbol_to_c_name_address="$lt_cv_global_symbol_to_c_name_address"
fi
if test -z "$global_symbol_pipe$global_symbol_to_cdec$global_symbol_to_c_name_address";
then
  AC_MSG_RESULT(failed)
else
  AC_MSG_RESULT(ok)
fi
]) # AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE

# _LT_AC_LIBTOOL_SYS_PATH_SEPARATOR
# ---------------------------------
AC_DEFUN([_LT_AC_LIBTOOL_SYS_PATH_SEPARATOR],
[# Find the correct PATH separator.  Usually this is `:', but
# DJGPP uses `;' like DOS.
if test "X${PATH_SEPARATOR+set}" != Xset; then
  UNAME=${UNAME-`uname 2>/dev/null`}
  case X$UNAME in
    *-DOS) lt_cv_sys_path_separator=';' ;;
    *)     lt_cv_sys_path_separator=':' ;;
  esac
  PATH_SEPARATOR=$lt_cv_sys_path_separator
fi
])# _LT_AC_LIBTOOL_SYS_PATH_SEPARATOR

# _LT_AC_PROG_ECHO_BACKSLASH
# --------------------------
# Add some code to the start of the generated configure script which
# will find an echo command which doesn't interpret backslashes.
AC_DEFUN([_LT_AC_PROG_ECHO_BACKSLASH],
[ifdef([AC_DIVERSION_NOTICE], [AC_DIVERT_PUSH(AC_DIVERSION_NOTICE)],
			      [AC_DIVERT_PUSH(NOTICE)])
_LT_AC_LIBTOOL_SYS_PATH_SEPARATOR

# Check that we are running under the correct shell.
SHELL=${CONFIG_SHELL-/bin/sh}

case X$ECHO in
X*--fallback-echo)
  # Remove one level of quotation (which was required for Make).
  ECHO=`echo "$ECHO" | sed 's,\\\\\[$]\\[$]0,'[$]0','`
  ;;
esac

echo=${ECHO-echo}
if test "X[$]1" = X--no-reexec; then
  # Discard the --no-reexec flag, and continue.
  shift
elif test "X[$]1" = X--fallback-echo; then
  # Avoid inline document here, it may be left over
  :
elif test "X`($echo '\t') 2>/dev/null`" = 'X\t'; then
  # Yippee, $echo works!
  :
else
  # Restart under the correct shell.
  exec $SHELL "[$]0" --no-reexec ${1+"[$]@"}
fi

if test "X[$]1" = X--fallback-echo; then
  # used as fallback echo
  shift
  cat <<EOF
$*
EOF
  exit 0
fi

# The HP-UX ksh and POSIX shell print the target directory to stdout
# if CDPATH is set.
if test "X${CDPATH+set}" = Xset; then CDPATH=:; export CDPATH; fi

if test -z "$ECHO"; then
if test "X${echo_test_string+set}" != Xset; then
# find a string as large as possible, as long as the shell can cope with it
  for cmd in 'sed 50q "[$]0"' 'sed 20q "[$]0"' 'sed 10q "[$]0"' 'sed 2q "[$]0"' 'echo test'; do
    # expected sizes: less than 2Kb, 1Kb, 512 bytes, 16 bytes, ...
    if (echo_test_string="`eval $cmd`") 2>/dev/null &&
       echo_test_string="`eval $cmd`" &&
       (test "X$echo_test_string" = "X$echo_test_string") 2>/dev/null
    then
      break
    fi
  done
fi

if test "X`($echo '\t') 2>/dev/null`" = 'X\t' &&
   echo_testing_string=`($echo "$echo_test_string") 2>/dev/null` &&
   test "X$echo_testing_string" = "X$echo_test_string"; then
  :
else
  # The Solaris, AIX, and Digital Unix default echo programs unquote
  # backslashes.  This makes it impossible to quote backslashes using
  #   echo "$something" | sed 's/\\/\\\\/g'
  #
  # So, first we look for a working echo in the user's PATH.

  IFS="${IFS= 	}"; save_ifs="$IFS"; IFS=$PATH_SEPARATOR
  for dir in $PATH /usr/ucb; do
    if (test -f $dir/echo || test -f $dir/echo$ac_exeext) &&
       test "X`($dir/echo '\t') 2>/dev/null`" = 'X\t' &&
       echo_testing_string=`($dir/echo "$echo_test_string") 2>/dev/null` &&
       test "X$echo_testing_string" = "X$echo_test_string"; then
      echo="$dir/echo"
      break
    fi
  done
  IFS="$save_ifs"

  if test "X$echo" = Xecho; then
    # We didn't find a better echo, so look for alternatives.
    if test "X`(print -r '\t') 2>/dev/null`" = 'X\t' &&
       echo_testing_string=`(print -r "$echo_test_string") 2>/dev/null` &&
       test "X$echo_testing_string" = "X$echo_test_string"; then
      # This shell has a builtin print -r that does the trick.
      echo='print -r'
    elif (test -f /bin/ksh || test -f /bin/ksh$ac_exeext) &&
	 test "X$CONFIG_SHELL" != X/bin/ksh; then
      # If we have ksh, try running configure again with it.
      ORIGINAL_CONFIG_SHELL=${CONFIG_SHELL-/bin/sh}
      export ORIGINAL_CONFIG_SHELL
      CONFIG_SHELL=/bin/ksh
      export CONFIG_SHELL
      exec $CONFIG_SHELL "[$]0" --no-reexec ${1+"[$]@"}
    else
      # Try using printf.
      echo='printf %s\n'
      if test "X`($echo '\t') 2>/dev/null`" = 'X\t' &&
	 echo_testing_string=`($echo "$echo_test_string") 2>/dev/null` &&
	 test "X$echo_testing_string" = "X$echo_test_string"; then
	# Cool, printf works
	:
      elif echo_testing_string=`($ORIGINAL_CONFIG_SHELL "[$]0" --fallback-echo '\t') 2>/dev/null` &&
	   test "X$echo_testing_string" = 'X\t' &&
	   echo_testing_string=`($ORIGINAL_CONFIG_SHELL "[$]0" --fallback-echo "$echo_test_string") 2>/dev/null` &&
	   test "X$echo_testing_string" = "X$echo_test_string"; then
	CONFIG_SHELL=$ORIGINAL_CONFIG_SHELL
	export CONFIG_SHELL
	SHELL="$CONFIG_SHELL"
	export SHELL
	echo="$CONFIG_SHELL [$]0 --fallback-echo"
      elif echo_testing_string=`($CONFIG_SHELL "[$]0" --fallback-echo '\t') 2>/dev/null` &&
	   test "X$echo_testing_string" = 'X\t' &&
	   echo_testing_string=`($CONFIG_SHELL "[$]0" --fallback-echo "$echo_test_string") 2>/dev/null` &&
	   test "X$echo_testing_string" = "X$echo_test_string"; then
	echo="$CONFIG_SHELL [$]0 --fallback-echo"
      else
	# maybe with a smaller string...
	prev=:

	for cmd in 'echo test' 'sed 2q "[$]0"' 'sed 10q "[$]0"' 'sed 20q "[$]0"' 'sed 50q "[$]0"'; do
	  if (test "X$echo_test_string" = "X`eval $cmd`") 2>/dev/null
	  then
	    break
	  fi
	  prev="$cmd"
	done

	if test "$prev" != 'sed 50q "[$]0"'; then
	  echo_test_string=`eval $prev`
	  export echo_test_string
	  exec ${ORIGINAL_CONFIG_SHELL-${CONFIG_SHELL-/bin/sh}} "[$]0" ${1+"[$]@"}
	else
	  # Oops.  We lost completely, so just stick with echo.
	  echo=echo
	fi
      fi
    fi
  fi
fi
fi

# Copy echo and quote the copy suitably for passing to libtool from
# the Makefile, instead of quoting the original, which is used later.
ECHO=$echo
if test "X$ECHO" = "X$CONFIG_SHELL [$]0 --fallback-echo"; then
   ECHO="$CONFIG_SHELL \\\$\[$]0 --fallback-echo"
fi

AC_SUBST(ECHO)
AC_DIVERT_POP
])# _LT_AC_PROG_ECHO_BACKSLASH

# _LT_AC_TRY_DLOPEN_SELF (ACTION-IF-TRUE, ACTION-IF-TRUE-W-USCORE,
#                           ACTION-IF-FALSE, ACTION-IF-CROSS-COMPILING)
# ------------------------------------------------------------------
AC_DEFUN([_LT_AC_TRY_DLOPEN_SELF],
[if test "$cross_compiling" = yes; then :
  [$4]
else
  AC_REQUIRE([_LT_AC_CHECK_DLFCN])dnl
  lt_dlunknown=0; lt_dlno_uscore=1; lt_dlneed_uscore=2
  lt_status=$lt_dlunknown
  cat > conftest.$ac_ext <<EOF
[#line __oline__ "configure"
#include "confdefs.h"

#if HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#include <stdio.h>

#ifdef RTLD_GLOBAL
#  define LT_DLGLOBAL		RTLD_GLOBAL
#else
#  ifdef DL_GLOBAL
#    define LT_DLGLOBAL		DL_GLOBAL
#  else
#    define LT_DLGLOBAL		0
#  endif
#endif

/* We may have to define LT_DLLAZY_OR_NOW in the command line if we
   find out it does not work in some platform. */
#ifndef LT_DLLAZY_OR_NOW
#  ifdef RTLD_LAZY
#    define LT_DLLAZY_OR_NOW		RTLD_LAZY
#  else
#    ifdef DL_LAZY
#      define LT_DLLAZY_OR_NOW		DL_LAZY
#    else
#      ifdef RTLD_NOW
#        define LT_DLLAZY_OR_NOW	RTLD_NOW
#      else
#        ifdef DL_NOW
#          define LT_DLLAZY_OR_NOW	DL_NOW
#        else
#          define LT_DLLAZY_OR_NOW	0
#        endif
#      endif
#    endif
#  endif
#endif

#ifdef __cplusplus
extern "C" void exit (int);
#endif

void fnord() { int i=42;}
int main ()
{
  void *self = dlopen (0, LT_DLGLOBAL|LT_DLLAZY_OR_NOW);
  int status = $lt_dlunknown;

  if (self)
    {
      if (dlsym (self,"fnord"))       status = $lt_dlno_uscore;
      else if (dlsym( self,"_fnord")) status = $lt_dlneed_uscore;
      /* dlclose (self); */
    }

    exit (status);
}]
EOF
  if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext} 2>/dev/null; then
    (./conftest; exit; ) 2>/dev/null
    lt_status=$?
    case x$lt_status in
      x$lt_dlno_uscore) $1 ;;
      x$lt_dlneed_uscore) $2 ;;
      x$lt_unknown|x*) $3 ;;
    esac
  else :
    # compilation failed
    $3
  fi
fi
rm -fr conftest*
])# _LT_AC_TRY_DLOPEN_SELF

# AC_LIBTOOL_DLOPEN_SELF
# -------------------
AC_DEFUN([AC_LIBTOOL_DLOPEN_SELF],
[if test "x$enable_dlopen" != xyes; then
  enable_dlopen=unknown
  enable_dlopen_self=unknown
  enable_dlopen_self_static=unknown
else
  lt_cv_dlopen=no
  lt_cv_dlopen_libs=

  case $host_os in
  beos*)
    lt_cv_dlopen="load_add_on"
    lt_cv_dlopen_libs=
    lt_cv_dlopen_self=yes
    ;;

  cygwin* | mingw* | pw32*)
    lt_cv_dlopen="LoadLibrary"
    lt_cv_dlopen_libs=
   ;;

  *)
    AC_CHECK_FUNC([shl_load],
          [lt_cv_dlopen="shl_load"],
      [AC_CHECK_LIB([dld], [shl_load],
            [lt_cv_dlopen="shl_load" lt_cv_dlopen_libs="-dld"],
	[AC_CHECK_FUNC([dlopen],
	      [lt_cv_dlopen="dlopen"],
	  [AC_CHECK_LIB([dl], [dlopen],
	        [lt_cv_dlopen="dlopen" lt_cv_dlopen_libs="-ldl"],
	    [AC_CHECK_LIB([svld], [dlopen],
	          [lt_cv_dlopen="dlopen" lt_cv_dlopen_libs="-lsvld"],
	      [AC_CHECK_LIB([dld], [dld_link],
	            [lt_cv_dlopen="dld_link" lt_cv_dlopen_libs="-dld"])
	      ])
	    ])
	  ])
	])
      ])
    ;;
  esac

  if test "x$lt_cv_dlopen" != xno; then
    enable_dlopen=yes
  else
    enable_dlopen=no
  fi

  case $lt_cv_dlopen in
  dlopen)
    save_CPPFLAGS="$CPPFLAGS"
    AC_REQUIRE([_LT_AC_CHECK_DLFCN])dnl
    test "x$ac_cv_header_dlfcn_h" = xyes && CPPFLAGS="$CPPFLAGS -DHAVE_DLFCN_H"

    save_LDFLAGS="$LDFLAGS"
    eval LDFLAGS=\"\$LDFLAGS $export_dynamic_flag_spec\"

    save_LIBS="$LIBS"
    LIBS="$lt_cv_dlopen_libs $LIBS"

    AC_CACHE_CHECK([whether a program can dlopen itself],
	  lt_cv_dlopen_self, [dnl
	  _LT_AC_TRY_DLOPEN_SELF(
	    lt_cv_dlopen_self=yes, lt_cv_dlopen_self=yes,
	    lt_cv_dlopen_self=no, lt_cv_dlopen_self=cross)
    ])

    if test "x$lt_cv_dlopen_self" = xyes; then
      LDFLAGS="$LDFLAGS $link_static_flag"
      AC_CACHE_CHECK([whether a statically linked program can dlopen itself],
    	  lt_cv_dlopen_self_static, [dnl
	  _LT_AC_TRY_DLOPEN_SELF(
	    lt_cv_dlopen_self_static=yes, lt_cv_dlopen_self_static=yes,
	    lt_cv_dlopen_self_static=no,  lt_cv_dlopen_self_static=cross)
      ])
    fi

    CPPFLAGS="$save_CPPFLAGS"
    LDFLAGS="$save_LDFLAGS"
    LIBS="$save_LIBS"
    ;;
  esac

  case $lt_cv_dlopen_self in
  yes|no) enable_dlopen_self=$lt_cv_dlopen_self ;;
  *) enable_dlopen_self=unknown ;;
  esac

  case $lt_cv_dlopen_self_static in
  yes|no) enable_dlopen_self_static=$lt_cv_dlopen_self_static ;;
  *) enable_dlopen_self_static=unknown ;;
  esac
fi
])# AC_LIBTOOL_DLOPEN_SELF

AC_DEFUN([_LT_AC_LTCONFIG_HACK],
[AC_REQUIRE([AC_LIBTOOL_SYS_GLOBAL_SYMBOL_PIPE])dnl
# Sed substitution that helps us do robust quoting.  It backslashifies
# metacharacters that are still active within double-quoted strings.
Xsed='sed -e s/^X//'
sed_quote_subst='s/\([[\\"\\`$\\\\]]\)/\\\1/g'

# Same as above, but do not quote variable references.
double_quote_subst='s/\([[\\"\\`\\\\]]\)/\\\1/g'

# Sed substitution to delay expansion of an escaped shell variable in a
# double_quote_subst'ed string.
delay_variable_subst='s/\\\\\\\\\\\$/\\\\\\$/g'

# Constants:
rm="rm -f"

# Global variables:
default_ofile=libtool
can_build_shared=yes

# All known linkers require a `.a' archive for static linking (except M$VC,
# which needs '.lib').
libext=a
ltmain="$ac_aux_dir/ltmain.sh"
ofile="$default_ofile"
with_gnu_ld="$lt_cv_prog_gnu_ld"
need_locks="$enable_libtool_lock"

old_CC="$CC"
old_CFLAGS="$CFLAGS"

# Set sane defaults for various variables
test -z "$AR" && AR=ar
test -z "$AR_FLAGS" && AR_FLAGS=cru
test -z "$AS" && AS=as
test -z "$CC" && CC=cc
test -z "$DLLTOOL" && DLLTOOL=dlltool
test -z "$LD" && LD=ld
test -z "$LN_S" && LN_S="ln -s"
test -z "$MAGIC_CMD" && MAGIC_CMD=file
test -z "$NM" && NM=nm
test -z "$OBJDUMP" && OBJDUMP=objdump
test -z "$RANLIB" && RANLIB=:
test -z "$STRIP" && STRIP=:
test -z "$ac_objext" && ac_objext=o

if test x"$host" != x"$build"; then
  ac_tool_prefix=${host_alias}-
else
  ac_tool_prefix=
fi

# Transform linux* to *-*-linux-gnu*, to support old configure scripts.
case $host_os in
linux-gnu*) ;;
linux*) host=`echo $host | sed 's/^\(.*-.*-linux\)\(.*\)$/\1-gnu\2/'`
esac

case $host_os in
aix3*)
  # AIX sometimes has problems with the GCC collect2 program.  For some
  # reason, if we set the COLLECT_NAMES environment variable, the problems
  # vanish in a puff of smoke.
  if test "X${COLLECT_NAMES+set}" != Xset; then
    COLLECT_NAMES=
    export COLLECT_NAMES
  fi
  ;;
esac

# Determine commands to create old-style static archives.
old_archive_cmds='$AR $AR_FLAGS $oldlib$oldobjs$old_deplibs'
old_postinstall_cmds='chmod 644 $oldlib'
old_postuninstall_cmds=

if test -n "$RANLIB"; then
  case $host_os in
  openbsd*)
    old_postinstall_cmds="\$RANLIB -t \$oldlib~$old_postinstall_cmds"
    ;;
  *)
    old_postinstall_cmds="\$RANLIB \$oldlib~$old_postinstall_cmds"
    ;;
  esac
  old_archive_cmds="$old_archive_cmds~\$RANLIB \$oldlib"
fi

# Allow CC to be a program name with arguments.
set dummy $CC
compiler="[$]2"

AC_MSG_CHECKING([for objdir])
rm -f .libs 2>/dev/null
mkdir .libs 2>/dev/null
if test -d .libs; then
  objdir=.libs
else
  # MS-DOS does not allow filenames that begin with a dot.
  objdir=_libs
fi
rmdir .libs 2>/dev/null
AC_MSG_RESULT($objdir)


AC_ARG_WITH(pic,
[  --with-pic              try to use only PIC/non-PIC objects [default=use both]],
pic_mode="$withval", pic_mode=default)
test -z "$pic_mode" && pic_mode=default

# We assume here that the value for lt_cv_prog_cc_pic will not be cached
# in isolation, and that seeing it set (from the cache) indicates that
# the associated values are set (in the cache) correctly too.
AC_MSG_CHECKING([for $compiler option to produce PIC])
AC_CACHE_VAL(lt_cv_prog_cc_pic,
[ lt_cv_prog_cc_pic=
  lt_cv_prog_cc_shlib=
  lt_cv_prog_cc_wl=
  lt_cv_prog_cc_static=
  lt_cv_prog_cc_no_builtin=
  lt_cv_prog_cc_can_build_shared=$can_build_shared

  if test "$GCC" = yes; then
    lt_cv_prog_cc_wl='-Wl,'
    lt_cv_prog_cc_static='-static'

    case $host_os in
    aix*)
      # Below there is a dirty hack to force normal static linking with -ldl
      # The problem is because libdl dynamically linked with both libc and
      # libC (AIX C++ library), which obviously doesn't included in libraries
      # list by gcc. This cause undefined symbols with -static flags.
      # This hack allows C programs to be linked with "-static -ldl", but
      # not sure about C++ programs.
      lt_cv_prog_cc_static="$lt_cv_prog_cc_static ${lt_cv_prog_cc_wl}-lC"
      ;;
    amigaos*)
      # FIXME: we need at least 68020 code to build shared libraries, but
      # adding the `-m68020' flag to GCC prevents building anything better,
      # like `-m68040'.
      lt_cv_prog_cc_pic='-m68020 -resident32 -malways-restore-a4'
      ;;
    beos* | irix5* | irix6* | nonstopux* | osf3* | osf4* | osf5*)
      # PIC is the default for these OSes.
      ;;
    darwin* | rhapsody*)
      # PIC is the default on this platform
      # Common symbols not allowed in MH_DYLIB files
      lt_cv_prog_cc_pic='-fno-common'
      ;;
    cygwin* | mingw* | pw32* | os2*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      lt_cv_prog_cc_pic='-DDLL_EXPORT'
      ;;
    sysv4*MP*)
      if test -d /usr/nec; then
	 lt_cv_prog_cc_pic=-Kconform_pic
      fi
      ;;
    *)
      lt_cv_prog_cc_pic='-fPIC'
      ;;
    esac
  else
    # PORTME Check for PIC flags for the system compiler.
    case $host_os in
    aix3* | aix4* | aix5*)
      lt_cv_prog_cc_wl='-Wl,'
      # All AIX code is PIC.
      if test "$host_cpu" = ia64; then
	# AIX 5 now supports IA64 processor
	lt_cv_prog_cc_static='-Bstatic'
      else
	lt_cv_prog_cc_static='-bnso -bI:/lib/syscalls.exp'
      fi
      ;;

    hpux9* | hpux10* | hpux11*)
      # Is there a better lt_cv_prog_cc_static that works with the bundled CC?
      lt_cv_prog_cc_wl='-Wl,'
      lt_cv_prog_cc_static="${lt_cv_prog_cc_wl}-a ${lt_cv_prog_cc_wl}archive"
      lt_cv_prog_cc_pic='+Z'
      ;;

    irix5* | irix6* | nonstopux*)
      lt_cv_prog_cc_wl='-Wl,'
      lt_cv_prog_cc_static='-non_shared'
      # PIC (with -KPIC) is the default.
      ;;

    cygwin* | mingw* | pw32* | os2*)
      # This hack is so that the source file can tell whether it is being
      # built for inclusion in a dll (and should export symbols for example).
      lt_cv_prog_cc_pic='-DDLL_EXPORT'
      ;;

    newsos6)
      lt_cv_prog_cc_pic='-KPIC'
      lt_cv_prog_cc_static='-Bstatic'
      ;;

    osf3* | osf4* | osf5*)
      # All OSF/1 code is PIC.
      lt_cv_prog_cc_wl='-Wl,'
      lt_cv_prog_cc_static='-non_shared'
      ;;

    sco3.2v5*)
      lt_cv_prog_cc_pic='-Kpic'
      lt_cv_prog_cc_static='-dn'
      lt_cv_prog_cc_shlib='-belf'
      ;;

    solaris*)
      lt_cv_prog_cc_pic='-KPIC'
      lt_cv_prog_cc_static='-Bstatic'
      lt_cv_prog_cc_wl='-Wl,'
      ;;

    sunos4*)
      lt_cv_prog_cc_pic='-PIC'
      lt_cv_prog_cc_static='-Bstatic'
      lt_cv_prog_cc_wl='-Qoption ld '
      ;;

    sysv4 | sysv4.2uw2* | sysv4.3* | sysv5*)
      lt_cv_prog_cc_pic='-KPIC'
      lt_cv_prog_cc_static='-Bstatic'
      lt_cv_prog_cc_wl='-Wl,'
      ;;

    uts4*)
      lt_cv_prog_cc_pic='-pic'
      lt_cv_prog_cc_static='-Bstatic'
      ;;

    sysv4*MP*)
      if test -d /usr/nec ;then
	lt_cv_prog_cc_pic='-Kconform_pic'
	lt_cv_prog_cc_static='-Bstatic'
      fi
      ;;

    *)
      lt_cv_prog_cc_can_build_shared=no
      ;;
    esac
  fi
])
if test -z "$lt_cv_prog_cc_pic"; then
  AC_MSG_RESULT([none])
else
  AC_MSG_RESULT([$lt_cv_prog_cc_pic])

  # Check to make sure the pic_flag actually works.
  AC_MSG_CHECKING([if $compiler PIC flag $lt_cv_prog_cc_pic works])
  AC_CACHE_VAL(lt_cv_prog_cc_pic_works, [dnl
    save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS $lt_cv_prog_cc_pic -DPIC"
    AC_TRY_COMPILE([], [], [dnl
      case $host_os in
      hpux9* | hpux10* | hpux11*)
	# On HP-UX, both CC and GCC only warn that PIC is supported... then
	# they create non-PIC objects.  So, if there were any warnings, we
	# assume that PIC is not supported.
	if test -s conftest.err; then
	  lt_cv_prog_cc_pic_works=no
	else
	  lt_cv_prog_cc_pic_works=yes
	fi
	;;
      *)
	lt_cv_prog_cc_pic_works=yes
	;;
      esac
    ], [dnl
      lt_cv_prog_cc_pic_works=no
    ])
    CFLAGS="$save_CFLAGS"
  ])

  if test "X$lt_cv_prog_cc_pic_works" = Xno; then
    lt_cv_prog_cc_pic=
    lt_cv_prog_cc_can_build_shared=no
  else
    lt_cv_prog_cc_pic=" $lt_cv_prog_cc_pic"
  fi

  AC_MSG_RESULT([$lt_cv_prog_cc_pic_works])
fi

# Check for any special shared library compilation flags.
if test -n "$lt_cv_prog_cc_shlib"; then
  AC_MSG_WARN([\`$CC' requires \`$lt_cv_prog_cc_shlib' to build shared libraries])
  if echo "$old_CC $old_CFLAGS " | egrep -e "[[ 	]]$lt_cv_prog_cc_shlib[[ 	]]" >/dev/null; then :
  else
   AC_MSG_WARN([add \`$lt_cv_prog_cc_shlib' to the CC or CFLAGS env variable and reconfigure])
    lt_cv_prog_cc_can_build_shared=no
  fi
fi

AC_MSG_CHECKING([if $compiler static flag $lt_cv_prog_cc_static works])
AC_CACHE_VAL([lt_cv_prog_cc_static_works], [dnl
  lt_cv_prog_cc_static_works=no
  save_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS $lt_cv_prog_cc_static"
  AC_TRY_LINK([], [], [lt_cv_prog_cc_static_works=yes])
  LDFLAGS="$save_LDFLAGS"
])

# Belt *and* braces to stop my trousers falling down:
test "X$lt_cv_prog_cc_static_works" = Xno && lt_cv_prog_cc_static=
AC_MSG_RESULT([$lt_cv_prog_cc_static_works])

pic_flag="$lt_cv_prog_cc_pic"
special_shlib_compile_flags="$lt_cv_prog_cc_shlib"
wl="$lt_cv_prog_cc_wl"
link_static_flag="$lt_cv_prog_cc_static"
no_builtin_flag="$lt_cv_prog_cc_no_builtin"
can_build_shared="$lt_cv_prog_cc_can_build_shared"


# Check to see if options -o and -c are simultaneously supported by compiler
AC_MSG_CHECKING([if $compiler supports -c -o file.$ac_objext])
AC_CACHE_VAL([lt_cv_compiler_c_o], [
$rm -r conftest 2>/dev/null
mkdir conftest
cd conftest
echo "int some_variable = 0;" > conftest.$ac_ext
mkdir out
# According to Tom Tromey, Ian Lance Taylor reported there are C compilers
# that will create temporary files in the current directory regardless of
# the output directory.  Thus, making CWD read-only will cause this test
# to fail, enabling locking or at least warning the user not to do parallel
# builds.
chmod -w .
save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -o out/conftest2.$ac_objext"
compiler_c_o=no
if { (eval echo configure:__oline__: \"$ac_compile\") 1>&5; (eval $ac_compile) 2>out/conftest.err; } && test -s out/conftest2.$ac_objext; then
  # The compiler can only warn and ignore the option if not recognized
  # So say no if there are warnings
  if test -s out/conftest.err; then
    lt_cv_compiler_c_o=no
  else
    lt_cv_compiler_c_o=yes
  fi
else
  # Append any errors to the config.log.
  cat out/conftest.err 1>&AC_FD_CC
  lt_cv_compiler_c_o=no
fi
CFLAGS="$save_CFLAGS"
chmod u+w .
$rm conftest* out/*
rmdir out
cd ..
rmdir conftest
$rm -r conftest 2>/dev/null
])
compiler_c_o=$lt_cv_compiler_c_o
AC_MSG_RESULT([$compiler_c_o])

if test x"$compiler_c_o" = x"yes"; then
  # Check to see if we can write to a .lo
  AC_MSG_CHECKING([if $compiler supports -c -o file.lo])
  AC_CACHE_VAL([lt_cv_compiler_o_lo], [
  lt_cv_compiler_o_lo=no
  save_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -c -o conftest.lo"
  save_objext="$ac_objext"
  ac_objext=lo
  AC_TRY_COMPILE([], [int some_variable = 0;], [dnl
    # The compiler can only warn and ignore the option if not recognized
    # So say no if there are warnings
    if test -s conftest.err; then
      lt_cv_compiler_o_lo=no
    else
      lt_cv_compiler_o_lo=yes
    fi
  ])
  ac_objext="$save_objext"
  CFLAGS="$save_CFLAGS"
  ])
  compiler_o_lo=$lt_cv_compiler_o_lo
  AC_MSG_RESULT([$compiler_o_lo])
else
  compiler_o_lo=no
fi

# Check to see if we can do hard links to lock some files if needed
hard_links="nottested"
if test "$compiler_c_o" = no && test "$need_locks" != no; then
  # do not overwrite the value of need_locks provided by the user
  AC_MSG_CHECKING([if we can lock with hard links])
  hard_links=yes
  $rm conftest*
  ln conftest.a conftest.b 2>/dev/null && hard_links=no
  touch conftest.a
  ln conftest.a conftest.b 2>&5 || hard_links=no
  ln conftest.a conftest.b 2>/dev/null && hard_links=no
  AC_MSG_RESULT([$hard_links])
  if test "$hard_links" = no; then
    AC_MSG_WARN([\`$CC' does not support \`-c -o', so \`make -j' may be unsafe])
    need_locks=warn
  fi
else
  need_locks=no
fi

if test "$GCC" = yes; then
  # Check to see if options -fno-rtti -fno-exceptions are supported by compiler
  AC_MSG_CHECKING([if $compiler supports -fno-rtti -fno-exceptions])
  echo "int some_variable = 0;" > conftest.$ac_ext
  save_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS -fno-rtti -fno-exceptions -c conftest.$ac_ext"
  compiler_rtti_exceptions=no
  AC_TRY_COMPILE([], [int some_variable = 0;], [dnl
    # The compiler can only warn and ignore the option if not recognized
    # So say no if there are warnings
    if test -s conftest.err; then
      compiler_rtti_exceptions=no
    else
      compiler_rtti_exceptions=yes
    fi
  ])
  CFLAGS="$save_CFLAGS"
  AC_MSG_RESULT([$compiler_rtti_exceptions])

  if test "$compiler_rtti_exceptions" = "yes"; then
    no_builtin_flag=' -fno-builtin -fno-rtti -fno-exceptions'
  else
    no_builtin_flag=' -fno-builtin'
  fi
fi

# See if the linker supports building shared libraries.
AC_MSG_CHECKING([whether the linker ($LD) supports shared libraries])

allow_undefined_flag=
no_undefined_flag=
need_lib_prefix=unknown
need_version=unknown
# when you set need_version to no, make sure it does not cause -set_version
# flags to be left without arguments
archive_cmds=
archive_expsym_cmds=
old_archive_from_new_cmds=
old_archive_from_expsyms_cmds=
export_dynamic_flag_spec=
whole_archive_flag_spec=
thread_safe_flag_spec=
hardcode_into_libs=no
hardcode_libdir_flag_spec=
hardcode_libdir_separator=
hardcode_direct=no
hardcode_minus_L=no
hardcode_shlibpath_var=unsupported
runpath_var=
link_all_deplibs=unknown
always_export_symbols=no
export_symbols_cmds='$NM $libobjs $convenience | $global_symbol_pipe | sed '\''s/.* //'\'' | sort | uniq > $export_symbols'
# include_expsyms should be a list of space-separated symbols to be *always*
# included in the symbol list
include_expsyms=
# exclude_expsyms can be an egrep regular expression of symbols to exclude
# it will be wrapped by ` (' and `)$', so one must not match beginning or
# end of line.  Example: `a|bc|.*d.*' will exclude the symbols `a' and `bc',
# as well as any symbol that contains `d'.
exclude_expsyms="_GLOBAL_OFFSET_TABLE_"
# Although _GLOBAL_OFFSET_TABLE_ is a valid symbol C name, most a.out
# platforms (ab)use it in PIC code, but their linkers get confused if
# the symbol is explicitly referenced.  Since portable code cannot
# rely on this symbol name, it's probably fine to never include it in
# preloaded symbol tables.
extract_expsyms_cmds=

case $host_os in
cygwin* | mingw* | pw32*)
  # FIXME: the MSVC++ port hasn't been tested in a loooong time
  # When not using gcc, we currently assume that we are using
  # Microsoft Visual C++.
  if test "$GCC" != yes; then
    with_gnu_ld=no
  fi
  ;;
openbsd*)
  with_gnu_ld=no
  ;;
esac

ld_shlibs=yes
if test "$with_gnu_ld" = yes; then
  # If archive_cmds runs LD, not CC, wlarc should be empty
  wlarc='${wl}'

  # See if GNU ld supports shared libraries.
  case $host_os in
  aix3* | aix4* | aix5*)
    # On AIX, the GNU linker is very broken
    # Note:Check GNU linker on AIX 5-IA64 when/if it becomes available.
    ld_shlibs=no
    cat <<EOF 1>&2

*** Warning: the GNU linker, at least up to release 2.9.1, is reported
*** to be unable to reliably create shared libraries on AIX.
*** Therefore, libtool is disabling shared libraries support.  If you
*** really care for shared libraries, you may want to modify your PATH
*** so that a non-GNU linker is found, and then restart.

EOF
    ;;

  amigaos*)
    archive_cmds='$rm $output_objdir/a2ixlibrary.data~$echo "#define NAME $libname" > $output_objdir/a2ixlibrary.data~$echo "#define LIBRARY_ID 1" >> $output_objdir/a2ixlibrary.data~$echo "#define VERSION $major" >> $output_objdir/a2ixlibrary.data~$echo "#define REVISION $revision" >> $output_objdir/a2ixlibrary.data~$AR $AR_FLAGS $lib $libobjs~$RANLIB $lib~(cd $output_objdir && a2ixlibrary -32)'
    hardcode_libdir_flag_spec='-L$libdir'
    hardcode_minus_L=yes

    # Samuel A. Falvo II <kc5tja@dolphin.openprojects.net> reports
    # that the semantics of dynamic libraries on AmigaOS, at least up
    # to version 4, is to share data among multiple programs linked
    # with the same dynamic library.  Since this doesn't match the
    # behavior of shared libraries on other platforms, we can use
    # them.
    ld_shlibs=no
    ;;

  beos*)
    if $LD --help 2>&1 | egrep ': supported targets:.* elf' > /dev/null; then
      allow_undefined_flag=unsupported
      # Joseph Beckenbach <jrb3@best.com> says some releases of gcc
      # support --undefined.  This deserves some investigation.  FIXME
      archive_cmds='$CC -nostart $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
    else
      ld_shlibs=no
    fi
    ;;

  cygwin* | mingw* | pw32*)
    # hardcode_libdir_flag_spec is actually meaningless, as there is
    # no search path for DLLs.
    hardcode_libdir_flag_spec='-L$libdir'
    allow_undefined_flag=unsupported
    always_export_symbols=yes

    extract_expsyms_cmds='test -f $output_objdir/impgen.c || \
      sed -e "/^# \/\* impgen\.c starts here \*\//,/^# \/\* impgen.c ends here \*\// { s/^# //;s/^# *$//; p; }" -e d < $''0 > $output_objdir/impgen.c~
      test -f $output_objdir/impgen.exe || (cd $output_objdir && \
      if test "x$HOST_CC" != "x" ; then $HOST_CC -o impgen impgen.c ; \
      else $CC -o impgen impgen.c ; fi)~
      $output_objdir/impgen $dir/$soroot > $output_objdir/$soname-def'

    old_archive_from_expsyms_cmds='$DLLTOOL --as=$AS --dllname $soname --def $output_objdir/$soname-def --output-lib $output_objdir/$newlib'

    # cygwin and mingw dlls have different entry points and sets of symbols
    # to exclude.
    # FIXME: what about values for MSVC?
    dll_entry=__cygwin_dll_entry@12
    dll_exclude_symbols=DllMain@12,_cygwin_dll_entry@12,_cygwin_noncygwin_dll_entry@12~
    case $host_os in
    mingw*)
      # mingw values
      dll_entry=_DllMainCRTStartup@12
      dll_exclude_symbols=DllMain@12,DllMainCRTStartup@12,DllEntryPoint@12~
      ;;
    esac

    # mingw and cygwin differ, and it's simplest to just exclude the union
    # of the two symbol sets.
    dll_exclude_symbols=DllMain@12,_cygwin_dll_entry@12,_cygwin_noncygwin_dll_entry@12,DllMainCRTStartup@12,DllEntryPoint@12

    # recent cygwin and mingw systems supply a stub DllMain which the user
    # can override, but on older systems we have to supply one (in ltdll.c)
    if test "x$lt_cv_need_dllmain" = "xyes"; then
      ltdll_obj='$output_objdir/$soname-ltdll.'"$ac_objext "
      ltdll_cmds='test -f $output_objdir/$soname-ltdll.c || sed -e "/^# \/\* ltdll\.c starts here \*\//,/^# \/\* ltdll.c ends here \*\// { s/^# //; p; }" -e d < $''0 > $output_objdir/$soname-ltdll.c~
	test -f $output_objdir/$soname-ltdll.$ac_objext || (cd $output_objdir && $CC -c $soname-ltdll.c)~'
    else
      ltdll_obj=
      ltdll_cmds=
    fi

    # Extract the symbol export list from an `--export-all' def file,
    # then regenerate the def file from the symbol export list, so that
    # the compiled dll only exports the symbol export list.
    # Be careful not to strip the DATA tag left be newer dlltools.
    export_symbols_cmds="$ltdll_cmds"'
      $DLLTOOL --export-all --exclude-symbols '$dll_exclude_symbols' --output-def $output_objdir/$soname-def '$ltdll_obj'$libobjs $convenience~
      sed -e "1,/EXPORTS/d" -e "s/ @ [[0-9]]*//" -e "s/ *;.*$//" < $output_objdir/$soname-def > $export_symbols'

    # If the export-symbols file already is a .def file (1st line
    # is EXPORTS), use it as is.
    # If DATA tags from a recent dlltool are present, honour them!
    archive_expsym_cmds='if test "x`sed 1q $export_symbols`" = xEXPORTS; then
	cp $export_symbols $output_objdir/$soname-def;
      else
	echo EXPORTS > $output_objdir/$soname-def;
	_lt_hint=1;
	cat $export_symbols | while read symbol; do
	 set dummy \$symbol;
	 case \[$]# in
	   2) echo "   \[$]2 @ \$_lt_hint ; " >> $output_objdir/$soname-def;;
	   4) echo "   \[$]2 \[$]3 \[$]4 ; " >> $output_objdir/$soname-def; _lt_hint=`expr \$_lt_hint - 1`;;
	   *) echo "     \[$]2 @ \$_lt_hint \[$]3 ; " >> $output_objdir/$soname-def;;
	 esac;
	 _lt_hint=`expr 1 + \$_lt_hint`;
	done;
      fi~
      '"$ltdll_cmds"'
      $CC -Wl,--base-file,$output_objdir/$soname-base '$lt_cv_cc_dll_switch' -Wl,-e,'$dll_entry' -o $output_objdir/$soname '$ltdll_obj'$libobjs $deplibs $compiler_flags~
      $DLLTOOL --as=$AS --dllname $soname --exclude-symbols '$dll_exclude_symbols' --def $output_objdir/$soname-def --base-file $output_objdir/$soname-base --output-exp $output_objdir/$soname-exp~
      $CC -Wl,--base-file,$output_objdir/$soname-base $output_objdir/$soname-exp '$lt_cv_cc_dll_switch' -Wl,-e,'$dll_entry' -o $output_objdir/$soname '$ltdll_obj'$libobjs $deplibs $compiler_flags~
      $DLLTOOL --as=$AS --dllname $soname --exclude-symbols '$dll_exclude_symbols' --def $output_objdir/$soname-def --base-file $output_objdir/$soname-base --output-exp $output_objdir/$soname-exp --output-lib $output_objdir/$libname.dll.a~
      $CC $output_objdir/$soname-exp '$lt_cv_cc_dll_switch' -Wl,-e,'$dll_entry' -o $output_objdir/$soname '$ltdll_obj'$libobjs $deplibs $compiler_flags'
    ;;

  netbsd*)
    if echo __ELF__ | $CC -E - | grep __ELF__ >/dev/null; then
      archive_cmds='$LD -Bshareable $libobjs $deplibs $linker_flags -o $lib'
      wlarc=
    else
      archive_cmds='$CC -shared -nodefaultlibs $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
      archive_expsym_cmds='$CC -shared -nodefaultlibs $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
    fi
    ;;

  solaris* | sysv5*)
    if $LD -v 2>&1 | egrep 'BFD 2\.8' > /dev/null; then
      ld_shlibs=no
      cat <<EOF 1>&2

*** Warning: The releases 2.8.* of the GNU linker cannot reliably
*** create shared libraries on Solaris systems.  Therefore, libtool
*** is disabling shared libraries support.  We urge you to upgrade GNU
*** binutils to release 2.9.1 or newer.  Another option is to modify
*** your PATH or compiler configuration so that the native linker is
*** used, and then restart.

EOF
    elif $LD --help 2>&1 | egrep ': supported targets:.* elf' > /dev/null; then
      archive_cmds='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
      archive_expsym_cmds='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
    else
      ld_shlibs=no
    fi
    ;;

  sunos4*)
    archive_cmds='$LD -assert pure-text -Bshareable -o $lib $libobjs $deplibs $linker_flags'
    wlarc=
    hardcode_direct=yes
    hardcode_shlibpath_var=no
    ;;

  *)
    if $LD --help 2>&1 | egrep ': supported targets:.* elf' > /dev/null; then
      archive_cmds='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname -o $lib'
      archive_expsym_cmds='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname $wl$soname ${wl}-retain-symbols-file $wl$export_symbols -o $lib'
    else
      ld_shlibs=no
    fi
    ;;
  esac

  if test "$ld_shlibs" = yes; then
    runpath_var=LD_RUN_PATH
    hardcode_libdir_flag_spec='${wl}--rpath ${wl}$libdir'
    export_dynamic_flag_spec='${wl}--export-dynamic'
    case $host_os in
    cygwin* | mingw* | pw32*)
      # dlltool doesn't understand --whole-archive et. al.
      whole_archive_flag_spec=
      ;;
    *)
      # ancient GNU ld didn't support --whole-archive et. al.
      if $LD --help 2>&1 | egrep 'no-whole-archive' > /dev/null; then
	whole_archive_flag_spec="$wlarc"'--whole-archive$convenience '"$wlarc"'--no-whole-archive'
      else
	whole_archive_flag_spec=
      fi
      ;;
    esac
  fi
else
  # PORTME fill in a description of your system's linker (not GNU ld)
  case $host_os in
  aix3*)
    allow_undefined_flag=unsupported
    always_export_symbols=yes
    archive_expsym_cmds='$LD -o $output_objdir/$soname $libobjs $deplibs $linker_flags -bE:$export_symbols -T512 -H512 -bM:SRE~$AR $AR_FLAGS $lib $output_objdir/$soname'
    # Note: this linker hardcodes the directories in LIBPATH if there
    # are no directories specified by -L.
    hardcode_minus_L=yes
    if test "$GCC" = yes && test -z "$link_static_flag"; then
      # Neither direct hardcoding nor static linking is supported with a
      # broken collect2.
      hardcode_direct=unsupported
    fi
    ;;

  aix4* | aix5*)
    if test "$host_cpu" = ia64; then
      # On IA64, the linker does run time linking by default, so we don't
      # have to do anything special.
      aix_use_runtimelinking=no
      exp_sym_flag='-Bexport'
      no_entry_flag=""
    else
      aix_use_runtimelinking=no

      # Test if we are trying to use run time linking or normal
      # AIX style linking. If -brtl is somewhere in LDFLAGS, we
      # need to do runtime linking.
      case $host_os in aix4.[[23]]|aix4.[[23]].*|aix5*)
	for ld_flag in $LDFLAGS; do
	  case $ld_flag in
	  *-brtl*)
	    aix_use_runtimelinking=yes
	    break
	  ;;
	  esac
	done
      esac

      exp_sym_flag='-bexport'
      no_entry_flag='-bnoentry'
    fi

    # When large executables or shared objects are built, AIX ld can
    # have problems creating the table of contents.  If linking a library
    # or program results in "error TOC overflow" add -mminimal-toc to
    # CXXFLAGS/CFLAGS for g++/gcc.  In the cases where that is not
    # enough to fix the problem, add -Wl,-bbigtoc to LDFLAGS.

    hardcode_direct=yes
    archive_cmds=''
    hardcode_libdir_separator=':'
    if test "$GCC" = yes; then
      case $host_os in aix4.[[012]]|aix4.[[012]].*)
	collect2name=`${CC} -print-prog-name=collect2`
	if test -f "$collect2name" && \
	  strings "$collect2name" | grep resolve_lib_name >/dev/null
	then
	  # We have reworked collect2
	  hardcode_direct=yes
	else
	  # We have old collect2
	  hardcode_direct=unsupported
	  # It fails to find uninstalled libraries when the uninstalled
	  # path is not listed in the libpath.  Setting hardcode_minus_L
	  # to unsupported forces relinking
	  hardcode_minus_L=yes
	  hardcode_libdir_flag_spec='-L$libdir'
	  hardcode_libdir_separator=
	fi
      esac

      shared_flag='-shared'
    else
      # not using gcc
      if test "$host_cpu" = ia64; then
	shared_flag='${wl}-G'
      else
	if test "$aix_use_runtimelinking" = yes; then
	  shared_flag='${wl}-G'
	else
	  shared_flag='${wl}-bM:SRE'
	fi
      fi
    fi

    # It seems that -bexpall can do strange things, so it is better to
    # generate a list of symbols to export.
    always_export_symbols=yes
    if test "$aix_use_runtimelinking" = yes; then
      # Warning - without using the other runtime loading flags (-brtl),
      # -berok will link without error, but may produce a broken library.
      allow_undefined_flag='-berok'
      hardcode_libdir_flag_spec='${wl}-blibpath:$libdir:/usr/lib:/lib'
      archive_expsym_cmds="\$CC"' -o $output_objdir/$soname $libobjs $deplibs $compiler_flags `if test "x${allow_undefined_flag}" != "x"; then echo "${wl}${allow_undefined_flag}"; else :; fi` '"\${wl}$no_entry_flag \${wl}$exp_sym_flag:\$export_symbols $shared_flag"
    else
      if test "$host_cpu" = ia64; then
	hardcode_libdir_flag_spec='${wl}-R $libdir:/usr/lib:/lib'
	allow_undefined_flag="-z nodefs"
	archive_expsym_cmds="\$CC $shared_flag"' -o $output_objdir/$soname ${wl}-h$soname $libobjs $deplibs $compiler_flags ${wl}${allow_undefined_flag} '"\${wl}$no_entry_flag \${wl}$exp_sym_flag:\$export_symbols"
      else
	hardcode_libdir_flag_spec='${wl}-bnolibpath ${wl}-blibpath:$libdir:/usr/lib:/lib'
	# Warning - without using the other run time loading flags,
	# -berok will link without error, but may produce a broken library.
	allow_undefined_flag='${wl}-berok'
	# This is a bit strange, but is similar to how AIX traditionally builds
	# it's shared libraries.
	archive_expsym_cmds="\$CC $shared_flag"' -o $output_objdir/$soname $libobjs $deplibs $compiler_flags ${allow_undefined_flag} '"\${wl}$no_entry_flag \${wl}$exp_sym_flag:\$export_symbols"' ~$AR -crlo $objdir/$libname$release.a $objdir/$soname'
      fi
    fi
    ;;

  amigaos*)
    archive_cmds='$rm $output_objdir/a2ixlibrary.data~$echo "#define NAME $libname" > $output_objdir/a2ixlibrary.data~$echo "#define LIBRARY_ID 1" >> $output_objdir/a2ixlibrary.data~$echo "#define VERSION $major" >> $output_objdir/a2ixlibrary.data~$echo "#define REVISION $revision" >> $output_objdir/a2ixlibrary.data~$AR $AR_FLAGS $lib $libobjs~$RANLIB $lib~(cd $output_objdir && a2ixlibrary -32)'
    hardcode_libdir_flag_spec='-L$libdir'
    hardcode_minus_L=yes
    # see comment about different semantics on the GNU ld section
    ld_shlibs=no
    ;;

  cygwin* | mingw* | pw32*)
    # When not using gcc, we currently assume that we are using
    # Microsoft Visual C++.
    # hardcode_libdir_flag_spec is actually meaningless, as there is
    # no search path for DLLs.
    hardcode_libdir_flag_spec=' '
    allow_undefined_flag=unsupported
    # Tell ltmain to make .lib files, not .a files.
    libext=lib
    # FIXME: Setting linknames here is a bad hack.
    archive_cmds='$CC -o $lib $libobjs $compiler_flags `echo "$deplibs" | sed -e '\''s/ -lc$//'\''` -link -dll~linknames='
    # The linker will automatically build a .lib file if we build a DLL.
    old_archive_from_new_cmds='true'
    # FIXME: Should let the user specify the lib program.
    old_archive_cmds='lib /OUT:$oldlib$oldobjs$old_deplibs'
    fix_srcfile_path='`cygpath -w "$srcfile"`'
    ;;

  darwin* | rhapsody*)
    case "$host_os" in
    rhapsody* | darwin1.[[012]])
      allow_undefined_flag='-undefined suppress'
      ;;
    *) # Darwin 1.3 on
      allow_undefined_flag='-flat_namespace -undefined suppress'
      ;;
    esac
    # FIXME: Relying on posixy $() will cause problems for
    #        cross-compilation, but unfortunately the echo tests do not
    #        yet detect zsh echo's removal of \ escapes.  Also zsh mangles
    #	     `"' quotes if we put them in here... so don't!
    archive_cmds='$CC -r -keep_private_externs -nostdlib -o ${lib}-master.o $libobjs && $CC $(test .$module = .yes && echo -bundle || echo -dynamiclib) $allow_undefined_flag -o $lib ${lib}-master.o $deplibs$linker_flags $(test .$module != .yes && echo -install_name $rpath/$soname $verstring)'
    # We need to add '_' to the symbols in $export_symbols first
    #archive_expsym_cmds="$archive_cmds"' && strip -s $export_symbols'
    hardcode_direct=yes
    hardcode_shlibpath_var=no
    whole_archive_flag_spec='-all_load $convenience'
    ;;

  freebsd1*)
    ld_shlibs=no
    ;;

  # FreeBSD 2.2.[012] allows us to include c++rt0.o to get C++ constructor
  # support.  Future versions do this automatically, but an explicit c++rt0.o
  # does not break anything, and helps significantly (at the cost of a little
  # extra space).
  freebsd2.2*)
    archive_cmds='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags /usr/lib/c++rt0.o'
    hardcode_libdir_flag_spec='-R$libdir'
    hardcode_direct=yes
    hardcode_shlibpath_var=no
    ;;

  # Unfortunately, older versions of FreeBSD 2 do not have this feature.
  freebsd2*)
    archive_cmds='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'
    hardcode_direct=yes
    hardcode_minus_L=yes
    hardcode_shlibpath_var=no
    ;;

  # FreeBSD 3 and greater uses gcc -shared to do shared libraries.
  freebsd*)
    archive_cmds='$CC -shared -o $lib $libobjs $deplibs $compiler_flags'
    hardcode_libdir_flag_spec='-R$libdir'
    hardcode_direct=yes
    hardcode_shlibpath_var=no
    ;;

  hpux9* | hpux10* | hpux11*)
    case $host_os in
    hpux9*) archive_cmds='$rm $output_objdir/$soname~$LD -b +b $install_libdir -o $output_objdir/$soname $libobjs $deplibs $linker_flags~test $output_objdir/$soname = $lib || mv $output_objdir/$soname $lib' ;;
    *) archive_cmds='$LD -b +h $soname +b $install_libdir -o $lib $libobjs $deplibs $linker_flags' ;;
    esac
    hardcode_libdir_flag_spec='${wl}+b ${wl}$libdir'
    hardcode_libdir_separator=:
    hardcode_direct=yes
    hardcode_minus_L=yes # Not in the search PATH, but as the default
			 # location of the library.
    export_dynamic_flag_spec='${wl}-E'
    ;;

  irix5* | irix6* | nonstopux*)
    if test "$GCC" = yes; then
      archive_cmds='$CC -shared $libobjs $deplibs $compiler_flags ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'
      hardcode_libdir_flag_spec='${wl}-rpath ${wl}$libdir'
    else
      archive_cmds='$LD -shared $libobjs $deplibs $linker_flags -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'
      hardcode_libdir_flag_spec='-rpath $libdir'
    fi
    hardcode_libdir_separator=:
    link_all_deplibs=yes
    ;;

  netbsd*)
    if echo __ELF__ | $CC -E - | grep __ELF__ >/dev/null; then
      archive_cmds='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'  # a.out
    else
      archive_cmds='$LD -shared -o $lib $libobjs $deplibs $linker_flags'      # ELF
    fi
    hardcode_libdir_flag_spec='-R$libdir'
    hardcode_direct=yes
    hardcode_shlibpath_var=no
    ;;

  newsos6)
    archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
    hardcode_direct=yes
    hardcode_libdir_flag_spec='${wl}-rpath ${wl}$libdir'
    hardcode_libdir_separator=:
    hardcode_shlibpath_var=no
    ;;

  openbsd*)
    hardcode_direct=yes
    hardcode_shlibpath_var=no
    if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`" || test "$host_os-$host_cpu" = "openbsd2.8-powerpc"; then
      archive_cmds='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
      hardcode_libdir_flag_spec='${wl}-rpath,$libdir'
      export_dynamic_flag_spec='${wl}-E'
    else
      case "$host_os" in
      openbsd[[01]].* | openbsd2.[[0-7]] | openbsd2.[[0-7]].*)
	archive_cmds='$LD -Bshareable -o $lib $libobjs $deplibs $linker_flags'
	hardcode_libdir_flag_spec='-R$libdir'
        ;;
      *)
        archive_cmds='$CC -shared $pic_flag -o $lib $libobjs $deplibs $compiler_flags'
        hardcode_libdir_flag_spec='${wl}-rpath,$libdir'
        ;;
      esac
    fi
    ;;

  os2*)
    hardcode_libdir_flag_spec='-L$libdir'
    hardcode_minus_L=yes
    allow_undefined_flag=unsupported
    archive_cmds='$echo "LIBRARY $libname INITINSTANCE" > $output_objdir/$libname.def~$echo "DESCRIPTION \"$libname\"" >> $output_objdir/$libname.def~$echo DATA >> $output_objdir/$libname.def~$echo " SINGLE NONSHARED" >> $output_objdir/$libname.def~$echo EXPORTS >> $output_objdir/$libname.def~emxexp $libobjs >> $output_objdir/$libname.def~$CC -Zdll -Zcrtdll -o $lib $libobjs $deplibs $compiler_flags $output_objdir/$libname.def'
    old_archive_from_new_cmds='emximp -o $output_objdir/$libname.a $output_objdir/$libname.def'
    ;;

  osf3*)
    if test "$GCC" = yes; then
      allow_undefined_flag=' ${wl}-expect_unresolved ${wl}\*'
      archive_cmds='$CC -shared${allow_undefined_flag} $libobjs $deplibs $compiler_flags ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'
    else
      allow_undefined_flag=' -expect_unresolved \*'
      archive_cmds='$LD -shared${allow_undefined_flag} $libobjs $deplibs $linker_flags -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'
    fi
    hardcode_libdir_flag_spec='${wl}-rpath ${wl}$libdir'
    hardcode_libdir_separator=:
    ;;

  osf4* | osf5*)	# as osf3* with the addition of -msym flag
    if test "$GCC" = yes; then
      allow_undefined_flag=' ${wl}-expect_unresolved ${wl}\*'
      archive_cmds='$CC -shared${allow_undefined_flag} $libobjs $deplibs $compiler_flags ${wl}-msym ${wl}-soname ${wl}$soname `test -n "$verstring" && echo ${wl}-set_version ${wl}$verstring` ${wl}-update_registry ${wl}${output_objdir}/so_locations -o $lib'
      hardcode_libdir_flag_spec='${wl}-rpath ${wl}$libdir'
    else
      allow_undefined_flag=' -expect_unresolved \*'
      archive_cmds='$LD -shared${allow_undefined_flag} $libobjs $deplibs $linker_flags -msym -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${output_objdir}/so_locations -o $lib'
      archive_expsym_cmds='for i in `cat $export_symbols`; do printf "-exported_symbol " >> $lib.exp; echo "\$i" >> $lib.exp; done; echo "-hidden">> $lib.exp~
      $LD -shared${allow_undefined_flag} -input $lib.exp $linker_flags $libobjs $deplibs -soname $soname `test -n "$verstring" && echo -set_version $verstring` -update_registry ${objdir}/so_locations -o $lib~$rm $lib.exp'

      #Both c and cxx compiler support -rpath directly
      hardcode_libdir_flag_spec='-rpath $libdir'
    fi
    hardcode_libdir_separator=:
    ;;

  sco3.2v5*)
    archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
    hardcode_shlibpath_var=no
    runpath_var=LD_RUN_PATH
    hardcode_runpath_var=yes
    export_dynamic_flag_spec='${wl}-Bexport'
    ;;

  solaris*)
    # gcc --version < 3.0 without binutils cannot create self contained
    # shared libraries reliably, requiring libgcc.a to resolve some of
    # the object symbols generated in some cases.  Libraries that use
    # assert need libgcc.a to resolve __eprintf, for example.  Linking
    # a copy of libgcc.a into every shared library to guarantee resolving
    # such symbols causes other problems:  According to Tim Van Holder
    # <tim.van.holder@pandora.be>, C++ libraries end up with a separate
    # (to the application) exception stack for one thing.
    no_undefined_flag=' -z defs'
    if test "$GCC" = yes; then
      case `$CC --version 2>/dev/null` in
      [[12]].*)
	cat <<EOF 1>&2

*** Warning: Releases of GCC earlier than version 3.0 cannot reliably
*** create self contained shared libraries on Solaris systems, without
*** introducing a dependency on libgcc.a.  Therefore, libtool is disabling
*** -no-undefined support, which will at least allow you to build shared
*** libraries.  However, you may find that when you link such libraries
*** into an application without using GCC, you have to manually add
*** \`gcc --print-libgcc-file-name\` to the link command.  We urge you to
*** upgrade to a newer version of GCC.  Another option is to rebuild your
*** current GCC to use the GNU linker from GNU binutils 2.9.1 or newer.

EOF
        no_undefined_flag=
	;;
      esac
    fi
    # $CC -shared without GNU ld will not create a library from C++
    # object files and a static libstdc++, better avoid it by now
    archive_cmds='$LD -G${allow_undefined_flag} -h $soname -o $lib $libobjs $deplibs $linker_flags'
    archive_expsym_cmds='$echo "{ global:" > $lib.exp~cat $export_symbols | sed -e "s/\(.*\)/\1;/" >> $lib.exp~$echo "local: *; };" >> $lib.exp~
		$LD -G${allow_undefined_flag} -M $lib.exp -h $soname -o $lib $libobjs $deplibs $linker_flags~$rm $lib.exp'
    hardcode_libdir_flag_spec='-R$libdir'
    hardcode_shlibpath_var=no
    case $host_os in
    solaris2.[[0-5]] | solaris2.[[0-5]].*) ;;
    *) # Supported since Solaris 2.6 (maybe 2.5.1?)
      whole_archive_flag_spec='-z allextract$convenience -z defaultextract' ;;
    esac
    link_all_deplibs=yes
    ;;

  sunos4*)
    if test "x$host_vendor" = xsequent; then
      # Use $CC to link under sequent, because it throws in some extra .o
      # files that make .init and .fini sections work.
      archive_cmds='$CC -G ${wl}-h $soname -o $lib $libobjs $deplibs $compiler_flags'
    else
      archive_cmds='$LD -assert pure-text -Bstatic -o $lib $libobjs $deplibs $linker_flags'
    fi
    hardcode_libdir_flag_spec='-L$libdir'
    hardcode_direct=yes
    hardcode_minus_L=yes
    hardcode_shlibpath_var=no
    ;;

  sysv4)
    case $host_vendor in
      sni)
        archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
        hardcode_direct=yes # is this really true???
        ;;
      siemens)
        ## LD is ld it makes a PLAMLIB
        ## CC just makes a GrossModule.
        archive_cmds='$LD -G -o $lib $libobjs $deplibs $linker_flags'
        reload_cmds='$CC -r -o $output$reload_objs'
        hardcode_direct=no
        ;;
      motorola)
        archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
        hardcode_direct=no #Motorola manual says yes, but my tests say they lie
        ;;
    esac
    runpath_var='LD_RUN_PATH'
    hardcode_shlibpath_var=no
    ;;

  sysv4.3*)
    archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
    hardcode_shlibpath_var=no
    export_dynamic_flag_spec='-Bexport'
    ;;

  sysv5*)
    no_undefined_flag=' -z text'
    # $CC -shared without GNU ld will not create a library from C++
    # object files and a static libstdc++, better avoid it by now
    archive_cmds='$LD -G${allow_undefined_flag} -h $soname -o $lib $libobjs $deplibs $linker_flags'
    archive_expsym_cmds='$echo "{ global:" > $lib.exp~cat $export_symbols | sed -e "s/\(.*\)/\1;/" >> $lib.exp~$echo "local: *; };" >> $lib.exp~
		$LD -G${allow_undefined_flag} -M $lib.exp -h $soname -o $lib $libobjs $deplibs $linker_flags~$rm $lib.exp'
    hardcode_libdir_flag_spec=
    hardcode_shlibpath_var=no
    runpath_var='LD_RUN_PATH'
    ;;

  uts4*)
    archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
    hardcode_libdir_flag_spec='-L$libdir'
    hardcode_shlibpath_var=no
    ;;

  dgux*)
    archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
    hardcode_libdir_flag_spec='-L$libdir'
    hardcode_shlibpath_var=no
    ;;

  sysv4*MP*)
    if test -d /usr/nec; then
      archive_cmds='$LD -G -h $soname -o $lib $libobjs $deplibs $linker_flags'
      hardcode_shlibpath_var=no
      runpath_var=LD_RUN_PATH
      hardcode_runpath_var=yes
      ld_shlibs=yes
    fi
    ;;

  sysv4.2uw2*)
    archive_cmds='$LD -G -o $lib $libobjs $deplibs $linker_flags'
    hardcode_direct=yes
    hardcode_minus_L=no
    hardcode_shlibpath_var=no
    hardcode_runpath_var=yes
    runpath_var=LD_RUN_PATH
    ;;

  sysv5uw7* | unixware7*)
    no_undefined_flag='${wl}-z ${wl}text'
    if test "$GCC" = yes; then
      archive_cmds='$CC -shared ${wl}-h ${wl}$soname -o $lib $libobjs $deplibs $compiler_flags'
    else
      archive_cmds='$CC -G ${wl}-h ${wl}$soname -o $lib $libobjs $deplibs $compiler_flags'
    fi
    runpath_var='LD_RUN_PATH'
    hardcode_shlibpath_var=no
    ;;

  *)
    ld_shlibs=no
    ;;
  esac
fi
AC_MSG_RESULT([$ld_shlibs])
test "$ld_shlibs" = no && can_build_shared=no

# Check hardcoding attributes.
AC_MSG_CHECKING([how to hardcode library paths into programs])
hardcode_action=
if test -n "$hardcode_libdir_flag_spec" || \
   test -n "$runpath_var"; then

  # We can hardcode non-existant directories.
  if test "$hardcode_direct" != no &&
     # If the only mechanism to avoid hardcoding is shlibpath_var, we
     # have to relink, otherwise we might link with an installed library
     # when we should be linking with a yet-to-be-installed one
     ## test "$hardcode_shlibpath_var" != no &&
     test "$hardcode_minus_L" != no; then
    # Linking always hardcodes the temporary library directory.
    hardcode_action=relink
  else
    # We can link without hardcoding, and we can hardcode nonexisting dirs.
    hardcode_action=immediate
  fi
else
  # We cannot hardcode anything, or else we can only hardcode existing
  # directories.
  hardcode_action=unsupported
fi
AC_MSG_RESULT([$hardcode_action])

striplib=
old_striplib=
AC_MSG_CHECKING([whether stripping libraries is possible])
if test -n "$STRIP" && $STRIP -V 2>&1 | grep "GNU strip" >/dev/null; then
  test -z "$old_striplib" && old_striplib="$STRIP --strip-debug"
  test -z "$striplib" && striplib="$STRIP --strip-unneeded"
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

reload_cmds='$LD$reload_flag -o $output$reload_objs'
test -z "$deplibs_check_method" && deplibs_check_method=unknown

# PORTME Fill in your ld.so characteristics
AC_MSG_CHECKING([dynamic linker characteristics])
library_names_spec=
libname_spec='lib$name'
soname_spec=
postinstall_cmds=
postuninstall_cmds=
finish_cmds=
finish_eval=
shlibpath_var=
shlibpath_overrides_runpath=unknown
version_type=none
dynamic_linker="$host_os ld.so"
sys_lib_dlsearch_path_spec="/lib /usr/lib"
sys_lib_search_path_spec="/lib /usr/lib /usr/local/lib"

case $host_os in
aix3*)
  version_type=linux
  library_names_spec='${libname}${release}.so$versuffix $libname.a'
  shlibpath_var=LIBPATH

  # AIX has no versioning support, so we append a major version to the name.
  soname_spec='${libname}${release}.so$major'
  ;;

aix4* | aix5*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  hardcode_into_libs=yes
  if test "$host_cpu" = ia64; then
    # AIX 5 supports IA64
    library_names_spec='${libname}${release}.so$major ${libname}${release}.so$versuffix $libname.so'
    shlibpath_var=LD_LIBRARY_PATH
  else
    # With GCC up to 2.95.x, collect2 would create an import file
    # for dependence libraries.  The import file would start with
    # the line `#! .'.  This would cause the generated library to
    # depend on `.', always an invalid library.  This was fixed in
    # development snapshots of GCC prior to 3.0.
    case $host_os in
      aix4 | aix4.[[01]] | aix4.[[01]].*)
	if { echo '#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 97)'
	     echo ' yes '
	     echo '#endif'; } | ${CC} -E - | grep yes > /dev/null; then
	  :
	else
	  can_build_shared=no
	fi
	;;
    esac
    # AIX (on Power*) has no versioning support, so currently we can
    # not hardcode correct soname into executable. Probably we can
    # add versioning support to collect2, so additional links can
    # be useful in future.
    if test "$aix_use_runtimelinking" = yes; then
      # If using run time linking (on AIX 4.2 or later) use lib<name>.so
      # instead of lib<name>.a to let people know that these are not
      # typical AIX shared libraries.
      library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
    else
      # We preserve .a as extension for shared libraries through AIX4.2
      # and later when we are not doing run time linking.
      library_names_spec='${libname}${release}.a $libname.a'
      soname_spec='${libname}${release}.so$major'
    fi
    shlibpath_var=LIBPATH
  fi
  hardcode_into_libs=yes
  ;;

amigaos*)
  library_names_spec='$libname.ixlibrary $libname.a'
  # Create ${libname}_ixlibrary.a entries in /sys/libs.
  finish_eval='for lib in `ls $libdir/*.ixlibrary 2>/dev/null`; do libname=`$echo "X$lib" | $Xsed -e '\''s%^.*/\([[^/]]*\)\.ixlibrary$%\1%'\''`; test $rm /sys/libs/${libname}_ixlibrary.a; $show "(cd /sys/libs && $LN_S $lib ${libname}_ixlibrary.a)"; (cd /sys/libs && $LN_S $lib ${libname}_ixlibrary.a) || exit 1; done'
  ;;

beos*)
  library_names_spec='${libname}.so'
  dynamic_linker="$host_os ld.so"
  shlibpath_var=LIBRARY_PATH
  ;;

bsdi4*)
  version_type=linux
  need_version=no
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  soname_spec='${libname}${release}.so$major'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  sys_lib_search_path_spec="/shlib /usr/lib /usr/X11/lib /usr/contrib/lib /lib /usr/local/lib"
  sys_lib_dlsearch_path_spec="/shlib /usr/lib /usr/local/lib"
  export_dynamic_flag_spec=-rdynamic
  # the default ld.so.conf also contains /usr/contrib/lib and
  # /usr/X11R6/lib (/usr/X11 is a link to /usr/X11R6), but let us allow
  # libtool to hard-code these into programs
  ;;

cygwin* | mingw* | pw32*)
  version_type=windows
  need_version=no
  need_lib_prefix=no
  case $GCC,$host_os in
  yes,cygwin*)
    library_names_spec='$libname.dll.a'
    soname_spec='`echo ${libname} | sed -e 's/^lib/cyg/'``echo ${release} | sed -e 's/[[.]]/-/g'`${versuffix}.dll'
    postinstall_cmds='dlpath=`bash 2>&1 -c '\''. $dir/${file}i;echo \$dlname'\''`~
      dldir=$destdir/`dirname \$dlpath`~
      test -d \$dldir || mkdir -p \$dldir~
      $install_prog .libs/$dlname \$dldir/$dlname'
    postuninstall_cmds='dldll=`bash 2>&1 -c '\''. $file; echo \$dlname'\''`~
      dlpath=$dir/\$dldll~
       $rm \$dlpath'
    ;;
  yes,mingw*)
    library_names_spec='${libname}`echo ${release} | sed -e 's/[[.]]/-/g'`${versuffix}.dll'
    sys_lib_search_path_spec=`$CC -print-search-dirs | grep "^libraries:" | sed -e "s/^libraries://" -e "s/;/ /g" -e "s,=/,/,g"`
    ;;
  yes,pw32*)
    library_names_spec='`echo ${libname} | sed -e 's/^lib/pw/'``echo ${release} | sed -e 's/[.]/-/g'`${versuffix}.dll'
    ;;
  *)
    library_names_spec='${libname}`echo ${release} | sed -e 's/[[.]]/-/g'`${versuffix}.dll $libname.lib'
    ;;
  esac
  dynamic_linker='Win32 ld.exe'
  # FIXME: first we should search . and the directory the executable is in
  shlibpath_var=PATH
  ;;

darwin* | rhapsody*)
  dynamic_linker="$host_os dyld"
  version_type=darwin
  need_lib_prefix=no
  need_version=no
  # FIXME: Relying on posixy $() will cause problems for
  #        cross-compilation, but unfortunately the echo tests do not
  #        yet detect zsh echo's removal of \ escapes.
  library_names_spec='${libname}${release}${versuffix}.$(test .$module = .yes && echo so || echo dylib) ${libname}${release}${major}.$(test .$module = .yes && echo so || echo dylib) ${libname}.$(test .$module = .yes && echo so || echo dylib)'
  soname_spec='${libname}${release}${major}.$(test .$module = .yes && echo so || echo dylib)'
  shlibpath_overrides_runpath=yes
  shlibpath_var=DYLD_LIBRARY_PATH
  ;;

freebsd1*)
  dynamic_linker=no
  ;;

freebsd*)
  objformat=`test -x /usr/bin/objformat && /usr/bin/objformat || echo aout`
  version_type=freebsd-$objformat
  case $version_type in
    freebsd-elf*)
      library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so $libname.so'
      need_version=no
      need_lib_prefix=no
      ;;
    freebsd-*)
      library_names_spec='${libname}${release}.so$versuffix $libname.so$versuffix'
      need_version=yes
      ;;
  esac
  shlibpath_var=LD_LIBRARY_PATH
  case $host_os in
  freebsd2*)
    shlibpath_overrides_runpath=yes
    ;;
  *)
    shlibpath_overrides_runpath=no
    hardcode_into_libs=yes
    ;;
  esac
  ;;

gnu*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so${major} ${libname}.so'
  soname_spec='${libname}${release}.so$major'
  shlibpath_var=LD_LIBRARY_PATH
  hardcode_into_libs=yes
  ;;

hpux9* | hpux10* | hpux11*)
  # Give a soname corresponding to the major version so that dld.sl refuses to
  # link against other versions.
  dynamic_linker="$host_os dld.sl"
  version_type=sunos
  need_lib_prefix=no
  need_version=no
  shlibpath_var=SHLIB_PATH
  shlibpath_overrides_runpath=no # +s is required to enable SHLIB_PATH
  library_names_spec='${libname}${release}.sl$versuffix ${libname}${release}.sl$major $libname.sl'
  soname_spec='${libname}${release}.sl$major'
  # HP-UX runs *really* slowly unless shared libraries are mode 555.
  postinstall_cmds='chmod 555 $lib'
  ;;

irix5* | irix6* | nonstopux*)
  case $host_os in
    nonstopux*) version_type=nonstopux ;;
    *)          version_type=irix ;;
  esac
  need_lib_prefix=no
  need_version=no
  soname_spec='${libname}${release}.so$major'
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major ${libname}${release}.so $libname.so'
  case $host_os in
  irix5* | nonstopux*)
    libsuff= shlibsuff=
    ;;
  *)
    case $LD in # libtool.m4 will add one of these switches to LD
    *-32|*"-32 ") libsuff= shlibsuff= libmagic=32-bit;;
    *-n32|*"-n32 ") libsuff=32 shlibsuff=N32 libmagic=N32;;
    *-64|*"-64 ") libsuff=64 shlibsuff=64 libmagic=64-bit;;
    *) libsuff= shlibsuff= libmagic=never-match;;
    esac
    ;;
  esac
  shlibpath_var=LD_LIBRARY${shlibsuff}_PATH
  shlibpath_overrides_runpath=no
  sys_lib_search_path_spec="/usr/lib${libsuff} /lib${libsuff} /usr/local/lib${libsuff}"
  sys_lib_dlsearch_path_spec="/usr/lib${libsuff} /lib${libsuff}"
  ;;

# No shared lib support for Linux oldld, aout, or coff.
linux-gnuoldld* | linux-gnuaout* | linux-gnucoff*)
  dynamic_linker=no
  ;;

# This must be Linux ELF.
linux-gnu*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  soname_spec='${libname}${release}.so$major'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig -n $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=no
  # This implies no fast_install, which is unacceptable.
  # Some rework will be needed to allow for fast_install
  # before this can be enabled.
  hardcode_into_libs=yes

  # We used to test for /lib/ld.so.1 and disable shared libraries on
  # powerpc, because MkLinux only supported shared libraries with the
  # GNU dynamic linker.  Since this was broken with cross compilers,
  # most powerpc-linux boxes support dynamic linking these days and
  # people can always --disable-shared, the test was removed, and we
  # assume the GNU/Linux dynamic linker is in use.
  dynamic_linker='GNU/Linux ld.so'

  # Find out which ABI we are using (multilib Linux x86_64 hack).
  libsuff=
  case "$host_cpu" in
  x86_64*|s390x*)
    echo '[#]line __oline__ "configure"' > conftest.$ac_ext
    if AC_TRY_EVAL(ac_compile); then
      case `/usr/bin/file conftest.$ac_objext` in
      *64-bit*)
        libsuff=64
        ;;
      esac
    fi
    rm -rf conftest*
    ;;
  *)
    ;;
  esac
  sys_lib_dlsearch_path_spec="/lib${libsuff} /usr/lib${libsuff}"
  sys_lib_search_path_spec="/lib${libsuff} /usr/lib${libsuff} /usr/local/lib${libsuff}"
  ;;

netbsd*)
  version_type=sunos
  need_lib_prefix=no
  need_version=no
  if echo __ELF__ | $CC -E - | grep __ELF__ >/dev/null; then
    library_names_spec='${libname}${release}.so$versuffix ${libname}.so$versuffix'
    finish_cmds='PATH="\$PATH:/sbin" ldconfig -m $libdir'
    dynamic_linker='NetBSD (a.out) ld.so'
  else
    library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major ${libname}${release}.so ${libname}.so'
    soname_spec='${libname}${release}.so$major'
    dynamic_linker='NetBSD ld.elf_so'
  fi
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  hardcode_into_libs=yes
  ;;

newsos6)
  version_type=linux
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  ;;

openbsd*)
  version_type=sunos
  need_lib_prefix=no
  need_version=no
  if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`" || test "$host_os-$host_cpu" = "openbsd2.8-powerpc"; then
    case "$host_os" in
    openbsd2.[[89]] | openbsd2.[[89]].*)
      shlibpath_overrides_runpath=no
      ;;
    *)
      shlibpath_overrides_runpath=yes
      ;;
    esac
  else
    shlibpath_overrides_runpath=yes
  fi
  library_names_spec='${libname}${release}.so$versuffix ${libname}.so$versuffix'
  finish_cmds='PATH="\$PATH:/sbin" ldconfig -m $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

os2*)
  libname_spec='$name'
  need_lib_prefix=no
  library_names_spec='$libname.dll $libname.a'
  dynamic_linker='OS/2 ld.exe'
  shlibpath_var=LIBPATH
  ;;

osf3* | osf4* | osf5*)
  version_type=osf
  need_version=no
  soname_spec='${libname}${release}.so$major'
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  shlibpath_var=LD_LIBRARY_PATH
  sys_lib_search_path_spec="/usr/shlib /usr/ccs/lib /usr/lib/cmplrs/cc /usr/lib /usr/local/lib /var/shlib"
  sys_lib_dlsearch_path_spec="$sys_lib_search_path_spec"
  hardcode_into_libs=yes
  ;;

sco3.2v5*)
  version_type=osf
  soname_spec='${libname}${release}.so$major'
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

solaris*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  soname_spec='${libname}${release}.so$major'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  hardcode_into_libs=yes
  # ldd complains unless libraries are executable
  postinstall_cmds='chmod +x $lib'
  ;;

sunos4*)
  version_type=sunos
  library_names_spec='${libname}${release}.so$versuffix ${libname}.so$versuffix'
  finish_cmds='PATH="\$PATH:/usr/etc" ldconfig $libdir'
  shlibpath_var=LD_LIBRARY_PATH
  shlibpath_overrides_runpath=yes
  if test "$with_gnu_ld" = yes; then
    need_lib_prefix=no
  fi
  need_version=yes
  ;;

sysv4 | sysv4.2uw2* | sysv4.3* | sysv5*)
  version_type=linux
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  soname_spec='${libname}${release}.so$major'
  shlibpath_var=LD_LIBRARY_PATH
  case $host_vendor in
    sni)
      shlibpath_overrides_runpath=no
      need_lib_prefix=no
      export_dynamic_flag_spec='${wl}-Blargedynsym'
      runpath_var=LD_RUN_PATH
      ;;
    siemens)
      need_lib_prefix=no
      ;;
    motorola)
      need_lib_prefix=no
      need_version=no
      shlibpath_overrides_runpath=no
      sys_lib_search_path_spec='/lib /usr/lib /usr/ccs/lib'
      ;;
  esac
  ;;

uts4*)
  version_type=linux
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  soname_spec='${libname}${release}.so$major'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

dgux*)
  version_type=linux
  need_lib_prefix=no
  need_version=no
  library_names_spec='${libname}${release}.so$versuffix ${libname}${release}.so$major $libname.so'
  soname_spec='${libname}${release}.so$major'
  shlibpath_var=LD_LIBRARY_PATH
  ;;

sysv4*MP*)
  if test -d /usr/nec ;then
    version_type=linux
    library_names_spec='$libname.so.$versuffix $libname.so.$major $libname.so'
    soname_spec='$libname.so.$major'
    shlibpath_var=LD_LIBRARY_PATH
  fi
  ;;

*)
  dynamic_linker=no
  ;;
esac
AC_MSG_RESULT([$dynamic_linker])
test "$dynamic_linker" = no && can_build_shared=no

# Report the final consequences.
AC_MSG_CHECKING([if libtool supports shared libraries])
AC_MSG_RESULT([$can_build_shared])

AC_MSG_CHECKING([whether to build shared libraries])
test "$can_build_shared" = "no" && enable_shared=no

# On AIX, shared libraries and static libraries use the same namespace, and
# are all built from PIC.
case "$host_os" in
aix3*)
  test "$enable_shared" = yes && enable_static=no
  if test -n "$RANLIB"; then
    archive_cmds="$archive_cmds~\$RANLIB \$lib"
    postinstall_cmds='$RANLIB $lib'
  fi
  ;;

aix4*)
  if test "$host_cpu" != ia64 && test "$aix_use_runtimelinking" = no ; then
    test "$enable_shared" = yes && enable_static=no
  fi
  ;;
esac
AC_MSG_RESULT([$enable_shared])

AC_MSG_CHECKING([whether to build static libraries])
# Make sure either enable_shared or enable_static is yes.
test "$enable_shared" = yes || enable_static=yes
AC_MSG_RESULT([$enable_static])

if test "$hardcode_action" = relink; then
  # Fast installation is not supported
  enable_fast_install=no
elif test "$shlibpath_overrides_runpath" = yes ||
     test "$enable_shared" = no; then
  # Fast installation is not necessary
  enable_fast_install=needless
fi

variables_saved_for_relink="PATH $shlibpath_var $runpath_var"
if test "$GCC" = yes; then
  variables_saved_for_relink="$variables_saved_for_relink GCC_EXEC_PREFIX COMPILER_PATH LIBRARY_PATH"
fi

AC_LIBTOOL_DLOPEN_SELF

if test "$enable_shared" = yes && test "$GCC" = yes; then
  case $archive_cmds in
  *'~'*)
    # FIXME: we may have to deal with multi-command sequences.
    ;;
  '$CC '*)
    # Test whether the compiler implicitly links with -lc since on some
    # systems, -lgcc has to come before -lc. If gcc already passes -lc
    # to ld, don't add -lc before -lgcc.
    AC_MSG_CHECKING([whether -lc should be explicitly linked in])
    AC_CACHE_VAL([lt_cv_archive_cmds_need_lc],
    [$rm conftest*
    echo 'static int dummy;' > conftest.$ac_ext

    if AC_TRY_EVAL(ac_compile); then
      soname=conftest
      lib=conftest
      libobjs=conftest.$ac_objext
      deplibs=
      wl=$lt_cv_prog_cc_wl
      compiler_flags=-v
      linker_flags=-v
      verstring=
      output_objdir=.
      libname=conftest
      save_allow_undefined_flag=$allow_undefined_flag
      allow_undefined_flag=
      if AC_TRY_EVAL(archive_cmds 2\>\&1 \| grep \" -lc \" \>/dev/null 2\>\&1)
      then
	lt_cv_archive_cmds_need_lc=no
      else
	lt_cv_archive_cmds_need_lc=yes
      fi
      allow_undefined_flag=$save_allow_undefined_flag
    else
      cat conftest.err 1>&5
    fi])
    AC_MSG_RESULT([$lt_cv_archive_cmds_need_lc])
    ;;
  esac
fi
need_lc=${lt_cv_archive_cmds_need_lc-yes}

# The second clause should only fire when bootstrapping the
# libtool distribution, otherwise you forgot to ship ltmain.sh
# with your package, and you will get complaints that there are
# no rules to generate ltmain.sh.
if test -f "$ltmain"; then
  :
else
  # If there is no Makefile yet, we rely on a make rule to execute
  # `config.status --recheck' to rerun these tests and create the
  # libtool script then.
  test -f Makefile && make "$ltmain"
fi

if test -f "$ltmain"; then
  trap "$rm \"${ofile}T\"; exit 1" 1 2 15
  $rm -f "${ofile}T"

  echo creating $ofile

  # Now quote all the things that may contain metacharacters while being
  # careful not to overquote the AC_SUBSTed values.  We take copies of the
  # variables and quote the copies for generation of the libtool script.
  for var in echo old_CC old_CFLAGS SED \
    AR AR_FLAGS CC LD LN_S NM SHELL \
    reload_flag reload_cmds wl \
    pic_flag link_static_flag no_builtin_flag export_dynamic_flag_spec \
    thread_safe_flag_spec whole_archive_flag_spec libname_spec \
    library_names_spec soname_spec \
    RANLIB old_archive_cmds old_archive_from_new_cmds old_postinstall_cmds \
    old_postuninstall_cmds archive_cmds archive_expsym_cmds postinstall_cmds \
    postuninstall_cmds extract_expsyms_cmds old_archive_from_expsyms_cmds \
    old_striplib striplib file_magic_cmd export_symbols_cmds \
    deplibs_check_method allow_undefined_flag no_undefined_flag \
    finish_cmds finish_eval global_symbol_pipe global_symbol_to_cdecl \
    global_symbol_to_c_name_address \
    hardcode_libdir_flag_spec hardcode_libdir_separator  \
    sys_lib_search_path_spec sys_lib_dlsearch_path_spec \
    compiler_c_o compiler_o_lo need_locks exclude_expsyms include_expsyms; do

    case $var in
    reload_cmds | old_archive_cmds | old_archive_from_new_cmds | \
    old_postinstall_cmds | old_postuninstall_cmds | \
    export_symbols_cmds | archive_cmds | archive_expsym_cmds | \
    extract_expsyms_cmds | old_archive_from_expsyms_cmds | \
    postinstall_cmds | postuninstall_cmds | \
    finish_cmds | sys_lib_search_path_spec | sys_lib_dlsearch_path_spec)
      # Double-quote double-evaled strings.
      eval "lt_$var=\\\"\`\$echo \"X\$$var\" | \$Xsed -e \"\$double_quote_subst\" -e \"\$sed_quote_subst\" -e \"\$delay_variable_subst\"\`\\\""
      ;;
    *)
      eval "lt_$var=\\\"\`\$echo \"X\$$var\" | \$Xsed -e \"\$sed_quote_subst\"\`\\\""
      ;;
    esac
  done

  cat <<__EOF__ > "${ofile}T"
#! $SHELL

# `$echo "$ofile" | sed 's%^.*/%%'` - Provide generalized library-building support services.
# Generated automatically by $PROGRAM (GNU $PACKAGE $VERSION$TIMESTAMP)
# NOTE: Changes made to this file will be lost: look at ltmain.sh.
#
# Copyright (C) 1996-2000 Free Software Foundation, Inc.
# Originally by Gordon Matzigkeit <gord@gnu.ai.mit.edu>, 1996
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
# As a special exception to the GNU General Public License, if you
# distribute this file as part of a program that contains a
# configuration script generated by Autoconf, you may include it under
# the same distribution terms that you use for the rest of that program.

# A sed that does not truncate output.
SED=$lt_SED

# Sed that helps us avoid accidentally triggering echo(1) options like -n.
Xsed="${SED} -e s/^X//"

# The HP-UX ksh and POSIX shell print the target directory to stdout
# if CDPATH is set.
if test "X\${CDPATH+set}" = Xset; then CDPATH=:; export CDPATH; fi

# ### BEGIN LIBTOOL CONFIG

# Libtool was configured on host `(hostname || uname -n) 2>/dev/null | sed 1q`:

# Shell to use when invoking shell scripts.
SHELL=$lt_SHELL

# Whether or not to build shared libraries.
build_libtool_libs=$enable_shared

# Whether or not to build static libraries.
build_old_libs=$enable_static

# Whether or not to add -lc for building shared libraries.
build_libtool_need_lc=$need_lc

# Whether or not to optimize for fast installation.
fast_install=$enable_fast_install

# The host system.
host_alias=$host_alias
host=$host

# An echo program that does not interpret backslashes.
echo=$lt_echo

# The archiver.
AR=$lt_AR
AR_FLAGS=$lt_AR_FLAGS

# The default C compiler.
CC=$lt_CC

# Is the compiler the GNU C compiler?
with_gcc=$GCC

# The linker used to build libraries.
LD=$lt_LD

# Whether we need hard or soft links.
LN_S=$lt_LN_S

# A BSD-compatible nm program.
NM=$lt_NM

# A symbol stripping program
STRIP=$STRIP

# Used to examine libraries when file_magic_cmd begins "file"
MAGIC_CMD=$MAGIC_CMD

# Used on cygwin: DLL creation program.
DLLTOOL="$DLLTOOL"

# Used on cygwin: object dumper.
OBJDUMP="$OBJDUMP"

# Used on cygwin: assembler.
AS="$AS"

# The name of the directory that contains temporary libtool files.
objdir=$objdir

# How to create reloadable object files.
reload_flag=$lt_reload_flag
reload_cmds=$lt_reload_cmds

# How to pass a linker flag through the compiler.
wl=$lt_wl

# Object file suffix (normally "o").
objext="$ac_objext"

# Old archive suffix (normally "a").
libext="$libext"

# Executable file suffix (normally "").
exeext="$exeext"

# Additional compiler flags for building library objects.
pic_flag=$lt_pic_flag
pic_mode=$pic_mode

# Does compiler simultaneously support -c and -o options?
compiler_c_o=$lt_compiler_c_o

# Can we write directly to a .lo ?
compiler_o_lo=$lt_compiler_o_lo

# Must we lock files when doing compilation ?
need_locks=$lt_need_locks

# Do we need the lib prefix for modules?
need_lib_prefix=$need_lib_prefix

# Do we need a version for libraries?
need_version=$need_version

# Whether dlopen is supported.
dlopen_support=$enable_dlopen

# Whether dlopen of programs is supported.
dlopen_self=$enable_dlopen_self

# Whether dlopen of statically linked programs is supported.
dlopen_self_static=$enable_dlopen_self_static

# Compiler flag to prevent dynamic linking.
link_static_flag=$lt_link_static_flag

# Compiler flag to turn off builtin functions.
no_builtin_flag=$lt_no_builtin_flag

# Compiler flag to allow reflexive dlopens.
export_dynamic_flag_spec=$lt_export_dynamic_flag_spec

# Compiler flag to generate shared objects directly from archives.
whole_archive_flag_spec=$lt_whole_archive_flag_spec

# Compiler flag to generate thread-safe objects.
thread_safe_flag_spec=$lt_thread_safe_flag_spec

# Library versioning type.
version_type=$version_type

# Format of library name prefix.
libname_spec=$lt_libname_spec

# List of archive names.  First name is the real one, the rest are links.
# The last name is the one that the linker finds with -lNAME.
library_names_spec=$lt_library_names_spec

# The coded name of the library, if different from the real name.
soname_spec=$lt_soname_spec

# Commands used to build and install an old-style archive.
RANLIB=$lt_RANLIB
old_archive_cmds=$lt_old_archive_cmds
old_postinstall_cmds=$lt_old_postinstall_cmds
old_postuninstall_cmds=$lt_old_postuninstall_cmds

# Create an old-style archive from a shared archive.
old_archive_from_new_cmds=$lt_old_archive_from_new_cmds

# Create a temporary old-style archive to link instead of a shared archive.
old_archive_from_expsyms_cmds=$lt_old_archive_from_expsyms_cmds

# Commands used to build and install a shared archive.
archive_cmds=$lt_archive_cmds
archive_expsym_cmds=$lt_archive_expsym_cmds
postinstall_cmds=$lt_postinstall_cmds
postuninstall_cmds=$lt_postuninstall_cmds

# Commands to strip libraries.
old_striplib=$lt_old_striplib
striplib=$lt_striplib

# Method to check whether dependent libraries are shared objects.
deplibs_check_method=$lt_deplibs_check_method

# Command to use when deplibs_check_method == file_magic.
file_magic_cmd=$lt_file_magic_cmd

# Flag that allows shared libraries with undefined symbols to be built.
allow_undefined_flag=$lt_allow_undefined_flag

# Flag that forces no undefined symbols.
no_undefined_flag=$lt_no_undefined_flag

# Commands used to finish a libtool library installation in a directory.
finish_cmds=$lt_finish_cmds

# Same as above, but a single script fragment to be evaled but not shown.
finish_eval=$lt_finish_eval

# Take the output of nm and produce a listing of raw symbols and C names.
global_symbol_pipe=$lt_global_symbol_pipe

# Transform the output of nm in a proper C declaration
global_symbol_to_cdecl=$lt_global_symbol_to_cdecl

# Transform the output of nm in a C name address pair
global_symbol_to_c_name_address=$lt_global_symbol_to_c_name_address

# This is the shared library runtime path variable.
runpath_var=$runpath_var

# This is the shared library path variable.
shlibpath_var=$shlibpath_var

# Is shlibpath searched before the hard-coded library search path?
shlibpath_overrides_runpath=$shlibpath_overrides_runpath

# How to hardcode a shared library path into an executable.
hardcode_action=$hardcode_action

# Whether we should hardcode library paths into libraries.
hardcode_into_libs=$hardcode_into_libs

# Flag to hardcode \$libdir into a binary during linking.
# This must work even if \$libdir does not exist.
hardcode_libdir_flag_spec=$lt_hardcode_libdir_flag_spec

# Whether we need a single -rpath flag with a separated argument.
hardcode_libdir_separator=$lt_hardcode_libdir_separator

# Set to yes if using DIR/libNAME.so during linking hardcodes DIR into the
# resulting binary.
hardcode_direct=$hardcode_direct

# Set to yes if using the -LDIR flag during linking hardcodes DIR into the
# resulting binary.
hardcode_minus_L=$hardcode_minus_L

# Set to yes if using SHLIBPATH_VAR=DIR during linking hardcodes DIR into
# the resulting binary.
hardcode_shlibpath_var=$hardcode_shlibpath_var

# Variables whose values should be saved in libtool wrapper scripts and
# restored at relink time.
variables_saved_for_relink="$variables_saved_for_relink"

# Whether libtool must link a program against all its dependency libraries.
link_all_deplibs=$link_all_deplibs

# Compile-time system search path for libraries
sys_lib_search_path_spec=$lt_sys_lib_search_path_spec

# Run-time system search path for libraries
sys_lib_dlsearch_path_spec=$lt_sys_lib_dlsearch_path_spec

# Fix the shell variable \$srcfile for the compiler.
fix_srcfile_path="$fix_srcfile_path"

# Set to yes if exported symbols are required.
always_export_symbols=$always_export_symbols

# The commands to list exported symbols.
export_symbols_cmds=$lt_export_symbols_cmds

# The commands to extract the exported symbol list from a shared archive.
extract_expsyms_cmds=$lt_extract_expsyms_cmds

# Symbols that should not be listed in the preloaded symbols.
exclude_expsyms=$lt_exclude_expsyms

# Symbols that must always be exported.
include_expsyms=$lt_include_expsyms

# ### END LIBTOOL CONFIG

__EOF__

  case $host_os in
  aix3*)
    cat <<\EOF >> "${ofile}T"

# AIX sometimes has problems with the GCC collect2 program.  For some
# reason, if we set the COLLECT_NAMES environment variable, the problems
# vanish in a puff of smoke.
if test "X${COLLECT_NAMES+set}" != Xset; then
  COLLECT_NAMES=
  export COLLECT_NAMES
fi
EOF
    ;;
  esac

  case $host_os in
  cygwin* | mingw* | pw32* | os2*)
    cat <<'EOF' >> "${ofile}T"
      # This is a source program that is used to create dlls on Windows
      # Don't remove nor modify the starting and closing comments
# /* ltdll.c starts here */
# #define WIN32_LEAN_AND_MEAN
# #include <windows.h>
# #undef WIN32_LEAN_AND_MEAN
# #include <stdio.h>
#
# #ifndef __CYGWIN__
# #  ifdef __CYGWIN32__
# #    define __CYGWIN__ __CYGWIN32__
# #  endif
# #endif
#
# #ifdef __cplusplus
# extern "C" {
# #endif
# BOOL APIENTRY DllMain (HINSTANCE hInst, DWORD reason, LPVOID reserved);
# #ifdef __cplusplus
# }
# #endif
#
# #ifdef __CYGWIN__
# #include <cygwin/cygwin_dll.h>
# DECLARE_CYGWIN_DLL( DllMain );
# #endif
# HINSTANCE __hDllInstance_base;
#
# BOOL APIENTRY
# DllMain (HINSTANCE hInst, DWORD reason, LPVOID reserved)
# {
#   __hDllInstance_base = hInst;
#   return TRUE;
# }
# /* ltdll.c ends here */
	# This is a source program that is used to create import libraries
	# on Windows for dlls which lack them. Don't remove nor modify the
	# starting and closing comments
# /* impgen.c starts here */
# /*   Copyright (C) 1999-2000 Free Software Foundation, Inc.
#
#  This file is part of GNU libtool.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#  */
#
# #include <stdio.h>		/* for printf() */
# #include <unistd.h>		/* for open(), lseek(), read() */
# #include <fcntl.h>		/* for O_RDONLY, O_BINARY */
# #include <string.h>		/* for strdup() */
#
# /* O_BINARY isn't required (or even defined sometimes) under Unix */
# #ifndef O_BINARY
# #define O_BINARY 0
# #endif
#
# static unsigned int
# pe_get16 (fd, offset)
#      int fd;
#      int offset;
# {
#   unsigned char b[2];
#   lseek (fd, offset, SEEK_SET);
#   read (fd, b, 2);
#   return b[0] + (b[1]<<8);
# }
#
# static unsigned int
# pe_get32 (fd, offset)
#     int fd;
#     int offset;
# {
#   unsigned char b[4];
#   lseek (fd, offset, SEEK_SET);
#   read (fd, b, 4);
#   return b[0] + (b[1]<<8) + (b[2]<<16) + (b[3]<<24);
# }
#
# static unsigned int
# pe_as32 (ptr)
#      void *ptr;
# {
#   unsigned char *b = ptr;
#   return b[0] + (b[1]<<8) + (b[2]<<16) + (b[3]<<24);
# }
#
# int
# main (argc, argv)
#     int argc;
#     char *argv[];
# {
#     int dll;
#     unsigned long pe_header_offset, opthdr_ofs, num_entries, i;
#     unsigned long export_rva, export_size, nsections, secptr, expptr;
#     unsigned long name_rvas, nexp;
#     unsigned char *expdata, *erva;
#     char *filename, *dll_name;
#
#     filename = argv[1];
#
#     dll = open(filename, O_RDONLY|O_BINARY);
#     if (dll < 1)
# 	return 1;
#
#     dll_name = filename;
#
#     for (i=0; filename[i]; i++)
# 	if (filename[i] == '/' || filename[i] == '\\'  || filename[i] == ':')
# 	    dll_name = filename + i +1;
#
#     pe_header_offset = pe_get32 (dll, 0x3c);
#     opthdr_ofs = pe_header_offset + 4 + 20;
#     num_entries = pe_get32 (dll, opthdr_ofs + 92);
#
#     if (num_entries < 1) /* no exports */
# 	return 1;
#
#     export_rva = pe_get32 (dll, opthdr_ofs + 96);
#     export_size = pe_get32 (dll, opthdr_ofs + 100);
#     nsections = pe_get16 (dll, pe_header_offset + 4 +2);
#     secptr = (pe_header_offset + 4 + 20 +
# 	      pe_get16 (dll, pe_header_offset + 4 + 16));
#
#     expptr = 0;
#     for (i = 0; i < nsections; i++)
#     {
# 	char sname[8];
# 	unsigned long secptr1 = secptr + 40 * i;
# 	unsigned long vaddr = pe_get32 (dll, secptr1 + 12);
# 	unsigned long vsize = pe_get32 (dll, secptr1 + 16);
# 	unsigned long fptr = pe_get32 (dll, secptr1 + 20);
# 	lseek(dll, secptr1, SEEK_SET);
# 	read(dll, sname, 8);
# 	if (vaddr <= export_rva && vaddr+vsize > export_rva)
# 	{
# 	    expptr = fptr + (export_rva - vaddr);
# 	    if (export_rva + export_size > vaddr + vsize)
# 		export_size = vsize - (export_rva - vaddr);
# 	    break;
# 	}
#     }
#
#     expdata = (unsigned char*)malloc(export_size);
#     lseek (dll, expptr, SEEK_SET);
#     read (dll, expdata, export_size);
#     erva = expdata - export_rva;
#
#     nexp = pe_as32 (expdata+24);
#     name_rvas = pe_as32 (expdata+32);
#
#     printf ("EXPORTS\n");
#     for (i = 0; i<nexp; i++)
#     {
# 	unsigned long name_rva = pe_as32 (erva+name_rvas+i*4);
# 	printf ("\t%s @ %ld ;\n", erva+name_rva, 1+ i);
#     }
#
#     return 0;
# }
# /* impgen.c ends here */

EOF
    ;;
  esac

  # We use sed instead of cat because bash on DJGPP gets confused if
  # if finds mixed CR/LF and LF-only lines.  Since sed operates in
  # text mode, it properly converts lines to CR/LF.  This bash problem
  # is reportedly fixed, but why not run on old versions too?
  sed '$q' "$ltmain" >> "${ofile}T" || (rm -f "${ofile}T"; exit 1)

  mv -f "${ofile}T" "$ofile" || \
    (rm -f "$ofile" && cp "${ofile}T" "$ofile" && rm -f "${ofile}T")
  chmod +x "$ofile"
fi

])# _LT_AC_LTCONFIG_HACK

# AC_LIBTOOL_DLOPEN - enable checks for dlopen support
AC_DEFUN([AC_LIBTOOL_DLOPEN], [AC_BEFORE([$0],[AC_LIBTOOL_SETUP])])

# AC_LIBTOOL_WIN32_DLL - declare package support for building win32 dll's
AC_DEFUN([AC_LIBTOOL_WIN32_DLL], [AC_BEFORE([$0], [AC_LIBTOOL_SETUP])])

# AC_ENABLE_SHARED - implement the --enable-shared flag
# Usage: AC_ENABLE_SHARED[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN([AC_ENABLE_SHARED],
[define([AC_ENABLE_SHARED_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(shared,
changequote(<<, >>)dnl
<<  --enable-shared[=PKGS]  build shared libraries [default=>>AC_ENABLE_SHARED_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case $enableval in
yes) enable_shared=yes ;;
no) enable_shared=no ;;
*)
  enable_shared=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_shared=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_shared=AC_ENABLE_SHARED_DEFAULT)dnl
])

# AC_DISABLE_SHARED - set the default shared flag to --disable-shared
AC_DEFUN([AC_DISABLE_SHARED],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_SHARED(no)])

# AC_ENABLE_STATIC - implement the --enable-static flag
# Usage: AC_ENABLE_STATIC[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN([AC_ENABLE_STATIC],
[define([AC_ENABLE_STATIC_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(static,
changequote(<<, >>)dnl
<<  --enable-static[=PKGS]  build static libraries [default=>>AC_ENABLE_STATIC_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case $enableval in
yes) enable_static=yes ;;
no) enable_static=no ;;
*)
  enable_static=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_static=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_static=AC_ENABLE_STATIC_DEFAULT)dnl
])

# AC_DISABLE_STATIC - set the default static flag to --disable-static
AC_DEFUN([AC_DISABLE_STATIC],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_STATIC(no)])


# AC_ENABLE_FAST_INSTALL - implement the --enable-fast-install flag
# Usage: AC_ENABLE_FAST_INSTALL[(DEFAULT)]
#   Where DEFAULT is either `yes' or `no'.  If omitted, it defaults to
#   `yes'.
AC_DEFUN([AC_ENABLE_FAST_INSTALL],
[define([AC_ENABLE_FAST_INSTALL_DEFAULT], ifelse($1, no, no, yes))dnl
AC_ARG_ENABLE(fast-install,
changequote(<<, >>)dnl
<<  --enable-fast-install[=PKGS]  optimize for fast installation [default=>>AC_ENABLE_FAST_INSTALL_DEFAULT],
changequote([, ])dnl
[p=${PACKAGE-default}
case $enableval in
yes) enable_fast_install=yes ;;
no) enable_fast_install=no ;;
*)
  enable_fast_install=no
  # Look at the argument we got.  We use all the common list separators.
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS="${IFS}:,"
  for pkg in $enableval; do
    if test "X$pkg" = "X$p"; then
      enable_fast_install=yes
    fi
  done
  IFS="$ac_save_ifs"
  ;;
esac],
enable_fast_install=AC_ENABLE_FAST_INSTALL_DEFAULT)dnl
])

# AC_DISABLE_FAST_INSTALL - set the default to --disable-fast-install
AC_DEFUN([AC_DISABLE_FAST_INSTALL],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
AC_ENABLE_FAST_INSTALL(no)])

# AC_LIBTOOL_PICMODE - implement the --with-pic flag
# Usage: AC_LIBTOOL_PICMODE[(MODE)]
#   Where MODE is either `yes' or `no'.  If omitted, it defaults to
#   `both'.
AC_DEFUN([AC_LIBTOOL_PICMODE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
pic_mode=ifelse($#,1,$1,default)])


# AC_PATH_TOOL_PREFIX - find a file program which can recognise shared library
AC_DEFUN([AC_PATH_TOOL_PREFIX],
[AC_MSG_CHECKING([for $1])
AC_CACHE_VAL(lt_cv_path_MAGIC_CMD,
[case $MAGIC_CMD in
  /*)
  lt_cv_path_MAGIC_CMD="$MAGIC_CMD" # Let the user override the test with a path.
  ;;
  ?:/*)
  lt_cv_path_MAGIC_CMD="$MAGIC_CMD" # Let the user override the test with a dos path.
  ;;
  *)
  ac_save_MAGIC_CMD="$MAGIC_CMD"
  IFS="${IFS=   }"; ac_save_ifs="$IFS"; IFS=":"
dnl $ac_dummy forces splitting on constant user-supplied paths.
dnl POSIX.2 word splitting is done only on the output of word expansions,
dnl not every word.  This closes a longstanding sh security hole.
  ac_dummy="ifelse([$2], , $PATH, [$2])"
  for ac_dir in $ac_dummy; do
    test -z "$ac_dir" && ac_dir=.
    if test -f $ac_dir/$1; then
      lt_cv_path_MAGIC_CMD="$ac_dir/$1"
      if test -n "$file_magic_test_file"; then
	case $deplibs_check_method in
	"file_magic "*)
	  file_magic_regex="`expr \"$deplibs_check_method\" : \"file_magic \(.*\)\"`"
	  MAGIC_CMD="$lt_cv_path_MAGIC_CMD"
	  if eval $file_magic_cmd \$file_magic_test_file 2> /dev/null |
	    egrep "$file_magic_regex" > /dev/null; then
	    :
	  else
	    cat <<EOF 1>&2

*** Warning: the command libtool uses to detect shared libraries,
*** $file_magic_cmd, produces output that libtool cannot recognize.
*** The result is that libtool may fail to recognize shared libraries
*** as such.  This will affect the creation of libtool libraries that
*** depend on shared libraries, but programs linked with such libtool
*** libraries will work regardless of this problem.  Nevertheless, you
*** may want to report the problem to your system manager and/or to
*** bug-libtool@gnu.org

EOF
	  fi ;;
	esac
      fi
      break
    fi
  done
  IFS="$ac_save_ifs"
  MAGIC_CMD="$ac_save_MAGIC_CMD"
  ;;
esac])
MAGIC_CMD="$lt_cv_path_MAGIC_CMD"
if test -n "$MAGIC_CMD"; then
  AC_MSG_RESULT($MAGIC_CMD)
else
  AC_MSG_RESULT(no)
fi
])


# AC_PATH_MAGIC - find a file program which can recognise a shared library
AC_DEFUN([AC_PATH_MAGIC],
[AC_REQUIRE([AC_CHECK_TOOL_PREFIX])dnl
AC_PATH_TOOL_PREFIX(${ac_tool_prefix}file, /usr/bin:$PATH)
if test -z "$lt_cv_path_MAGIC_CMD"; then
  if test -n "$ac_tool_prefix"; then
    AC_PATH_TOOL_PREFIX(file, /usr/bin:$PATH)
  else
    MAGIC_CMD=:
  fi
fi
])


# AC_PROG_LD - find the path to the GNU or non-GNU linker
AC_DEFUN([AC_PROG_LD],
[AC_ARG_WITH(gnu-ld,
[  --with-gnu-ld           assume the C compiler uses GNU ld [default=no]],
test "$withval" = no || with_gnu_ld=yes, with_gnu_ld=no)
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_CANONICAL_HOST])dnl
AC_REQUIRE([AC_CANONICAL_BUILD])dnl
AC_REQUIRE([_LT_AC_LIBTOOL_SYS_PATH_SEPARATOR])dnl
ac_prog=ld
if test "$GCC" = yes; then
  # Check if gcc -print-prog-name=ld gives a path.
  AC_MSG_CHECKING([for ld used by GCC])
  case $host in
  *-*-mingw*)
    # gcc leaves a trailing carriage return which upsets mingw
    ac_prog=`($CC -print-prog-name=ld) 2>&5 | tr -d '\015'` ;;
  *)
    ac_prog=`($CC -print-prog-name=ld) 2>&5` ;;
  esac
  case $ac_prog in
    # Accept absolute paths.
    [[\\/]]* | [[A-Za-z]]:[[\\/]]*)
      re_direlt='/[[^/]][[^/]]*/\.\./'
      # Canonicalize the path of ld
      ac_prog=`echo $ac_prog| sed 's%\\\\%/%g'`
      while echo $ac_prog | grep "$re_direlt" > /dev/null 2>&1; do
	ac_prog=`echo $ac_prog| sed "s%$re_direlt%/%"`
      done
      test -z "$LD" && LD="$ac_prog"
      ;;
  "")
    # If it fails, then pretend we aren't using GCC.
    ac_prog=ld
    ;;
  *)
    # If it is relative, then search for the first ld in PATH.
    with_gnu_ld=unknown
    ;;
  esac
elif test "$with_gnu_ld" = yes; then
  AC_MSG_CHECKING([for GNU ld])
else
  AC_MSG_CHECKING([for non-GNU ld])
fi
AC_CACHE_VAL(lt_cv_path_LD,
[if test -z "$LD"; then
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=$PATH_SEPARATOR
  for ac_dir in $PATH; do
    test -z "$ac_dir" && ac_dir=.
    if test -f "$ac_dir/$ac_prog" || test -f "$ac_dir/$ac_prog$ac_exeext"; then
      lt_cv_path_LD="$ac_dir/$ac_prog"
      # Check to see if the program is GNU ld.  I'd rather use --version,
      # but apparently some GNU ld's only accept -v.
      # Break only if it was the GNU/non-GNU ld that we prefer.
      if "$lt_cv_path_LD" -v 2>&1 < /dev/null | egrep '(GNU|with BFD)' > /dev/null; then
	test "$with_gnu_ld" != no && break
      else
	test "$with_gnu_ld" != yes && break
      fi
    fi
  done
  IFS="$ac_save_ifs"
else
  lt_cv_path_LD="$LD" # Let the user override the test with a path.
fi])
LD="$lt_cv_path_LD"
if test -n "$LD"; then
  AC_MSG_RESULT($LD)
else
  AC_MSG_RESULT(no)
fi
test -z "$LD" && AC_MSG_ERROR([no acceptable ld found in \$PATH])
AC_PROG_LD_GNU
])

# AC_PROG_LD_GNU -
AC_DEFUN([AC_PROG_LD_GNU],
[AC_CACHE_CHECK([if the linker ($LD) is GNU ld], lt_cv_prog_gnu_ld,
[# I'd rather use --version here, but apparently some GNU ld's only accept -v.
if $LD -v 2>&1 </dev/null | egrep '(GNU|with BFD)' 1>&5; then
  lt_cv_prog_gnu_ld=yes
else
  lt_cv_prog_gnu_ld=no
fi])
with_gnu_ld=$lt_cv_prog_gnu_ld
])

# AC_PROG_LD_RELOAD_FLAG - find reload flag for linker
#   -- PORTME Some linkers may need a different reload flag.
AC_DEFUN([AC_PROG_LD_RELOAD_FLAG],
[AC_CACHE_CHECK([for $LD option to reload object files], lt_cv_ld_reload_flag,
[lt_cv_ld_reload_flag='-r'])
reload_flag=$lt_cv_ld_reload_flag
test -n "$reload_flag" && reload_flag=" $reload_flag"
])

# AC_DEPLIBS_CHECK_METHOD - how to check for library dependencies
#  -- PORTME fill in with the dynamic library characteristics
AC_DEFUN([AC_DEPLIBS_CHECK_METHOD],
[AC_CACHE_CHECK([how to recognise dependent libraries],
lt_cv_deplibs_check_method,
[lt_cv_file_magic_cmd='$MAGIC_CMD'
lt_cv_file_magic_test_file=
lt_cv_deplibs_check_method='unknown'
# Need to set the preceding variable on all platforms that support
# interlibrary dependencies.
# 'none' -- dependencies not supported.
# `unknown' -- same as none, but documents that we really don't know.
# 'pass_all' -- all dependencies passed with no checks.
# 'test_compile' -- check by making test program.
# 'file_magic [[regex]]' -- check by looking for files in library path
# which responds to the $file_magic_cmd with a given egrep regex.
# If you have `file' or equivalent on your system and you're not sure
# whether `pass_all' will *always* work, you probably want this one.

case $host_os in
aix4* | aix5*)
  lt_cv_deplibs_check_method=pass_all
  ;;

beos*)
  lt_cv_deplibs_check_method=pass_all
  ;;

bsdi4*)
  lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (shared object|dynamic lib)'
  lt_cv_file_magic_cmd='/usr/bin/file -L'
  lt_cv_file_magic_test_file=/shlib/libc.so
  ;;

cygwin* | mingw* | pw32*)
  lt_cv_deplibs_check_method='file_magic file format pei*-i386(.*architecture: i386)?'
  lt_cv_file_magic_cmd='$OBJDUMP -f'
  ;;

darwin* | rhapsody*)
  lt_cv_deplibs_check_method='file_magic Mach-O dynamically linked shared library'
  lt_cv_file_magic_cmd='/usr/bin/file -L'
  case "$host_os" in
  rhapsody* | darwin1.[[012]])
    lt_cv_file_magic_test_file=`echo /System/Library/Frameworks/System.framework/Versions/*/System | head -1`
    ;;
  *) # Darwin 1.3 on
    lt_cv_file_magic_test_file='/usr/lib/libSystem.dylib'
    ;;
  esac
  ;;

freebsd*)
  if echo __ELF__ | $CC -E - | grep __ELF__ > /dev/null; then
    case $host_cpu in
    i*86 )
      # Not sure whether the presence of OpenBSD here was a mistake.
      # Let's accept both of them until this is cleared up.
      lt_cv_deplibs_check_method='file_magic (FreeBSD|OpenBSD)/i[[3-9]]86 (compact )?demand paged shared library'
      lt_cv_file_magic_cmd=/usr/bin/file
      lt_cv_file_magic_test_file=`echo /usr/lib/libc.so.*`
      ;;
    esac
  else
    lt_cv_deplibs_check_method=pass_all
  fi
  ;;

gnu*)
  lt_cv_deplibs_check_method=pass_all
  ;;

hpux10.20*|hpux11*)
  lt_cv_deplibs_check_method='file_magic (s[[0-9]][[0-9]][[0-9]]|PA-RISC[[0-9]].[[0-9]]) shared library'
  lt_cv_file_magic_cmd=/usr/bin/file
  lt_cv_file_magic_test_file=/usr/lib/libc.sl
  ;;

irix5* | irix6* | nonstopux*)
  case $host_os in
  irix5* | nonstopux*)
    # this will be overridden with pass_all, but let us keep it just in case
    lt_cv_deplibs_check_method="file_magic ELF 32-bit MSB dynamic lib MIPS - version 1"
    ;;
  *)
    case $LD in
    *-32|*"-32 ") libmagic=32-bit;;
    *-n32|*"-n32 ") libmagic=N32;;
    *-64|*"-64 ") libmagic=64-bit;;
    *) libmagic=never-match;;
    esac
    # this will be overridden with pass_all, but let us keep it just in case
    lt_cv_deplibs_check_method="file_magic ELF ${libmagic} MSB mips-[[1234]] dynamic lib MIPS - version 1"
    ;;
  esac
  lt_cv_file_magic_test_file=`echo /lib${libsuff}/libc.so*`
  lt_cv_deplibs_check_method=pass_all
  ;;

# This must be Linux ELF.
linux-gnu*)
  case $host_cpu in
  alpha* | hppa* | i*86 | mips | mipsel | powerpc* | sparc* | ia64* | s390* | x86_64*)
    lt_cv_deplibs_check_method=pass_all ;;
  *)
    # glibc up to 2.1.1 does not perform some relocations on ARM
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB (shared object|dynamic lib )' ;;
  esac
  lt_cv_file_magic_test_file=`echo /lib/libc.so* /lib/libc-*.so`
  ;;

netbsd*)
  if echo __ELF__ | $CC -E - | grep __ELF__ > /dev/null; then
    lt_cv_deplibs_check_method='match_pattern /lib[[^/\.]]+\.so\.[[0-9]]+\.[[0-9]]+$'
  else
    lt_cv_deplibs_check_method='match_pattern /lib[[^/\.]]+\.so$'
  fi
  ;;

newos6*)
  lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (executable|dynamic lib)'
  lt_cv_file_magic_cmd=/usr/bin/file
  lt_cv_file_magic_test_file=/usr/lib/libnls.so
  ;;

openbsd*)
  lt_cv_file_magic_cmd=/usr/bin/file
  lt_cv_file_magic_test_file=`echo /usr/lib/libc.so.*`
  if test -z "`echo __ELF__ | $CC -E - | grep __ELF__`" || test "$host_os-$host_cpu" = "openbsd2.8-powerpc"; then
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB shared object'
  else
    lt_cv_deplibs_check_method='file_magic OpenBSD.* shared library'
  fi
  ;;

osf3* | osf4* | osf5*)
  # this will be overridden with pass_all, but let us keep it just in case
  lt_cv_deplibs_check_method='file_magic COFF format alpha shared library'
  lt_cv_file_magic_test_file=/shlib/libc.so
  lt_cv_deplibs_check_method=pass_all
  ;;

sco3.2v5*)
  lt_cv_deplibs_check_method=pass_all
  ;;

solaris*)
  lt_cv_deplibs_check_method=pass_all
  lt_cv_file_magic_test_file=/lib/libc.so
  ;;

sysv5uw[[78]]* | sysv4*uw2*)
  lt_cv_deplibs_check_method=pass_all
  ;;

sysv4 | sysv4.2uw2* | sysv4.3* | sysv5*)
  case $host_vendor in
  motorola)
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[ML]]SB (shared object|dynamic lib) M[[0-9]][[0-9]]* Version [[0-9]]'
    lt_cv_file_magic_test_file=`echo /usr/lib/libc.so*`
    ;;
  ncr)
    lt_cv_deplibs_check_method=pass_all
    ;;
  sequent)
    lt_cv_file_magic_cmd='/bin/file'
    lt_cv_deplibs_check_method='file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB (shared object|dynamic lib )'
    ;;
  sni)
    lt_cv_file_magic_cmd='/bin/file'
    lt_cv_deplibs_check_method="file_magic ELF [[0-9]][[0-9]]*-bit [[LM]]SB dynamic lib"
    lt_cv_file_magic_test_file=/lib/libc.so
    ;;
  siemens)
    lt_cv_deplibs_check_method=pass_all
    ;;
  esac
  ;;
esac
])
file_magic_cmd=$lt_cv_file_magic_cmd
deplibs_check_method=$lt_cv_deplibs_check_method
])


# AC_PROG_NM - find the path to a BSD-compatible name lister
AC_DEFUN([AC_PROG_NM],
[AC_REQUIRE([_LT_AC_LIBTOOL_SYS_PATH_SEPARATOR])dnl
AC_MSG_CHECKING([for BSD-compatible nm])
AC_CACHE_VAL(lt_cv_path_NM,
[if test -n "$NM"; then
  # Let the user override the test.
  lt_cv_path_NM="$NM"
else
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=$PATH_SEPARATOR
  for ac_dir in $PATH /usr/ccs/bin /usr/ucb /bin; do
    test -z "$ac_dir" && ac_dir=.
    tmp_nm=$ac_dir/${ac_tool_prefix}nm
    if test -f $tmp_nm || test -f $tmp_nm$ac_exeext ; then
      # Check to see if the nm accepts a BSD-compat flag.
      # Adding the `sed 1q' prevents false positives on HP-UX, which says:
      #   nm: unknown option "B" ignored
      # Tru64's nm complains that /dev/null is an invalid object file
      if ($tmp_nm -B /dev/null 2>&1 | sed '1q'; exit 0) | egrep '(/dev/null|Invalid file or object type)' >/dev/null; then
	lt_cv_path_NM="$tmp_nm -B"
	break
      elif ($tmp_nm -p /dev/null 2>&1 | sed '1q'; exit 0) | egrep /dev/null >/dev/null; then
	lt_cv_path_NM="$tmp_nm -p"
	break
      else
	lt_cv_path_NM=${lt_cv_path_NM="$tmp_nm"} # keep the first match, but
	continue # so that we can try to find one that supports BSD flags
      fi
    fi
  done
  IFS="$ac_save_ifs"
  test -z "$lt_cv_path_NM" && lt_cv_path_NM=nm
fi])
NM="$lt_cv_path_NM"
AC_MSG_RESULT([$NM])
])

# AC_CHECK_LIBM - check for math library
AC_DEFUN([AC_CHECK_LIBM],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl
LIBM=
case $host in
*-*-beos* | *-*-cygwin* | *-*-pw32*)
  # These system don't have libm
  ;;
*-ncr-sysv4.3*)
  AC_CHECK_LIB(mw, _mwvalidcheckl, LIBM="-lmw")
  AC_CHECK_LIB(m, main, LIBM="$LIBM -lm")
  ;;
*)
  AC_CHECK_LIB(m, main, LIBM="-lm")
  ;;
esac
])

# AC_LIBLTDL_CONVENIENCE[(dir)] - sets LIBLTDL to the link flags for
# the libltdl convenience library and LTDLINCL to the include flags for
# the libltdl header and adds --enable-ltdl-convenience to the
# configure arguments.  Note that LIBLTDL and LTDLINCL are not
# AC_SUBSTed, nor is AC_CONFIG_SUBDIRS called.  If DIR is not
# provided, it is assumed to be `libltdl'.  LIBLTDL will be prefixed
# with '${top_builddir}/' and LTDLINCL will be prefixed with
# '${top_srcdir}/' (note the single quotes!).  If your package is not
# flat and you're not using automake, define top_builddir and
# top_srcdir appropriately in the Makefiles.
AC_DEFUN([AC_LIBLTDL_CONVENIENCE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  case $enable_ltdl_convenience in
  no) AC_MSG_ERROR([this package needs a convenience libltdl]) ;;
  "") enable_ltdl_convenience=yes
      ac_configure_args="$ac_configure_args --enable-ltdl-convenience" ;;
  esac
  LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdlc.la
  LTDLINCL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
  # For backwards non-gettext consistent compatibility...
  INCLTDL="$LTDLINCL"
])

# AC_LIBLTDL_INSTALLABLE[(dir)] - sets LIBLTDL to the link flags for
# the libltdl installable library and LTDLINCL to the include flags for
# the libltdl header and adds --enable-ltdl-install to the configure
# arguments.  Note that LIBLTDL and LTDLINCL are not AC_SUBSTed, nor is
# AC_CONFIG_SUBDIRS called.  If DIR is not provided and an installed
# libltdl is not found, it is assumed to be `libltdl'.  LIBLTDL will
# be prefixed with '${top_builddir}/' and LTDLINCL will be prefixed
# with '${top_srcdir}/' (note the single quotes!).  If your package is
# not flat and you're not using automake, define top_builddir and
# top_srcdir appropriately in the Makefiles.
# In the future, this macro may have to be called after AC_PROG_LIBTOOL.
AC_DEFUN([AC_LIBLTDL_INSTALLABLE],
[AC_BEFORE([$0],[AC_LIBTOOL_SETUP])dnl
  AC_CHECK_LIB(ltdl, main,
  [test x"$enable_ltdl_install" != xyes && enable_ltdl_install=no],
  [if test x"$enable_ltdl_install" = xno; then
     AC_MSG_WARN([libltdl not installed, but installation disabled])
   else
     enable_ltdl_install=yes
   fi
  ])
  if test x"$enable_ltdl_install" = x"yes"; then
    ac_configure_args="$ac_configure_args --enable-ltdl-install"
    LIBLTDL='${top_builddir}/'ifelse($#,1,[$1],['libltdl'])/libltdl.la
    LTDLINCL='-I${top_srcdir}/'ifelse($#,1,[$1],['libltdl'])
  else
    ac_configure_args="$ac_configure_args --enable-ltdl-install=no"
    LIBLTDL="-lltdl"
    LTDLINCL=
  fi
  # For backwards non-gettext consistent compatibility...
  INCLTDL="$LTDLINCL"
])

# old names
AC_DEFUN([AM_PROG_LIBTOOL],   [AC_PROG_LIBTOOL])
AC_DEFUN([AM_ENABLE_SHARED],  [AC_ENABLE_SHARED($@)])
AC_DEFUN([AM_ENABLE_STATIC],  [AC_ENABLE_STATIC($@)])
AC_DEFUN([AM_DISABLE_SHARED], [AC_DISABLE_SHARED($@)])
AC_DEFUN([AM_DISABLE_STATIC], [AC_DISABLE_STATIC($@)])
AC_DEFUN([AM_PROG_LD],        [AC_PROG_LD])
AC_DEFUN([AM_PROG_NM],        [AC_PROG_NM])

# This is just to silence aclocal about the macro not being used
ifelse([AC_DISABLE_FAST_INSTALL])

# NOTE: This macro has been submitted for inclusion into   #
#  GNU Autoconf as AC_PROG_SED.  When it is available in   #
#  a released version of Autoconf we should remove this    #
#  macro and use it instead.                               #
# LT_AC_PROG_SED
# --------------
# Check for a fully-functional sed program, that truncates
# as few characters as possible.  Prefer GNU sed if found.
AC_DEFUN([LT_AC_PROG_SED],
[AC_MSG_CHECKING([for a sed that does not truncate output])
AC_CACHE_VAL(lt_cv_path_SED,
[# Loop through the user's path and test for sed and gsed.
# Then use that list of sed's as ones to test for truncation.
as_executable_p="test -f"
as_save_IFS=$IFS; IFS=$PATH_SEPARATOR
for as_dir in $PATH
do
  IFS=$as_save_IFS
  test -z "$as_dir" && as_dir=.
  for ac_prog in sed gsed; do
    for ac_exec_ext in '' $ac_executable_extensions; do
      if $as_executable_p "$as_dir/$ac_prog$ac_exec_ext"; then
        _sed_list="$_sed_list $as_dir/$ac_prog$ac_exec_ext"
      fi
    done
  done
done

  # Create a temporary directory, and hook for its removal unless debugging.
$debug ||
{
  trap 'exit_status=$?; rm -rf $tmp && exit $exit_status' 0
  trap '{ (exit 1); exit 1; }' 1 2 13 15
}

# Create a (secure) tmp directory for tmp files.
: ${TMPDIR=/tmp}
{
  tmp=`(umask 077 && mktemp -d -q "$TMPDIR/sedXXXXXX") 2>/dev/null` &&
  test -n "$tmp" && test -d "$tmp"
}  ||
{
  tmp=$TMPDIR/sed$$-$RANDOM
  (umask 077 && mkdir $tmp)
} ||
{
   echo "$me: cannot create a temporary directory in $TMPDIR" >&2
   { (exit 1); exit 1; }
}
  _max=0
  _count=0
  # Add /usr/xpg4/bin/sed as it is typically found on Solaris
  # along with /bin/sed that truncates output.
  for _sed in $_sed_list /usr/xpg4/bin/sed; do
    test ! -f ${_sed} && break
    cat /dev/null > "$tmp/sed.in"
    _count=0
    echo ${ECHO_N-$ac_n} "0123456789${ECHO_C-$ac_c}" >"$tmp/sed.in"
    # Check for GNU sed and select it if it is found.
    if "${_sed}" --version 2>&1 < /dev/null | egrep '(GNU)' > /dev/null; then
      lt_cv_path_SED=${_sed}
      break
    fi
    while true; do
      cat "$tmp/sed.in" "$tmp/sed.in" >"$tmp/sed.tmp"
      mv "$tmp/sed.tmp" "$tmp/sed.in"
      cp "$tmp/sed.in" "$tmp/sed.nl"
      echo >>"$tmp/sed.nl"
      ${_sed} -e 's/a$//' < "$tmp/sed.nl" >"$tmp/sed.out" || break
      cmp -s "$tmp/sed.out" "$tmp/sed.nl" || break
      # 40000 chars as input seems more than enough
      test $_count -gt 10 && break
      _count=`expr $_count + 1`
      if test $_count -gt $_max; then
        _max=$_count
        lt_cv_path_SED=$_sed
      fi
    done
  done
  rm -rf "$tmp"
])
if test "X$SED" != "X"; then
  lt_cv_path_SED=$SED
else
  SED=$lt_cv_path_SED
fi
AC_MSG_RESULT([$SED])
])

