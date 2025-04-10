/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.0
 *
 * This file is not intended to be easily readable and contains a number of
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG
 * interface file instead.
 * ----------------------------------------------------------------------------- */

/* ---------------------------------------------------------------
 * Programmer(s): Auto-generated by swig.
 * ---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -------------------------------------------------------------*/

/* -----------------------------------------------------------------------------
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 * ----------------------------------------------------------------------------- */

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
# if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#  define SWIGTEMPLATEDISAMBIGUATOR template
# elif defined(__HP_aCC)
/* Needed even with `aCC -AA' when `aCC -V' reports HP ANSI C++ B3910B A.03.55 */
/* If we find a maximum version that requires this, the test would be __HP_aCC <= 35500 for A.03.55 */
#  define SWIGTEMPLATEDISAMBIGUATOR template
# else
#  define SWIGTEMPLATEDISAMBIGUATOR
# endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__))
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__))
# else
#   define SWIGUNUSED
# endif
#endif

#ifndef SWIG_MSC_UNSUPPRESS_4505
# if defined(_MSC_VER)
#   pragma warning(disable : 4505) /* unreferenced local function has been removed */
# endif
#endif

#ifndef SWIGUNUSEDPARM
# ifdef __cplusplus
#   define SWIGUNUSEDPARM(p)
# else
#   define SWIGUNUSEDPARM(p) p SWIGUNUSED
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* qualifier for exported *const* global data variables*/
#ifndef SWIGEXTERN
# ifdef __cplusplus
#   define SWIGEXTERN extern
# else
#   define SWIGEXTERN
# endif
#endif

/* exporting methods */
#if defined(__GNUC__)
#  if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#    ifndef GCC_HASCLASSVISIBILITY
#      define GCC_HASCLASSVISIBILITY
#    endif
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif
#endif

/* Deal with Microsoft's attempt at deprecating C standard runtime functions */
#if !defined(SWIG_NO_CRT_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
# define _CRT_SECURE_NO_DEPRECATE
#endif

/* Deal with Microsoft's attempt at deprecating methods in the standard C++ library */
#if !defined(SWIG_NO_SCL_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_SCL_SECURE_NO_DEPRECATE)
# define _SCL_SECURE_NO_DEPRECATE
#endif

/* Deal with Apple's deprecated 'AssertMacros.h' from Carbon-framework */
#if defined(__APPLE__) && !defined(__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES)
# define __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES 0
#endif

/* Intel's compiler complains if a variable which was never initialised is
 * cast to void, which is a common idiom which we use to indicate that we
 * are aware a variable isn't used.  So we just silence that warning.
 * See: https://github.com/swig/swig/issues/192 for more discussion.
 */
#ifdef __INTEL_COMPILER
# pragma warning disable 592
#endif

/*  Errors in SWIG */
#define  SWIG_UnknownError    	   -1
#define  SWIG_IOError        	   -2
#define  SWIG_RuntimeError   	   -3
#define  SWIG_IndexError     	   -4
#define  SWIG_TypeError      	   -5
#define  SWIG_DivisionByZero 	   -6
#define  SWIG_OverflowError  	   -7
#define  SWIG_SyntaxError    	   -8
#define  SWIG_ValueError     	   -9
#define  SWIG_SystemError    	   -10
#define  SWIG_AttributeError 	   -11
#define  SWIG_MemoryError    	   -12
#define  SWIG_NullReferenceError   -13




#include <assert.h>
#define SWIG_exception_impl(DECL, CODE, MSG, RETURNNULL) \
 { printf("In " DECL ": " MSG); assert(0); RETURNNULL; }


enum {
    SWIG_MEM_OWN = 0x01,
    SWIG_MEM_RVALUE = 0x02,
    SWIG_MEM_CONST = 0x04
};


#define SWIG_check_mutable(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL) \
    if ((SWIG_CLASS_WRAPPER).cmemflags & SWIG_MEM_CONST) { \
        SWIG_exception_impl(FUNCNAME, SWIG_TypeError, \
            "Cannot pass const " TYPENAME " (class " FNAME ") " \
            "as a mutable reference", \
            RETURNNULL); \
    }


#define SWIG_check_nonnull(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL) \
  if (!(SWIG_CLASS_WRAPPER).cptr) { \
    SWIG_exception_impl(FUNCNAME, SWIG_TypeError, \
                        "Cannot pass null " TYPENAME " (class " FNAME ") " \
                        "as a reference", RETURNNULL); \
  }


#define SWIG_check_mutable_nonnull(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL) \
    SWIG_check_nonnull(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL); \
    SWIG_check_mutable(SWIG_CLASS_WRAPPER, TYPENAME, FNAME, FUNCNAME, RETURNNULL);


#include <stdio.h>
#if defined(_MSC_VER) || defined(__BORLANDC__) || defined(_WATCOM)
# ifndef snprintf
#  define snprintf _snprintf
# endif
#endif


/* Support for the `contract` feature.
 *
 * Note that RETURNNULL is first because it's inserted via a 'Replaceall' in
 * the fortran.cxx file.
 */
#define SWIG_contract_assert(RETURNNULL, EXPR, MSG) \
 if (!(EXPR)) { SWIG_exception_impl("$decl", SWIG_ValueError, MSG, RETURNNULL); } 


#define SWIGVERSION 0x040000 
#define SWIG_VERSION SWIGVERSION


#define SWIG_as_voidptr(a) (void *)((const void *)(a)) 
#define SWIG_as_voidptrptr(a) ((void)SWIG_as_voidptr(*a),(void**)(a)) 


#include "sundials/sundials_adaptcontroller.h"


#include "sunadaptcontroller/sunadaptcontroller_mrihtol.h"


typedef struct {
    void* cptr;
    int cmemflags;
} SwigClassWrapper;


SWIGINTERN SwigClassWrapper SwigClassWrapper_uninitialized() {
    SwigClassWrapper result;
    result.cptr = NULL;
    result.cmemflags = 0;
    return result;
}


#include <stdlib.h>
#ifdef _MSC_VER
# ifndef strtoull
#  define strtoull _strtoui64
# endif
# ifndef strtoll
#  define strtoll _strtoi64
# endif
#endif


#include <string.h>


SWIGINTERN void SWIG_assign(SwigClassWrapper* self, SwigClassWrapper other) {
  if (self->cptr == NULL) {
    /* LHS is unassigned */
    if (other.cmemflags & SWIG_MEM_RVALUE) {
      /* Capture pointer from RHS, clear 'moving' flag */
      self->cptr = other.cptr;
      self->cmemflags = other.cmemflags & (~SWIG_MEM_RVALUE);
    } else {
      /* Become a reference to the other object */
      self->cptr = other.cptr;
      self->cmemflags = other.cmemflags & (~SWIG_MEM_OWN);
    }
  } else if (other.cptr == NULL) {
    /* Replace LHS with a null pointer */
    free(self->cptr);
    *self = SwigClassWrapper_uninitialized();
  } else {
    if (self->cmemflags & SWIG_MEM_OWN) {
      free(self->cptr);
    }
    self->cptr = other.cptr;
    if (other.cmemflags & SWIG_MEM_RVALUE) {
      /* Capture RHS */
      self->cmemflags = other.cmemflags & ~SWIG_MEM_RVALUE;
    } else {
      /* Point to RHS */
      self->cmemflags = other.cmemflags & ~SWIG_MEM_OWN;
    }
  }
}

SWIGEXPORT void _wrap_SUNAdaptControllerContent_MRIHTol__HControl_set(SwigClassWrapper const *farg1, SUNAdaptController farg2) {
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  SUNAdaptController arg2 = (SUNAdaptController) 0 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::HControl", return );
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  arg2 = (SUNAdaptController)(farg2);
  if (arg1) (arg1)->HControl = arg2;
}


SWIGEXPORT SUNAdaptController _wrap_SUNAdaptControllerContent_MRIHTol__HControl_get(SwigClassWrapper const *farg1) {
  SUNAdaptController fresult ;
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  SUNAdaptController result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::HControl", return 0);
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  result = (SUNAdaptController) ((arg1)->HControl);
  fresult = result;
  return fresult;
}


SWIGEXPORT void _wrap_SUNAdaptControllerContent_MRIHTol__TolControl_set(SwigClassWrapper const *farg1, SUNAdaptController farg2) {
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  SUNAdaptController arg2 = (SUNAdaptController) 0 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::TolControl", return );
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  arg2 = (SUNAdaptController)(farg2);
  if (arg1) (arg1)->TolControl = arg2;
}


SWIGEXPORT SUNAdaptController _wrap_SUNAdaptControllerContent_MRIHTol__TolControl_get(SwigClassWrapper const *farg1) {
  SUNAdaptController fresult ;
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  SUNAdaptController result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::TolControl", return 0);
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  result = (SUNAdaptController) ((arg1)->TolControl);
  fresult = result;
  return fresult;
}


SWIGEXPORT void _wrap_SUNAdaptControllerContent_MRIHTol__inner_max_relch_set(SwigClassWrapper const *farg1, double const *farg2) {
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  sunrealtype arg2 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::inner_max_relch", return );
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  arg2 = (sunrealtype)(*farg2);
  if (arg1) (arg1)->inner_max_relch = arg2;
}


SWIGEXPORT double _wrap_SUNAdaptControllerContent_MRIHTol__inner_max_relch_get(SwigClassWrapper const *farg1) {
  double fresult ;
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  sunrealtype result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::inner_max_relch", return 0);
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  result = (sunrealtype) ((arg1)->inner_max_relch);
  fresult = (sunrealtype)(result);
  return fresult;
}


SWIGEXPORT void _wrap_SUNAdaptControllerContent_MRIHTol__inner_min_tolfac_set(SwigClassWrapper const *farg1, double const *farg2) {
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  sunrealtype arg2 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::inner_min_tolfac", return );
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  arg2 = (sunrealtype)(*farg2);
  if (arg1) (arg1)->inner_min_tolfac = arg2;
}


SWIGEXPORT double _wrap_SUNAdaptControllerContent_MRIHTol__inner_min_tolfac_get(SwigClassWrapper const *farg1) {
  double fresult ;
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  sunrealtype result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::inner_min_tolfac", return 0);
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  result = (sunrealtype) ((arg1)->inner_min_tolfac);
  fresult = (sunrealtype)(result);
  return fresult;
}


SWIGEXPORT void _wrap_SUNAdaptControllerContent_MRIHTol__inner_max_tolfac_set(SwigClassWrapper const *farg1, double const *farg2) {
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  sunrealtype arg2 ;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::inner_max_tolfac", return );
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  arg2 = (sunrealtype)(*farg2);
  if (arg1) (arg1)->inner_max_tolfac = arg2;
}


SWIGEXPORT double _wrap_SUNAdaptControllerContent_MRIHTol__inner_max_tolfac_get(SwigClassWrapper const *farg1) {
  double fresult ;
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  sunrealtype result;
  
  SWIG_check_mutable_nonnull(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::inner_max_tolfac", return 0);
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  result = (sunrealtype) ((arg1)->inner_max_tolfac);
  fresult = (sunrealtype)(result);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_new_SUNAdaptControllerContent_MRIHTol_() {
  SwigClassWrapper fresult ;
  struct SUNAdaptControllerContent_MRIHTol_ *result = 0 ;
  
  result = (struct SUNAdaptControllerContent_MRIHTol_ *)calloc(1, sizeof(struct SUNAdaptControllerContent_MRIHTol_));
  fresult.cptr = result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (1 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT void _wrap_delete_SUNAdaptControllerContent_MRIHTol_(SwigClassWrapper *farg1) {
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  
  SWIG_check_mutable(*farg1, "struct SUNAdaptControllerContent_MRIHTol_ *", "SUNAdaptControllerContent_MRIHTol_", "SUNAdaptControllerContent_MRIHTol_::~SUNAdaptControllerContent_MRIHTol_()", return );
  arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *)(farg1->cptr);
  free((char *) arg1);
}


SWIGEXPORT void _wrap_SUNAdaptControllerContent_MRIHTol__op_assign__(SwigClassWrapper *farg1, SwigClassWrapper const *farg2) {
  struct SUNAdaptControllerContent_MRIHTol_ *arg1 = (struct SUNAdaptControllerContent_MRIHTol_ *) 0 ;
  struct SUNAdaptControllerContent_MRIHTol_ *arg2 = 0 ;
  
  (void)sizeof(arg1);
  (void)sizeof(arg2);
  SWIG_assign(farg1, *farg2);
  
}


SWIGEXPORT SUNAdaptController _wrap_FSUNAdaptController_MRIHTol(SUNAdaptController farg1, SUNAdaptController farg2, void *farg3) {
  SUNAdaptController fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  SUNAdaptController arg2 = (SUNAdaptController) 0 ;
  SUNContext arg3 = (SUNContext) 0 ;
  SUNAdaptController result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (SUNAdaptController)(farg2);
  arg3 = (SUNContext)(farg3);
  result = (SUNAdaptController)SUNAdaptController_MRIHTol(arg1,arg2,arg3);
  fresult = result;
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_SetParams_MRIHTol(SUNAdaptController farg1, double const *farg2, double const *farg3, double const *farg4) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  sunrealtype arg2 ;
  sunrealtype arg3 ;
  sunrealtype arg4 ;
  SUNErrCode result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (sunrealtype)(*farg2);
  arg3 = (sunrealtype)(*farg3);
  arg4 = (sunrealtype)(*farg4);
  result = (SUNErrCode)SUNAdaptController_SetParams_MRIHTol(arg1,arg2,arg3,arg4);
  fresult = (SUNErrCode)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_GetSlowController_MRIHTol(SUNAdaptController farg1, void *farg2) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  SUNAdaptController *arg2 = (SUNAdaptController *) 0 ;
  SUNErrCode result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (SUNAdaptController *)(farg2);
  result = (SUNErrCode)SUNAdaptController_GetSlowController_MRIHTol(arg1,arg2);
  fresult = (SUNErrCode)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_GetFastController_MRIHTol(SUNAdaptController farg1, void *farg2) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  SUNAdaptController *arg2 = (SUNAdaptController *) 0 ;
  SUNErrCode result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (SUNAdaptController *)(farg2);
  result = (SUNErrCode)SUNAdaptController_GetFastController_MRIHTol(arg1,arg2);
  fresult = (SUNErrCode)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_GetType_MRIHTol(SUNAdaptController farg1) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  SUNAdaptController_Type result;
  
  arg1 = (SUNAdaptController)(farg1);
  result = (SUNAdaptController_Type)SUNAdaptController_GetType_MRIHTol(arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_EstimateStepTol_MRIHTol(SUNAdaptController farg1, double const *farg2, double const *farg3, int const *farg4, double const *farg5, double const *farg6, double *farg7, double *farg8) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  sunrealtype arg2 ;
  sunrealtype arg3 ;
  int arg4 ;
  sunrealtype arg5 ;
  sunrealtype arg6 ;
  sunrealtype *arg7 = (sunrealtype *) 0 ;
  sunrealtype *arg8 = (sunrealtype *) 0 ;
  int result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (sunrealtype)(*farg2);
  arg3 = (sunrealtype)(*farg3);
  arg4 = (int)(*farg4);
  arg5 = (sunrealtype)(*farg5);
  arg6 = (sunrealtype)(*farg6);
  arg7 = (sunrealtype *)(farg7);
  arg8 = (sunrealtype *)(farg8);
  result = (int)SUNAdaptController_EstimateStepTol_MRIHTol(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_Reset_MRIHTol(SUNAdaptController farg1) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  int result;
  
  arg1 = (SUNAdaptController)(farg1);
  result = (int)SUNAdaptController_Reset_MRIHTol(arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_SetDefaults_MRIHTol(SUNAdaptController farg1) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  int result;
  
  arg1 = (SUNAdaptController)(farg1);
  result = (int)SUNAdaptController_SetDefaults_MRIHTol(arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_Write_MRIHTol(SUNAdaptController farg1, void *farg2) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  FILE *arg2 = (FILE *) 0 ;
  int result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (FILE *)(farg2);
  result = (int)SUNAdaptController_Write_MRIHTol(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_SetErrorBias_MRIHTol(SUNAdaptController farg1, double const *farg2) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  sunrealtype arg2 ;
  int result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (sunrealtype)(*farg2);
  result = (int)SUNAdaptController_SetErrorBias_MRIHTol(arg1,arg2);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_UpdateMRIHTol_MRIHTol(SUNAdaptController farg1, double const *farg2, double const *farg3, double const *farg4, double const *farg5) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  sunrealtype arg2 ;
  sunrealtype arg3 ;
  sunrealtype arg4 ;
  sunrealtype arg5 ;
  int result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (sunrealtype)(*farg2);
  arg3 = (sunrealtype)(*farg3);
  arg4 = (sunrealtype)(*farg4);
  arg5 = (sunrealtype)(*farg5);
  result = (int)SUNAdaptController_UpdateMRIHTol_MRIHTol(arg1,arg2,arg3,arg4,arg5);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_FSUNAdaptController_Space_MRIHTol(SUNAdaptController farg1, long *farg2, long *farg3) {
  int fresult ;
  SUNAdaptController arg1 = (SUNAdaptController) 0 ;
  long *arg2 = (long *) 0 ;
  long *arg3 = (long *) 0 ;
  int result;
  
  arg1 = (SUNAdaptController)(farg1);
  arg2 = (long *)(farg2);
  arg3 = (long *)(farg3);
  result = (int)SUNAdaptController_Space_MRIHTol(arg1,arg2,arg3);
  fresult = (int)(result);
  return fresult;
}



