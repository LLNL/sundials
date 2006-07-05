/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:32:33 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the Fortran interface include file for the rootfinding
 * feature of CVODE.
 * -----------------------------------------------------------------
 */

/*
 * ==============================================================================
 *
 *                   FCVROOT Interface Package
 *
 * The FCVROOT interface package allows programs written in FORTRAN to
 * use the rootfinding feature of the CVODE solver module.
 *
 * The user-callable functions constituting the FCVROOT package are the
 * following: FCVROOTINIT, FCVROOTINFO, and FCVROOTFREE. The corresponding
 * CVODE subroutine called by each interface function is given below.
 *
 *   -----------------      -----------------------
 *  | FCVROOT routine |    | CVODE function called |
 *   -----------------      -----------------------
 *      FCVROOTINIT     ->     CVodeRootInit
 *      FCVROOTINFO     ->     CVodeGetRootInfo
 *      FCVROOTFREE     ->     CVodeRootInit
 *
 * FCVROOTFN is a user-supplied subroutine defining the functions whose
 * roots are sought.
 *
 * ==============================================================================
 *
 *                     Usage of the FCVROOT Interface Package
 *
 * 1. In order to use the rootfinding feature of the CVODE package the user must
 * define the following subroutine:
 *
 *   SUBROUTINE FCVROOTFN (T, Y, G, IPAR, RPAR, IER)
 *   DIMENSION Y(*), G(*), IPAR(*), RPAR(*)
 *
 * The arguments are:
 *   T = independent variable value t  [input]
 *   Y = dependent variable vector y  [input]
 *   G = function values g(t,y)  [output]
 *   IPAR, RPAR = user (integer and real) data [input/output]
 *   IER = return flag (0 for success, a non-zero value if an error occurred.)
 *
 * 2. After calling FCVMALLOC but prior to calling FCVODE, the user must
 * allocate and initialize memory for the FCVROOT module by making the
 * following call:
 *
 *   CALL FCVROOTINIT (NRTFN, IER)
 *
 * The arguments are:
 *   NRTFN = total number of root functions  [input]
 *   IER   = return completion flag (0 = success, -1 = CVODE memory NULL and
 *           -11  memory allocation error)  [output]
 *
 * 3. After calling FCVODE, to see whether a root was found, test the FCVODE
 * return flag IER.  The value IER = 2 means one or more roots were found.
 *
 * 4. If a root was found, and if NRTFN > 1, then to determine which root
 * functions G(*) were found to have a root, make the following call:
 *     CALL FCVROOTINFO (NRTFN, INFO, IER)
 * The arguments are:
 *   NRTFN = total number of root functions  [input]
 *   INFO  = integer array of length NRTFN, with values 0 or 1 [output]
 *           For i = 1,...,NRTFN, G(i) was found to have a root if INFO(i) = 1.
 *   IER   = completion flag (0 = success,  negative = failure)
 *
 * 5. The total number of calls made to the root function (FCVROOTFN), NGE,
 * can be obtained from IOUT(12).
 *
 * If the FCVODE/CVODE memory block is reinitialized to solve a different
 * problem via a call to FCVREINIT, then the counter variable NGE is cleared
 * (reset to zero).
 *
 * 6. To free the memory resources allocated by a prior call to FCVROOTINIT make
 * the following call:
 *   CALL FCVROOTFREE
 * See the CVODE documentation for additional information.
 *
 * ==============================================================================
 */

#ifndef _FCVROOT_H
#define _FCVROOT_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* header files */

#include <sundials/sundials_nvector.h> /* definition of type N_Vector          */
#include <sundials/sundials_types.h>   /* definition of SUNDIALS type realtype */

/* Definitions of interface function names */

#if defined(F77_FUNC)

#define FCV_ROOTINIT F77_FUNC(fcvrootinit, FCVROOTINIT)
#define FCV_ROOTINFO F77_FUNC(fcvrootinfo, FCVROOTINFO)
#define FCV_ROOTFREE F77_FUNC(fcvrootfree, FCVROOTFREE)
#define FCV_ROOTFN   F77_FUNC(fcvrootfn, FCVROOTFN)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FCV_ROOTINIT fcvrootinit
#define FCV_ROOTINFO fcvrootinfo
#define FCV_ROOTFREE fcvrootfree
#define FCV_ROOTFN   fcvrootfn

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FCV_ROOTINIT FCVROOTINIT
#define FCV_ROOTINFO FCVROOTINFO
#define FCV_ROOTFREE FCVROOTFREE
#define FCV_ROOTFN   FCVROOTFN

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FCV_ROOTINIT fcvrootinit_
#define FCV_ROOTINFO fcvrootinfo_
#define FCV_ROOTFREE fcvrootfree_
#define FCV_ROOTFN   fcvrootfn_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FCV_ROOTINIT FCVROOTINIT_
#define FCV_ROOTINFO FCVROOTINFO_
#define FCV_ROOTFREE FCVROOTFREE_
#define FCV_ROOTFN   FCVROOTFN_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FCV_ROOTINIT fcvrootinit__
#define FCV_ROOTINFO fcvrootinfo__
#define FCV_ROOTFREE fcvrootfree__
#define FCV_ROOTFN   fcvrootfn__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FCV_ROOTINIT FCVROOTINIT__
#define FCV_ROOTINFO FCVROOTINFO__
#define FCV_ROOTFREE FCVROOTFREE__
#define FCV_ROOTFN   FCVROOTFN__

#endif

/* Prototypes of exported function */

void FCV_ROOTINIT(int *nrtfn, int *ier);
void FCV_ROOTINFO(int *nrtfn, int *info, int *ier);
void FCV_ROOTFREE(void);

/* Prototype of function called by CVODE module */

int FCVrootfunc(realtype t, N_Vector y, realtype *gout, void *g_data);

#ifdef __cplusplus
}
#endif


#endif
