/*******************************************************************
 *                                                                 *
 * File          : cvbandpre.h                                     *
 * Programmers   : Michael Wittman, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 31 July 2003                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvode/LICENSE                         *
 *-----------------------------------------------------------------*
 * This is the header file for the CVBANDPRE module, which         *
 * provides a banded difference quotient Jacobian-based            *
 * preconditioner and solver routines for use with CVSSPGMR.       *
 *                                                                 *
 * Summary:                                                        *
 * These routines provide a band matrix preconditioner based on    *
 * difference quotients of the ODE right-hand side function f.     *
 * The user supplies parameters                                    *
 *   mu = upper half-bandwidth (number of super-diagonals)         *
 *   ml = lower half-bandwidth (number of sub-diagonals)           *
 * The routines generate a band matrix of bandwidth ml + mu + 1    *
 * and use this to form a preconditioner for use with the Krylov   *
 * linear solver in CVSPGMR.  Although this matrix is intended     *
 * to approximate the Jacobian df/dy, it may be a very crude       *
 * approximation.  The true Jacobian need not be banded, or its    *
 * true bandwith may be larger than ml + mu + 1, as long as the    *
 * banded approximation generated here is sufficiently accurate    *
 * to speed convergence as a preconditioner.                       *
 *                                                                 *
 * Usage:                                                          *
 *   The following is a summary of the usage of this module.       *
 *   Details of the calls to CVodeCreate, CVodeMalloc, CVSpgmr,    *
 *   and CVode are available in the CVODE User Guide.              *
 *   To use these routines, the sequence of calls in the user      *
 *   main program should be as follows:                            *
 *                                                                 *
 *   #include "cvbandpre.h"                                        *
 *   #include "nvector_serial.h"                                   *
 *   ...                                                           *
 *   NV_Spec nvspec;                                               *
 *   void *bp_data;                                                *
 *   ...                                                           *
 *   nvspec = NV_SpecInit_Serial(...);                             *
 *   ...                                                           *
 *   cvode_mem = CVodeCreate(...);                                 *
 *   ier = CVodeMalloc(...);                                       *
 *   ...                                                           *
 *   bp_data = CVBandPrecAlloc(N, f, f_data, mu, ml, cvode_mem);   *
 *   ...                                                           *
 *   flag = CVSpgmr(cvode_mem, pretype, maxl);                     *
 *   flag = CVSpgmrSetPrecSetupFn(cvode_mem, CVBandPrecSetup);     *
 *   flag = CVSpgmrSetPrecSolveFn(cvode_mem, CVBandPrecSolve);     *
 *   flag = CVSpgmrSetPrecData(cvode_mem, bp_data);                *
 *   ...                                                           *
 *   flag = CVode(...);                                            *
 *   ...                                                           *
 *   CVBandPrecFree(bp_data);                                      *
 *   ...                                                           *
 *   NV_SpecFree_Serial(nvspec);                                   *
 *   ...                                                           *
 *   CVodeFree(cvode_mem);                                         *
 *                                                                 *
 * Notes:                                                          *
 * (1) Include this file for the CVBandPrecData type definition.   *
 * (2) In the CVBandPrecAlloc call, the arguments N, f, and f_data *
 *     are the same as in the call to CVodeMalloc.                 *
 * (3) In the CVSpgmr call, the user is free to specify the inputs *
 *     pretype and gstype, and the optional inputs maxl and delt.  *
 *     But the last three arguments must be as shown, with the     *
 *     last argument being the pointer returned by CVBandPreAlloc. *
 * (4) The CVBandPrecSetup and CVBandPrecSolve functions are never *
 *     called by the user explicitly; they are simply passed to    *
 *     the CVSpgmr function.                                       *
 *                                                                 *
 *******************************************************************/

 
#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvbandpre_h
#define _cvbandpre_h

#include "cvode.h"
#include "sundialstypes.h"
#include "nvector.h"
#include "band.h"


/************* CVBandPrecData type definition ************/

typedef struct {

  /* Data set by user in CVBandPrecAlloc: */
  integertype N;
  integertype ml, mu;

  /* Data set by CVBandPrecSetup: */
  BandMat savedJ;
  BandMat savedP;
  integertype *pivots;

  /* Rhs calls */
  int nfeBP;

  /* Pointer to cvode_mem */
  CVodeMem cv_mem;

} *CVBandPrecData;


/******************************************************************
 *                                                                *
 * Function : CVBandPrecAlloc                                     *
 *----------------------------------------------------------------*
 * CVBandPrecAlloc allocates and initializes a CVBandPrecData     *
 * structure to be passed to CVSpgmr (and subsequently used by    *
 * CVBandPrecSetup and CVBandPrecSolve).                          *
 *                                                                *
 * The parameters of CVBandPrecAlloc are as follows:              *
 *                                                                *
 * n       is the length of all vector arguments.                 *
 *                                                                *
 * mu      is the upper half bandwidth.                           *
 *                                                                *
 * ml      is the lower half bandwidth.                           *
 *                                                                *
 * CVBandPrecAlloc returns the storage pointer of type            *
 * CVBandPrecData or NULL if the request for storage cannot be    *
 * satisfied.                                                     *
 *                                                                *
 * NOTE: The band preconditioner assumes a serial implementation  *
 *       of the NVECTOR package. Therefore, CVBandPrecAlloc will  *
 *       first test for a compatible N_Vector internal            *
 *       representation by checking (1) the vector specification  *
 *       ID tag and (2) that the functions N_VGetData, and        *
 *       N_VSetData are implemented.                              *
 *                                                                *
 ******************************************************************/

void *CVBandPrecAlloc(void *cvode_mem, integertype n,
                      integertype mu, integertype ml);


/******************************************************************
 * Function : CVBandPrecFree                                      *
 *----------------------------------------------------------------*
 * CVBandPrecFree frees the memory allocated by CVBandPrecAlloc   *
 * in the argument pdata.                                         *
 *                                                                *
 ******************************************************************/

void CVBandPrecFree(void *bp_data);

/******************************************************************
 * Function : CVBandPrecGet*                                      *
 *----------------------------------------------------------------*
 *                                                                *
 ******************************************************************/

int CVBandPrecGetIntWorkSpace(void *bp_data, long int *leniwBP);
int CVBandPrecGetRealWorkSpace(void *bp_data, long int *lenrwBP);
int CVBandPrecGetNumRhsEvals(void *bp_data, int *nfevalsBP);

/* Return values for CVBandPrecGet* functions */
/* OKAY = 0 */
enum { BP_NO_PDATA = -1 };

/* Prototypes of CVBandPrecSetup and CVBandPrecSolve */
  
int CVBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                    booleantype jok, booleantype *jcurPtr, 
                    realtype gamma, void *bp_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


int CVBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                    N_Vector r, N_Vector z, 
                    realtype gamma, realtype delta,
                    int lr, void *bp_data, N_Vector tmp);


#endif

#ifdef __cplusplus
}
#endif
