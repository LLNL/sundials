/*******************************************************************
 *                                                                 *
 * File          : cvbandpre.h                                     *
 * Programmers   : Michael Wittman, Alan C. Hindmarsh, and         *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvode/LICENSE                         *
 *-----------------------------------------------------------------*
 * This is the header file for the CVBANDPRE module, which         *
 * provides a banded difference quotient Jacobian-based            *
 * preconditioner and solver routines for use with CVSPGMR.        *
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
 *   Details of the calls to CVodeMalloc, CVSpgmr, and CVode are   *
 *   available in the CVODE User Guide.                            *
 *   To use these routines, the sequence of calls in the user      *
 *   main program should be as follows:                            *
 *                                                                 *
 *   #include "cvbandpre.h"                                        *
 *   #include "nvector_serial.h"                                   *
 *   ...                                                           *
 *   M_Env machEnv;                                                *
 *   CVBandPreData bp_data;                                        *
 *   ...                                                           *
 *   machEnv = M_EnvInit_Serial(...);                              *
 *   ...                                                           *
 *   cvode_mem = CVodeMalloc(...);                                 *
 *   ...                                                           *
 *   bp_data = CVBandPreAlloc(N, f, f_data, mu, ml, cvode_mem);    *
 *   ...                                                           *
 *   flag = CVSpgmr(cvode_mem, pretype, gstype, maxl, delt,        *
 *           CVBandPrecond, CVBandPSolve, bp_data);                *
 *   ...                                                           *
 *   flag = CVode(...);                                            *
 *   ...                                                           *
 *   CVBandPreFree(bp_data);                                       *
 *   ...                                                           *
 *   M_EnvFree_Serial(machEnv);                                    *
 *   ...                                                           *
 *   CVodeFree(cvode_mem);                                         *
 *                                                                 *
 * Notes:                                                          *
 * (1) Include this file for the CVBandPreData type definition.    *
 * (2) In the CVBandPreAlloc call, the arguments N, f, and f_data  * 
 *     are the same as in the call to CVodeMalloc.                 *
 * (3) In the CVSpgmr call, the user is free to specify the inputs *
 *     pretype and gstype, and the optional inputs maxl and delt.  *
 *     But the last three arguments must be as shown, with the     *
 *     last argument being the pointer returned by CVBandPreAlloc. *
 * (4) The CVBandPrecond and CVBandPSolve functions are never      *
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


/************* CVBandPreData type definition ************/

typedef struct {
  /* Data set by user in CVBandPreAlloc: */
  RhsFn f;
  void *f_data;
  integertype ml, mu;

  /* Data set by CVBandPrecond: */
  BandMat savedJ;
  BandMat savedP;
  integertype *pivots;
} *CVBandPreData;


/******************************************************************
 *                                                                *
 * Function : CVBandPreAlloc                                      *
 *----------------------------------------------------------------*
 * CVBandPreAlloc allocates and initializes a CVBandPreData       *
 * structure to be passed to CVSpgmr (and subsequently used by    *
 * CVBandPrecond and CVBandPSolve).                               *
 *                                                                *
 * The parameters of CVBandPreAlloc are as follows:               *
 *                                                                *
 * N       is the length of all vector arguments.                 *
 *                                                                *
 * f       is the right hand side function.                       *
 *                                                                *
 * f_data  is a pointer to the optional extra data for f.         *
 *                                                                *
 * mu      is the upper half bandwidth.                           *
 *                                                                *
 * ml      is the lower half bandwidth.                           *
 *                                                                *
 * CVBandPreAlloc returns the storage pointer (type CVBandPreData)*
 * or NULL if the request for storage cannot be satisfied.        *
 *                                                                *
 * NOTE: The band preconditioner assumes a serial implementation  *
 *       of the NVECTOR package. Therefore, CVBandPreAlloc will   *
 *       first test for a compatible N_Vector internal            *
 *       representation by checking (1) the machine environment   *
 *       ID tag and (2) that the functions N_VMake, N_VDispose,   *
 *       N_VGetData, and N_VSetData are implemented.              *
 *                                                                *
 ******************************************************************/

CVBandPreData CVBandPreAlloc(integertype N, RhsFn f, void *f_data,
                             integertype mu, integertype ml, 
                             void *cvode_mem);


/******************************************************************
 *                                                                *
 * Function : CVReInitBandPre                                     *
 *----------------------------------------------------------------*
 * CVReInitBandPre re-initializes the CVBANDPRE module when       *
 * solving a  sequence of problems of the same size with          *
 * CVSPGMR/CVBANDPRE, provided there is no change in N, mu, or ml.*
 * After solving one problem, and after calling CVReInit to       *
 * re-initialize CVODE for a subsequent problem, call             *
 * CVReInitBandPre.  Then call CVReInitSpgmr or CVSpgmr if        *
 * necessary, depending on changes made in the CVSpgmr            *
 * parameters, before calling CVode.                              *
 *                                                                *
 * The first argument to CVReInitBandPre must be the pointer      *
 * bpdata that was returned by CVBandPreAlloc.  All other         *
 * arguments have the same names and meanings as in               *
 * CVBandPreAlloc.                                                *
 *                                                                *
 * The return value of CVReInitBandPre is 0, indicating success.  *
 ******************************************************************/

int CVReInitBandPre(CVBandPreData bpdata, integertype N, RhsFn f,
                    void *f_data, integertype mu, integertype ml);


/******************************************************************
 *                                                                *
 * Function : CVBandPreFree                                       *
 *----------------------------------------------------------------*
 * CVBandPreFree frees the memory allocated by CVBandPreAlloc in  *
 * the argument pdata.                                            *
 *                                                                *
 ******************************************************************/

void CVBandPreFree(CVBandPreData pdata);



/* Prototypes of CVBandPrecond and CVBandPSolve */

  
int CVBandPrecond(integertype N, realtype t, N_Vector y, N_Vector fy, 
                  booleantype jok, booleantype *jcurPtr, realtype gamma, 
                  N_Vector ewt, realtype h, realtype uround, 
                  long int *nfePtr, void *bp_data,
                  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);


int CVBandPSolve(integertype N, realtype t, N_Vector y, N_Vector fy,
               N_Vector vtemp, realtype gamma, N_Vector ewt, realtype delta,
               long int *nfePtr, N_Vector r, int lr, void *bp_data, N_Vector z);


#endif

#ifdef __cplusplus
}
#endif
