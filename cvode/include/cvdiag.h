/*******************************************************************
 *                                                                 *
 * File          : cvdiag.h                                        *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 07 February 2004                                *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvode/LICENSE                         *
 *-----------------------------------------------------------------*
 * This is the header file for the CVODE diagonal linear solver,   *
 * CVDIAG.                                                         *
 *                                                                 *
 *******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvdiag_h
#define _cvdiag_h

#include <stdio.h>
#include "cvode.h"
#include "sundialstypes.h"
#include "nvector.h"

 
/******************************************************************
 *                                                                *
 * Function : CVDiag                                              *
 *----------------------------------------------------------------*
 * A call to the CVDiag function links the main CVODE integrator  *
 * with the CVDIAG linear solver.                                 *
 *                                                                *
 * cvode_mem is the pointer to CVODE memory returned by           *
 *              CVodeCreate.                                      *
 *                                                                *
 * The return values of CVDiag are:                               *
 *    SUCCESS   = 0  if successful                                *
 *    LMEM_FAIL = -1 if there was a memory allocation failure     *
 *                                                                *
 ******************************************************************/

int CVDiag(void *cvode_mem);

/******************************************************************
 * Optional outputs from the CVDIAG linear solver                 *
 *----------------------------------------------------------------*
 *                                                                *
 * CVDiagGetIntWorkSpace returns the integer workspace used by    *
 *     CVDIAG.                                                    *
 * CVDiagGetRealWorkSpace returns the real workspace used by      *
 *     CVDIAG.                                                    *
 * CVDiagGetNumRhsEvals returns the number of calls to the user   *
 *     f routine due to finite difference Jacobian evaluation.    *
 * Note: the number of diagonal approximate Jacobians formed is   *
 * equal to the number of CVDiagSetup calls.                      *
 * This number is available through CVodeGetNumLinSolvSetups.     *
 *                                                                *
 ******************************************************************/

int CVDiagGetIntWorkSpace(void *cvode_mem, long int *leniwDI);
int CVDiagGetRealWorkSpace(void *cvode_mem, long int *lenrwDI);
int CVDiagGetNumRhsEvals(void *cvode_mem, int *nfevalsDI);


/******************************************************************
 *                                                                *           
 * Types : CVDiagMemRec, CVDiagMem                                *
 *----------------------------------------------------------------*
 * The type CVDiagMem is pointer to a CVDiagMemRec. This          *
 * structure contains CVDiag solver-specific data.                *
 *                                                                *
 ******************************************************************/

typedef struct {

  realtype di_gammasv; /* gammasv = gamma at the last call to setup */
                       /* or solve                                  */

  N_Vector di_M;       /* M = (I - gamma J)^{-1} , gamma = h / l1   */

  N_Vector di_bit;     /* temporary storage vector                  */

  N_Vector di_bitcomp; /* temporary storage vector                  */

  int di_nfeDI;        /* no. of calls to f due to difference 
                          quotient diagonal Jacobian approximation  */

} CVDiagMemRec, *CVDiagMem;

 
#endif

#ifdef __cplusplus
}
#endif
