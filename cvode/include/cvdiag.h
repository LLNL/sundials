/*******************************************************************
 *                                                                 *
 * File          : cvdiag.h                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvode/LICENSE                         *
 *-----------------------------------------------------------------*
 * This is the header file for the CVODE diagonal linear solver,   *
 * CVDIAG.                                                         *
 *                                                                 *
 * Note: The type integer must be large enough to store the value  *
 * of the linear system size N.                                    *
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
 * CVDIAG solver statistics indices                               *
 *----------------------------------------------------------------*
 * The following enumeration gives a symbolic name to each        *
 * CVDIAG statistic. The symbolic names are used as indices into  *
 * the iopt and ropt arrays passed to CVodeMalloc.                *
 * The CVDIAG statistics are:                                     *
 *                                                                *
 * iopt[DIAG_LRW] : size (in realtype words) of real workspace    *
 *                  vectors used by this solver.                  *
 *                                                                *
 * iopt[DIAG_LIW] : size (in integertype words) of integer        *
 *                  workspace vectors used by this solver.        *
 *                                                                *
 * The number of diagonal approximate Jacobians formed is equal   *
 * to the number of CVDiagSetup calls. This number is available   *
 * in cv_iopt[NSETUPS].                                           *
 *                                                                *
 ******************************************************************/
 
enum { DIAG_LRW=CVODE_IOPT_SIZE, DIAG_LIW };

 
/******************************************************************
 *                                                                *
 * Function : CVDiag                                              *
 *----------------------------------------------------------------*
 * A call to the CVDiag function links the main CVODE integrator  *
 * with the CVDIAG linear solver.                                 *
 *                                                                *
 * cvode_mem is the pointer to CVODE memory returned by           *
 *              CVodeMalloc.                                      *
 *                                                                *
 * The return values of CVDiag are:                               *
 *    SUCCESS   = 0  if successful                                *
 *    LMEM_FAIL = -1 if there was a memory allocation failure     *
 *                                                                *
 ******************************************************************/

int CVDiag(void *cvode_mem);
 
#endif

#ifdef __cplusplus
}
#endif
