/*******************************************************************
 *                                                                 *
 * File          : cvsdense.h                                      *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and          *
 *                 Radu Serban @ LLNL                              *
 * Version of    : 28 March 2003                                   *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California * 
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/cvodes/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for the CVODES dense linear solver,     *
 * CVSDENSE.                                                       *
 *                                                                 *
 * Note: The type integertype must be large enough to store the    *
 * value of the linear system size N.                              *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvsdense_h
#define _cvsdense_h


#include <stdio.h>
#include "cvodes.h"
#include "sundialstypes.h"
#include "dense.h"
#include "nvector.h"

 
/******************************************************************
 *                                                                *
 * CVDENSE solver statistics indices                              *
 *----------------------------------------------------------------*
 * The following enumeration gives a symbolic name to each        *
 * CVDENSE statistic. The symbolic names are used as indices into *
 * the iopt and ropt arrays passed to CVodeMalloc.                *
 * The CVDENSE statistics are:                                    *
 *                                                                *
 * iopt[DENSE_NJE] : number of Jacobian evaluations, i.e. of      *
 *                   calls made to the dense Jacobian routine     *
 *                   (default or user-supplied).                  *
 *                                                                *
 * iopt[DENSE_LRW] : size (in realtype words) of real workspace   *
 *                   matrices and vectors used by this solver.    *
 *                                                                *
 * iopt[DENSE_LIW] : size (in integertype words) of integer       *
 *                   workspace vectors used by this solver.       *
 *                                                                *
 ******************************************************************/
 
enum { DENSE_NJE=CVODE_IOPT_SIZE, DENSE_LRW, DENSE_LIW };


/******************************************************************
 *                                                                *
 * CVDENSE solver constants                                       *
 *----------------------------------------------------------------*
 * CVD_MSBJ  : maximum number of steps between dense Jacobian     *
 *             evaluations                                        *
 *                                                                *
 * CVD_DGMAX : maximum change in gamma between dense Jacobian     *
 *             evaluations                                        *
 *                                                                *
 ******************************************************************/

#define CVD_MSBJ  50   

#define CVD_DGMAX RCONST(0.2)  

 
/******************************************************************
 *                                                                *           
 * Type : CVDenseJacFn                                            *
 *----------------------------------------------------------------*
 * A dense Jacobian approximation function Jac must have the      *
 * prototype given below. Its parameters are:                     *
 *                                                                *
 * N is the length of all vector arguments.                       *
 *                                                                *
 * J is the dense matrix (of type DenseMat) that will be loaded   *
 * by a CVDenseJacFn with an approximation to the Jacobian matrix *
 * J = (df_i/dy_j) at the point (t,y).                            *
 * J is preset to zero, so only the nonzero elements need to be   *
 * loaded. Two efficient ways to load J are:                      *
 *                                                                *
 * (1) (with macros - no explicit data structure references)      *
 *     for (j=0; j < n; j++) {                                    *
 *       col_j = DENSE_COL(J,j);                                  *
 *       for (i=0; i < n; i++) {                                  *
 *         generate J_ij = the (i,j)th Jacobian element           *
 *         col_j[i] = J_ij;                                       *
 *       }                                                        *
 *     }                                                          *
 *                                                                *  
 * (2) (without macros - explicit data structure references)      *
 *     for (j=0; j < n; j++) {                                    *
 *       col_j = (J->data)[j];                                    *
 *       for (i=0; i < n; i++) {                                  *
 *         generate J_ij = the (i,j)th Jacobian element           *
 *         col_j[i] = J_ij;                                       *
 *       }                                                        *
 *     }                                                          *
 *                                                                *
 * The DENSE_ELEM(A,i,j) macro is appropriate for use in small    *
 * problems in which efficiency of access is NOT a major concern. *
 *                                                                *
 * t is the current value of the independent variable.            *
 *                                                                *
 * y is the current value of the dependent variable vector,       *
 *      namely the predicted value of y(t).                       *
 *                                                                *
 * fy is the vector f(t,y).                                       *
 *                                                                *
 * jac_data is a pointer to user data - the same as the jac_data  *
 *          parameter passed to CVDense.                          *
 *                                                                *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for      *
 * vectors of length N which can be used by a CVDenseJacFn        *
 * as temporary storage or work space.                            *
 *                                                                *
 ******************************************************************/
  
typedef void (*CVDenseJacFn)(integertype n, DenseMat J, realtype t, 
                             N_Vector y, N_Vector fy, void *jac_data,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
 
 
/******************************************************************
 *                                                                *
 * Function : CVDense                                             *
 *----------------------------------------------------------------*
 * A call to the CVDense function links the main CVODE integrator *
 * with the CVDENSE linear solver.                                *
 *                                                                *
 * cvode_mem is the pointer to CVODE memory returned by           *
 *              CVodeMalloc.                                      *
 *                                                                *
 * n is the length of all vector arguments.                       *
 *                                                                *
 * djac is the dense Jacobian approximation routine to be used.   *
 *         A user-supplied djac routine must be of type           *
 *         CVDenseJacFn. Pass NULL for djac to use the default    *
 *         difference quotient routine CVDenseDQJac supplied      *
 *         with this solver.                                      *
 *                                                                *
 * jac_data is a pointer to user data which is passed to the      *
 *         djac routine every time it is called.                  *
 *                                                                *
 * The return values of CVDense are:                              *
 *    SUCCESS   = 0  if successful                                *
 *    LMEM_FAIL = -1 if there was a memory allocation failure     *
 *                                                                *
 ******************************************************************/
  
int CVDense(void *cvode_mem, integertype n, 
            CVDenseJacFn djac, void *jac_data);

/******************************************************************
 *                                                                *
 * Function : CVReInitDense                                       *
 *----------------------------------------------------------------*
 * A call to the CVReInitDense function resets the link between   *
 * the main CVODE integrator and the CVDENSE linear solver.       *
 * After solving one problem using CVDENSE, call CVReInit and then*
 * CVReInitDense to solve another problem of the same size, if    *
 * there is a change in the CVDense parameters djac or jac_data.  *
 * If there is no change in parameters, it is not necessary to    *
 * call either CVReInitDense or CVDense for the new problem.      *
 *                                                                *
 * All arguments to CVReInitDense have the same names and meanings*
 * as those of CVDense.  The cvode_mem argument must be identical *
 * to its value in the previous CVDense call.                     *
 *                                                                *
 * The return values of CVReInitDense are:                        *
 *   SUCCESS   = 0      if successful, or                         *
 *   LMEM_FAIL = -1     if the cvode_mem argument is NULL         *
 *                                                                *
 * NOTE: CVReInitDense performs the same compatibility tests as   *
 *       CVDense.                                                 *
 *                                                                *
 ******************************************************************/
  
int CVReInitDense(void *cvode_mem, CVDenseJacFn djac, void *jac_data);

#endif

#ifdef __cplusplus
}
#endif
