/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the diagonal implementation of the 
 * SUNMATRIX module.
 * 
 * Part I contains declarations specific to the diagonal 
 * implementation of the supplied SUNMATRIX module.
 * 
 * Part II defines accessor macros that allow the user to 
 * efficiently use this SUNMatrix type without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor 
 * SUNDiagonalMatrix as well as implementation-specific prototypes 
 * for various useful matrix operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNMatrix structure can be found
 *     in the header file sundials_matrix.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype' and 'indextype'.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_DIAGONAL_H
#define _SUNMATRIX_DIAGONAL_H

#include <sundials/sundials_matrix.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: Diagonal implementation of SUNMatrix
 *
 * The diagonal implementation of the SUNMatrix 'content' structure
 * contains:
 *   D     - N_Vector containing the matrix diagonal
 * -----------------------------------------------------------------
 */
  
struct _SUNMatrixContent_Diagonal {
  N_Vector d;
};

typedef struct _SUNMatrixContent_Diagonal *SUNMatrixContent_Diagonal;

/*
 * -----------------------------------------------------------------
 * PART II: macros SM_CONTENT_DIAG, SM_DATA_DIAG
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * SUNMatrix A;
 * SUNMatrixContent_Diagonal A_cont;
 * N_Vector A_diag;
 *
 * (1) SM_CONTENT_DIAG
 *
 *     This macro gives access to the contents of the diagonal
 *     SUNMatrix.
 *
 *     The assignment A_cont = SM_CONTENT_DIAG(A) sets A_cont to 
 *     be a pointer to the diagonal SUNMatrix content structure.
 *
 * (2) SM_DATA_DIAG
 *
 *     This macro gives access to the N_Vector containing the 
 *     matrix diagonal.
 *
 *     The assignment A_diag = SM_DATA_DIAG(A) sets A_diag to be
 *     the N_Vector containing the matrix diagonal.
 *
 * -----------------------------------------------------------------
 */

#define SM_CONTENT_DIAG(A)  ( (SUNMatrixContent_Diagonal)(A->content) )

#define SM_DATA_DIAG(A)     ( SM_CONTENT_DIAG(A)->d )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunmatrix_diagonal
 * 
 * CONSTRUCTORS:
 *    SUNDiagonalMatrix
 * OTHER:
 *    SUNDiagonalMatrix_Diag 
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function: SUNDiagonalMatrix
 * -----------------------------------------------------------------
 * Creates and allocates memory for an diagonal SUNMatrix, using 
 * the input as a template N_Vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNMatrix SUNDiagonalMatrix(N_Vector tmpl);

/*
 * -----------------------------------------------------------------
 * Accessor Functions: 
 *
 * SUNDiagonalMatrix_Diag
 *    Returns an N_Vector containing the matrix diagonal
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector SUNDiagonalMatrix_Diag(SUNMatrix A);

/*
 * -----------------------------------------------------------------
 * diagonal implementations of various useful matrix operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_Diagonal(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_Diagonal(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_Diagonal(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_Diagonal(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_Diagonal(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_Diagonal(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_Diagonal(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_Diagonal(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_Diagonal(SUNMatrix A, long int *lenrw,
                                         long int *leniw);

  
#ifdef __cplusplus
}
#endif

#endif
