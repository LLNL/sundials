/*
 * -----------------------------------------------------------------
 * $Revision:  $
 * $Date:  $
 * -----------------------------------------------------------------
 * Programmer: Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This header file contains definitions and declarations for use by
 * sparse linear solvers for Ax = b. 
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_SPARSE_H
#define _SUNDIALS_SPARSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_types.h>

/*
 * ==================================================================
 * Type definitions
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : SlsMat
 * -----------------------------------------------------------------
 * The type SlsMat is defined to be a pointer to a structure
 * with various sizes, a data field, and arrays for the row and 
 * column information for the sparse matrix entries.
 * The M and N fields indicates the number 
 * of rows and columns, respectively. The data field is a one 
 * dimensional array used for component storage. The NNZ field indicates
 * the number of nonzero entries in the matrix. The integer array, asub, 
 * holds the row index for each of the matrix entries.  The integer
 * array, xa, holds the index entry for the starting value of each column.
 * -----------------------------------------------------------------
 * The relevant fields in DlsMat are:
 *    M     - number of rows
 *    N     - number of columns
 *    NNZ   - the number of nonzero entries in the matrix
 *    data  - pointer to a contiguous block of realtype variables
 *    rowvals - row indices of each nonzero entry
 *    colptrs - starting index of the first entry in data in each column
 *
 * The nonzero entries of the matrix are stored in
 * compressed column format.  Row indices of entries in 
 * column j are stored in rowvals[colptrs[j]] through rowvals[colptrs[j+i]-1]
 * and corresponding numerical values of the matrix are stored 
 * in the same entries of data.
 * -----------------------------------------------------------------
 */

typedef struct _SlsMat {
  int M;
  int N;
  int NNZ;
  realtype *data;
  int *rowvals;
  int *colptrs;
} *SlsMat;

/*
 * ==================================================================
 * Exported function prototypes (functions working on SlsMat)
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function: NewSparseMat
 * -----------------------------------------------------------------
 * NewSparseMat allocates memory for a compressed column sparse
 * matrix with M rows, N columns, and NNZ nonzeros. NewSparseMat
 * returns NULL if the request for matrix storage cannot be
 * satisfied. See the above documentation for the type SlsMat
 * for matrix storage details.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT SlsMat NewSparseMat(int M, int N, int NNZ);

/*
 * -----------------------------------------------------------------
 * Functions: DestroySparseMat
 * -----------------------------------------------------------------
 * DestroySparseMat frees the memory allocated by NewSparseMat
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void DestroySparseMat(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Function : SlsSetToZero
 * -----------------------------------------------------------------
 * SetToZero sets all the elements of the sparse matrix A to 0.0.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void SlsSetToZero(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: PrintSparseMat
 * -----------------------------------------------------------------
 * This function prints the compressed column matrix information for 
 * matrix A to standard output.
 * It is intended as a debugging tool with small values of NNZ.
 * The elements are printed using the %g/%lg/%Lg option. 
 * A blank line is printed before and after the matrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void PrintSparseMat(SlsMat A);


#ifdef __cplusplus
}
#endif

#endif
