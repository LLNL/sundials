/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2004-07-22 21:21:16 $
 * ----------------------------------------------------------------- 
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, and         
 *              Radu Serban @ LLNL                              
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California 
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvodes/LICENSE
 * -----------------------------------------------------------------
 * Implementation header file for the diagonal linear solver, CVDIAG. 
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _cvdiag_impl_h
#define _cvdiag_impl_h

#include <stdio.h>

#include "cvdiag.h"

#include "nvector.h"
#include "sundialstypes.h"

/*
 * -----------------------------------------------------------------
 * Types: CVDiagMemRec, CVDiagMem
 * -----------------------------------------------------------------
 * The type CVDiagMem is pointer to a CVDiagMemRec.
 * This structure contains CVDiag solver-specific data.
 * -----------------------------------------------------------------
 */

typedef struct {

  realtype di_gammasv; /* gammasv = gamma at the last call to setup */
                       /* or solve                                  */

  N_Vector di_M;       /* M = (I - gamma J)^{-1} , gamma = h / l1   */

  N_Vector di_bit;     /* temporary storage vector                  */

  N_Vector di_bitcomp; /* temporary storage vector                  */

  long int di_nfeDI;   /* no. of calls to f due to difference 
                          quotient diagonal Jacobian approximation  */

} CVDiagMemRec, *CVDiagMem;

#endif

#ifdef __cplusplus
}
#endif
