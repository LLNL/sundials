/******************************************************************
 *                                                                *
 * File          : sensidadense.h                                 *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh            *
 * Version of    : 3 July 2001                                    *
 *----------------------------------------------------------------*
 *                                                                *
 * This is the header file for the IDA band linear solver,        *
 * IDABAND, for sensitivity analysis. The linear system size is   *
 * Ny, where Ny is the number of equations contained in           *
 * F(t,y,y',p) = 0.                                               *
 *                                                                *
 * This file is simply the idadense.h header file, along with a   *
 * function prototype for SensIDADense().                         *
 *                                                                *
 * Note: The type integer must be large enough to store the value *
 * of the linear system size Ny.                                  *
 *                                                                *
 ******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _sensidadense_h
#define _sensidadense_h


#include <stdio.h>
#include "ida.h"
#include "llnltyps.h"
#include "dense.h"
#include "nvector.h"
#include "idadense.h"

/* Return values for SensIDADense: */

/* SUCCESS = 0 (defined in ida.h) */
enum {SensIDA_DENSE_FAIL = -1};

 
/******************************************************************
 *                                                                *
 * Function : SensIDADense                                        *
 *----------------------------------------------------------------*
 * A call to the SensIDADense function links the main IDA         *
 * integrator with the SENSIDADENSE linear solver module.         *
 *                                                                *
 * IDA_mem is the pointer to IDA memory returned by SensIDAMalloc.*
 *                                                                *
 * djac is the dense Jacobian approximation routine to be used.   *
 *         A user-supplied djac routine must be of type           *
 *         IDADenseJacFn (see above).  Pass NULL for djac if IDA  *
 *         is to use the default difference quotient routine      *
 *         IDADenseDQJac supplied with this module.               *
 *                                                                *
 * jdata is a pointer to user data which is passed to the djac    *
 *         routine every time it is called.                       *
 *                                                                *
 * SensIDADense returns either                                    *
 *     SUCCESS = 0             if successful, or                  *
 *     IDA_DENSE_FAIL = -1     if either IDA_mem was null, or a   *
 *                             malloc failure occurred.           *
 ******************************************************************/

int SensIDADense(void *IDA_mem, IDADenseJacFn djac, void *jdata);

#endif

#ifdef __cplusplus
}
#endif
