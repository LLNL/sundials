/******************************************************************
 *                                                                *
 * File          : sensidaband.h                                  *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh            *
 * Version of    : 3 July 2001                                    *
 *----------------------------------------------------------------*
 *                                                                *
 * This is the header file for the SensIDA band linear solver     *
 * module, SENSIDABAND. It interfaces between the band module and *
 * the IDA package when a banded linear solver is appropriate.    *
 *                                                                *
 * This is the header file for the IDA band linear solver,        *
 * IDABAND, for sensitivity analysis. The linear system size is   *
 * Ny, where Ny is the number of equations contained in           *
 * F(t,y,y',p) = 0.                                               *
 *                                                                *
 * This file is simply the idaband.h header file, along with a    *
 * function prototype for SensIDABand().                          *
 *                                                                *
 * Note: The type integer must be large enough to store the value *
 * of the linear system size Ny.                                  *
 *                                                                *
 ******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _sensidaband_h
#define _sensidaband_h


#include <stdio.h>
#include "ida.h"
#include "llnltyps.h"
#include "band.h"
#include "nvector.h"
#include "idaband.h"

/* Return values for SensIDABand: */

/* SUCCESS = 0 (defined in ida.h) */
enum {SensIDA_BAND_FAIL = -1, SensIDA_BAND_BAD_ARG = -2};


/******************************************************************
 *                                                                *
 * Function : SensIDABand                                         *
 *----------------------------------------------------------------*
 * A call to the SensIDABand function links the main IDA          *
 * integrator with the IDABAND linear solver module.              *
 *                                                                *
 * IDA_mem is the pointer to IDA memory returned by SensIDAMalloc.*
 *                                                                *
 * mupper is the upper bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * mlower is the lower bandwidth of the banded Jacobian matrix.   *
 *                                                                *
 * bjac is the banded Jacobian approximation routine to be used.  *
 *         A user-supplied bjac routine must be of type           *
 *         IDABandJacFn (see above).  Pass NULL for bjac if IDA   *
 *         is to use the default difference quotient routine      *
 *         IDABandDQJac supplied with this module.                *
 *                                                                *
 * jdata is a pointer to user data which is passed to the bjac    *
 *         routine every time it is called.                       *
 *                                                                *
 * IDABand returns either                                         *
 *    SUCCESS = 0            if successful, or                    *
 *    IDA_BAND_FAIL = -1     if either IDA_mem was NULL or a      *
 *                           malloc failure occurred, or          *
 *    IDA_BAND_BAD_ARG = -2  if mupper or mlower is illegal.      *
 ******************************************************************/

int SensIDABand(void *IDA_mem, integer mupper, integer mlower,
		IDABandJacFn bjac, void *jdata);

#endif

#ifdef __cplusplus
}
#endif
