/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006-01-11 21:13:53 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Alan C. Hindmarsh @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * The FIDAROOT module contains the routines necessary to use
 * the rootfinding feature of the IDA module and to interface
 * with the user-supplied Fortran subroutine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida.h"              /* IDA constants and prototypes          */
#include "fida.h"             /* actual function names, prototypes and
			         global variables                      */
#include "fidaroot.h"         /* prototypes of interfaces to IDA       */
#include "sundials_nvector.h" /* definition of type N_Vector           */
#include "sundials_types.h"   /* definition of SUNDIALS type realtype  */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FIDA_ROOTFN(realtype *, realtype*, realtype*, realtype*);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FIDA_ROOTINIT(int *nrtfn, int *ier)
{
  *ier = IDARootInit(IDA_idamem, *nrtfn, (IDARootFn) FIDArootfunc, NULL);
  IDA_nrtfn = *nrtfn;

  return;
}

/***************************************************************************/

void FIDA_ROOTINFO(int *nrtfn, int *info, int *ier)
{
  *ier = IDAGetRootInfo(IDA_idamem, info);
  return; 
}

/***************************************************************************/

void FIDA_ROOTFREE(void)
{
  IDARootInit(IDA_idamem, 0, NULL, NULL);

  return;
}

/***************************************************************************/

void FIDArootfunc(realtype t, N_Vector y, N_Vector yp, realtype *gout,
                  void *g_data)
{
  realtype *ydata, *ypdata;

  ydata  = N_VGetArrayPointer(y);
  ypdata = N_VGetArrayPointer(yp);

  FIDA_ROOTFN(&t, ydata, ypdata, gout);

  return;
}
