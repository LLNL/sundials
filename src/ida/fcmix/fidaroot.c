/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2007-04-23 23:37:21 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Alan C. Hindmarsh @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * The FIDAROOT module contains the routines necessary to use
 * the rootfinding feature of the IDA module and to interface
 * with the user-supplied Fortran subroutine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"       /* actual function names, prototypes and global vars.*/
#include "fidaroot.h"   /* prototypes of interfaces to IDA                   */
#include "ida_impl.h"   /* definition of IDAMeme type                        */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FIDA_ROOTFN(realtype*,  /* T    */ 
                          realtype*,  /* Y    */
                          realtype*,  /* YP   */
                          realtype*,  /* G    */
                          long int*,  /* IPAR */
                          realtype*,  /* RPAR */
                          int*);      /* IER  */
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FIDA_ROOTINIT(int *nrtfn, int *ier)
{
  *ier = IDARootInit(IDA_idamem, *nrtfn, (IDARootFn) FIDArootfunc);
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
  IDARootInit(IDA_idamem, 0, NULL);

  return;
}

/***************************************************************************/

int FIDArootfunc(realtype t, N_Vector y, N_Vector yp, realtype *gout,
                 void *res_data)
{
  int ier;
  realtype *ydata, *ypdata;
  FIDAUserData IDA_userdata;

  ydata  = N_VGetArrayPointer(y);
  ypdata = N_VGetArrayPointer(yp);

  IDA_userdata = (FIDAUserData) res_data;

  FIDA_ROOTFN(&t, ydata, ypdata, gout, IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}
