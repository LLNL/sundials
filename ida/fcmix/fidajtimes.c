/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2006-02-02 00:34:31 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/ida/LICENSE.
 * -----------------------------------------------------------------
 * The C function FIDAJtimes is to interface between the
 * IDASPILS modules and the user-supplied Jacobian-vector
 * product routine FIDAJTIMES. Note the use of the generic name
 * FIDA_JTIMES below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_spils.h"
#include "fida.h"             /* actual fn. names, prototypes and global vars.  */
#include "ida_impl.h"         /* definition of IDAMem type                      */
#include "sundials_nvector.h" /* definitions of type N_Vector and vector macros */
#include "sundials_types.h"   /* definition of type realtype                    */

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_JTIMES(realtype*, realtype*, realtype*,     /* T, Y, YP   */
                          realtype*, realtype*, realtype*,     /* R, V, FJV  */
                          realtype*, realtype*, realtype*,     /* CJ, EWT, H */
                          long int*, realtype*,                /* IPAR, RPAR */
                          realtype*, realtype*,                /* WK1, WK2   */
                          int*);                               /* IER        */

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_SPILSSETJAC(int *flag, int *ier)
{
  IDAMem ida_mem;
  *ier = 0;

  if (*flag == 0) {

    *ier = IDASpilsSetJacTimesVecFn(IDA_idamem, NULL, NULL);

  } else {

    if (F2C_IDA_ewtvec == NULL) {
      F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
      if (F2C_IDA_ewtvec == NULL) {
        *ier = -1;
        return;
      }
    }

    ida_mem = (IDAMem) IDA_idamem;
    *ier = IDASpilsSetJacTimesVecFn(IDA_idamem,
                                   (IDASpilsJacTimesVecFn) FIDAJtimes, ida_mem->ida_rdata);
  }

  return;
}

/*************************************************/

int FIDAJtimes(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	       N_Vector v, N_Vector Jv,
	       realtype c_j, void *jac_data,
	       N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *yy_data, *yp_data, *rr_data, *vdata, *Jvdata, *ewtdata;
  realtype *v1data, *v2data;
  realtype h;
  FIDAUserData IDA_userdata;
  int ier;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = vdata = Jvdata = ewtdata = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);
  IDAGetLastStep(IDA_idamem, &h);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  vdata   = N_VGetArrayPointer(v);
  Jvdata  = N_VGetArrayPointer(Jv);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);

  IDA_userdata = (FIDAUserData) jac_data;

  /* Call user-supplied routine */
  FIDA_JTIMES(&t, yy_data, yp_data, rr_data, vdata, Jvdata,
	      &c_j, ewtdata, &h, 
              IDA_userdata->ipar, IDA_userdata->rpar,
              v1data, v2data, &ier);

  return(ier);
}
