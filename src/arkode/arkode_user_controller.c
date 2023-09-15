/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the ARKUserControl
 * SUNControl module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "arkode_user_controller.h"


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_CONTENT(C)     ( (ARKUserControlContent)(C->content) )
#define SC_HP(C)          ( SC_CONTENT(C)->hp )
#define SC_HPP(C)         ( SC_CONTENT(C)->hpp )
#define SC_EP(C)          ( SC_CONTENT(C)->ep )
#define SC_EPP(C)         ( SC_CONTENT(C)->epp )
#define SC_P(C)           ( SC_CONTENT(C)->p )
#define SC_Q(C)           ( SC_CONTENT(C)->q )
#define SC_ARKMEM(C)      ( SC_CONTENT(C)->ark_mem )
#define SC_HADAPT(C)      ( SC_CONTENT(C)->hadapt )
#define SC_DATA(C)        ( SC_CONTENT(C)->hadapt_data )


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new ARKUserControl controller
 */

SUNControl ARKUserControl(SUNContext sunctx, void* arkode_mem,
                          ARKAdaptFn hadapt, void* hadapt_data)
{
  SUNControl C;
  ARKUserControlContent content;

  /* Return with failure if hadapt or arkode_mem are NULL */
  if ((hadapt == NULL) || (arkode_mem == NULL)) { return (NULL); }

  /* Create an empty controller object */
  C = NULL;
  C = SUNControl_NewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->gettype           = SUNControlGetType_ARKUserControl;
  C->ops->estimatestep      = SUNControlEstimateStep_ARKUserControl;
  C->ops->reset             = SUNControlReset_ARKUserControl;
  C->ops->write             = SUNControlWrite_ARKUserControl;
  C->ops->setmethodorder    = SUNControlSetMethodOrder_ARKUserControl;
  C->ops->setembeddingorder = SUNControlSetEmbeddingOrder_ARKUserControl;
  C->ops->update            = SUNControlUpdate_ARKUserControl;
  C->ops->space             = SUNControlSpace_ARKUserControl;

  /* Create content */
  content = NULL;
  content = (ARKUserControlContent)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNControl_Destroy(C);
    return (NULL);
  }

  /* Attach content */
  C->content = content;

  /* Attach ARKODE memory structure */
  content->ark_mem = (ARKodeMem) arkode_mem;

  /* Attach user-provided adaptivity function and data */
  content->hadapt = hadapt;
  content->hadapt_data = hadapt_data;

  /* Initialize method and embedding orders */
  content->p = 1;
  content->q = 1;

  /* Fill content with default/reset values */
  SUNControlReset_ARKUserControl(C);

  return (C);
}


/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_Type SUNControlGetType_ARKUserControl(SUNControl C)
{ return SUNDIALS_CONTROL_H; }

int SUNControlEstimateStep_ARKUserControl(SUNControl C, realtype h,
                                          realtype dsm, realtype* hnew)
{
  /* call user-provided function to compute new step */
  int retval = SC_HADAPT(C)(SC_ARKMEM(C)->ycur, SC_ARKMEM(C)->tn, h, SC_HP(C),
                            SC_HPP(C), dsm, SC_EP(C), SC_EPP(C), SC_Q(C),
                            SC_P(C), hnew, SC_DATA(C));
  if (retval != SUNCONTROL_SUCCESS) { return(SUNCONTROL_USER_FCN_FAIL); }
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_ARKUserControl(SUNControl C)
{
  SC_EP(C)  = RCONST(1.0);
  SC_EPP(C) = RCONST(1.0);
  SC_HP(C)  = RCONST(0.0);
  SC_HPP(C) = RCONST(0.0);
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_ARKUserControl(SUNControl C, FILE *fptr)
{
  fprintf(fptr, "ARKUserControl module:\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fptr, "  hp = %12Lg\n", SC_HP(C));
  fprintf(fptr, "  hpp = %12Lg\n", SC_HPP(C));
  fprintf(fptr, "  ep = %12Lg\n", SC_EP(C));
  fprintf(fptr, "  epp = %12Lg\n", SC_EPP(C));
#else
  fprintf(fptr, "  hp = %12g\n", SC_HP(C));
  fprintf(fptr, "  hpp = %12g\n", SC_HPP(C));
  fprintf(fptr, "  ep = %12g\n", SC_EP(C));
  fprintf(fptr, "  epp = %12g\n", SC_EPP(C));
#endif
  fprintf(fptr, "  p = %i\n", SC_P(C));
  fprintf(fptr, "  q = %i\n", SC_Q(C));
  fprintf(fptr, "  hadapt_data = %p\n", SC_DATA(C));
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_ARKUserControl(SUNControl C, int q)
{
  /* check for legal input */
  if (q <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store value and return */
  SC_Q(C) = q;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetEmbeddingOrder_ARKUserControl(SUNControl C, int p)
{
  /* check for legal input */
  if (p <= 0) { return SUNCONTROL_ILL_INPUT; }

  /* store value and return */
  SC_P(C) = p;
  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdate_ARKUserControl(SUNControl C, realtype h, realtype dsm)
{
  SC_HPP(C) = SC_HP(C);
  SC_HP(C) = h;
  SC_EPP(C) = SC_EP(C);
  SC_EP(C) = dsm;
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_ARKUserControl(SUNControl C, long int* lenrw,
                                   long int* leniw)
{
  *lenrw = 4;
  *leniw = 5;
  return SUNCONTROL_SUCCESS;
}
