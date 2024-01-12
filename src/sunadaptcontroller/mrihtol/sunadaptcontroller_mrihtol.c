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
 * This is the implementation file for the SUNControl_MRIHTol module.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <suncontrol/suncontrol_mrihtol.h>
#include <sundials/sundials_math.h>


/* ---------------
 * Macro accessors
 * --------------- */

#define SC_MRIHTOL_CONTENT(C) ( (SUNControlContent_MRIHTol)(C->content) )
#define SC_MRIHTOL_CSLOW(C)   ( SC_MRIHTOL_CONTENT(C)->HControl )
#define SC_MRIHTOL_CFAST(C)   ( SC_MRIHTOL_CONTENT(C)->TolControl )


/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/* -----------------------------------------------------------------
 * Function to create a new MRIHTol controller
 */

SUNControl SUNControlMRIHTol(SUNContext sunctx, SUNControl HControl,
                             SUNControl TolControl)
{
  SUNControl C;
  SUNControlContent_MRIHTol content;

  /* Create an empty controller object */
  C = NULL;
  C = SUNControlNewEmpty(sunctx);
  if (C == NULL) { return (NULL); }

  /* Attach operations */
  C->ops->getid             = SUNControlGetID_MRIHTol;
  C->ops->estimatesteptol   = SUNControlEstimateStepTol_MRIHTol;
  C->ops->reset             = SUNControlReset_MRIHTol;
  C->ops->setdefaults       = SUNControlSetDefaults_MRIHTol;
  C->ops->write             = SUNControlWrite_MRIHTol;
  C->ops->setmethodorder    = SUNControlSetMethodOrder_MRIHTol;
  C->ops->setembeddingorder = SUNControlSetEmbeddingOrder_MRIHTol;
  C->ops->seterrorbias      = SUNControlSetErrorBias_MRIHTol;
  C->ops->updatemritol      = SUNControlUpdateMRITol_MRIHTol;
  C->ops->space             = SUNControlSpace_MRIHTol;

  /* Create content */
  content = NULL;
  content = (SUNControlContent_MRIHTol)malloc(sizeof *content);
  if (content == NULL)
  {
    SUNControlDestroy(C);
    return (NULL);
  }

  /* Set method and embedding "order" into tolerance controller:
     no matter the order of the fast integrator, we expect the
     fast integrator error to be directly proportional to the
     relative tolerance factor */
  SUNControlSetMethodOrder(TolControl, 1);
  SUNControlSetEmbeddingOrder(TolControl, 1);

  /* Attach input controllers */
  content->HControl = HControl;
  content->TolControl = TolControl;

  /* Attach content */
  C->content = content;

  return (C);
}

/* -----------------------------------------------------------------
 * Function to get slow and fast sub-controllers
 */

SUNControl SUNControlMRIHTol_GetSlowController(SUNControl C)
{ return SC_MRIHTOL_CSLOW(C); }

SUNControl SUNControlMRIHTol_GetFastController(SUNControl C)
{ return SC_MRIHTOL_CFAST(C); }



/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

SUNControl_ID SUNControlGetID_MRIHTol(SUNControl C) { return SUNDIALS_CONTROL_MRI_TOL; }

int SUNControlEstimateStepTol_MRIHTol(SUNControl C, realtype H,
                                         realtype tolfac, realtype DSM,
                                         realtype dsm, realtype *Hnew,
                                         realtype* tolfacnew)
{
  /* Call sub-controllers to fill outputs, and return with success */
  int retval;
  retval = SUNControlEstimateStep(SC_MRIHTOL_CSLOW(C), H, DSM, Hnew);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  retval = SUNControlEstimateStep(SC_MRIHTOL_CFAST(C), tolfac, dsm, tolfacnew);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlReset_MRIHTol(SUNControl C)
{
  /* Reset sub-controllers, and return with success */
  int retval;
  retval = SUNControlReset(SC_MRIHTOL_CSLOW(C));
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  retval = SUNControlReset(SC_MRIHTOL_CFAST(C));
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetDefaults_MRIHTol(SUNControl C)
{
  /* Set defaults in sub-controllers, and return with success */
  int retval;
  retval = SUNControlSetDefaults(SC_MRIHTOL_CSLOW(C));
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  retval = SUNControlSetDefaults(SC_MRIHTOL_CFAST(C));
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlWrite_MRIHTol(SUNControl C, FILE *fptr)
{
  /* Write both sub-controllers, and return with success */
  fprintf(fptr, "Multirate H-Tol SUNControl module:\n");
  fprintf(fptr, "\nSlow step controller:\n");
  int retval= SUNControlWrite(SC_MRIHTOL_CSLOW(C), fptr);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  fprintf(fptr, "\nFast tolerance controller:\n");
  retval= SUNControlWrite(SC_MRIHTOL_CFAST(C), fptr);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSetMethodOrder_MRIHTol(SUNControl C, int q)
{
  /* Send method order to slow controller */
  return (SUNControlSetMethodOrder(SC_MRIHTOL_CSLOW(C), q));
}

int SUNControlSetEmbeddingOrder_MRIHTol(SUNControl C, int p)
{
  /* Send method order to slow controller */
  return (SUNControlSetEmbeddingOrder(SC_MRIHTOL_CSLOW(C), p));
}

int SUNControlSetErrorBias_MRIHTol(SUNControl C, realtype bias)
{
  /* Set defaults in sub-controllers, and return with success */
  int retval;
  retval = SUNControlSetErrorBias(SC_MRIHTOL_CSLOW(C), bias);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  retval = SUNControlSetErrorBias(SC_MRIHTOL_CFAST(C), bias);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlUpdateMRITol_MRIHTol(SUNControl C, realtype H,
                                   realtype tolfac, realtype DSM,
                                   realtype dsm)
{
  /* Update sub-controllers, and return with success */
  int retval;
  retval = SUNControlUpdate(SC_MRIHTOL_CSLOW(C), H, DSM);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  retval = SUNControlUpdate(SC_MRIHTOL_CFAST(C), tolfac, dsm);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  return SUNCONTROL_SUCCESS;
}

int SUNControlSpace_MRIHTol(SUNControl C, long int* lenrw, long int* leniw)
{
  /* Accumulate space from sub-controllers, and return with success */
  int retval;
  long int lrw, liw;
  retval = SUNControlSpace(SC_MRIHTOL_CSLOW(C), lenrw, leniw);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  retval = SUNControlSpace(SC_MRIHTOL_CFAST(C), &lrw, &liw);
  if (retval != SUNCONTROL_SUCCESS) { return retval; }
  *lenrw += lrw;
  *leniw += liw;
  return SUNCONTROL_SUCCESS;
}
