/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for ARKODE's time step
 * adaptivity utilities.
 *--------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"

/*---------------------------------------------------------------
  arkAdaptInit:

  This routine creates and sets default values in an
  ARKodeHAdaptMem structure.  This returns a non-NULL structure
  if no errors occurred, or a NULL value otherwise.
  ---------------------------------------------------------------*/
ARKodeHAdaptMem arkAdaptInit(void)
{
  ARKodeHAdaptMem hadapt_mem;

  /* allocate structure */
  hadapt_mem = (ARKodeHAdaptMem)malloc(sizeof(struct ARKodeHAdaptMemRec));
  if (hadapt_mem == NULL) { return (NULL); }

  /* initialize values (default parameters are set in arkSetDefaults) */
  memset(hadapt_mem, 0, sizeof(struct ARKodeHAdaptMemRec));
  hadapt_mem->nst_acc = 0;
  hadapt_mem->nst_exp = 0;
  return (hadapt_mem);
}

/*---------------------------------------------------------------
  arkPrintAdaptMem

  This routine outputs the time step adaptivity memory structure
  to a specified file pointer.
  ---------------------------------------------------------------*/
void arkPrintAdaptMem(ARKodeHAdaptMem hadapt_mem, FILE* outfile)
{
  if (hadapt_mem != NULL)
  {
    fprintf(outfile, "ark_hadapt: etamax = %" RSYM "\n", hadapt_mem->etamax);
    fprintf(outfile, "ark_hadapt: etamx1 = %" RSYM "\n", hadapt_mem->etamx1);
    fprintf(outfile, "ark_hadapt: etamxf = %" RSYM "\n", hadapt_mem->etamxf);
    fprintf(outfile, "ark_hadapt: etamin = %" RSYM "\n", hadapt_mem->etamin);
    fprintf(outfile, "ark_hadapt: small_nef = %i\n", hadapt_mem->small_nef);
    fprintf(outfile, "ark_hadapt: etacf = %" RSYM "\n", hadapt_mem->etacf);
    fprintf(outfile, "ark_hadapt: cfl = %" RSYM "\n", hadapt_mem->cfl);
    fprintf(outfile, "ark_hadapt: safety = %" RSYM "\n", hadapt_mem->safety);
    fprintf(outfile, "ark_hadapt: growth = %" RSYM "\n", hadapt_mem->growth);
    fprintf(outfile, "ark_hadapt: lbound = %" RSYM "\n", hadapt_mem->lbound);
    fprintf(outfile, "ark_hadapt: ubound = %" RSYM "\n", hadapt_mem->ubound);
    fprintf(outfile, "ark_hadapt: nst_acc = %li\n", hadapt_mem->nst_acc);
    fprintf(outfile, "ark_hadapt: nst_exp = %li\n", hadapt_mem->nst_exp);
    fprintf(outfile, "ark_hadapt: pq = %i\n", hadapt_mem->pq);
    fprintf(outfile, "ark_hadapt: p = %i\n", hadapt_mem->p);
    fprintf(outfile, "ark_hadapt: q = %i\n", hadapt_mem->q);
    fprintf(outfile, "ark_hadapt: adjust = %i\n", hadapt_mem->adjust);
    if (hadapt_mem->expstab == arkExpStab)
    {
      fprintf(outfile, "  ark_hadapt: Default explicit stability function\n");
    }
    else
    {
      fprintf(outfile,
              "  ark_hadapt: User provided explicit stability function\n");
      fprintf(outfile, "  ark_hadapt: stability function data pointer = %p\n",
              hadapt_mem->estab_data);
    }
    (void)SUNAdaptController_Write(hadapt_mem->hcontroller, outfile);
  }
}

/*---------------------------------------------------------------
  arkAdapt is the time step adaptivity wrapper function.  This
  computes and sets the value of ark_eta inside of the ARKodeMem
  data structure.
  ---------------------------------------------------------------*/
int arkAdapt(void* arkode_mem, ARKodeHAdaptMem hadapt_mem, N_Vector ycur,
             sunrealtype tcur, sunrealtype hcur, sunrealtype dsm, long int nst)
{
  int retval;
  sunrealtype h_acc, h_cfl, int_dir;
  ARKodeMem ark_mem;
  int controller_order;
  if (arkode_mem == NULL)
  {
    arkProcessError(NULL, ARK_MEM_NULL, __LINE__, __func__, __FILE__,
                    MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem)arkode_mem;

  /* Request error-based step size from adaptivity controller */
  if (hadapt_mem->pq == 0)
  {
    controller_order = hadapt_mem->p + hadapt_mem->adjust;
  }
  else if (hadapt_mem->pq == 1)
  {
    controller_order = hadapt_mem->q + hadapt_mem->adjust;
  }
  else
  {
    controller_order = SUNMIN(hadapt_mem->p, hadapt_mem->q) + hadapt_mem->adjust;
  }
  retval = SUNAdaptController_EstimateStep(hadapt_mem->hcontroller, hcur,
                                           controller_order, dsm, &h_acc);
  if (retval != SUNADAPTCONTROLLER_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_CONTROLLER_ERR, __LINE__, __func__, __FILE__,
                    "SUNAdaptController_EstimateStep failure.");
    return (ARK_CONTROLLER_ERR);
  }

  /* determine direction of integration */
  int_dir = hcur / SUNRabs(hcur);

  /* Call explicit stability function */
  retval = hadapt_mem->expstab(ycur, tcur, &h_cfl, hadapt_mem->estab_data);
  if (retval != ARK_SUCCESS)
  {
    arkProcessError(ark_mem, ARK_ILL_INPUT, __LINE__, __func__, __FILE__,
                    "Error in explicit stability function.");
    return (ARK_ILL_INPUT);
  }
  if (h_cfl <= ZERO) { h_cfl = SUN_RCONST(1.0e30) * SUNRabs(hcur); }

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO, "ARKODE::arkAdapt",
                     "new-step-before-bounds",
                     "h_acc = %" RSYM ", h_cfl = %" RSYM, h_acc, h_cfl);
#endif

  /* enforce safety factors */
  h_acc *= hadapt_mem->safety;
  h_cfl *= hadapt_mem->cfl * int_dir;

  /* enforce maximum bound on time step growth */
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(hadapt_mem->etamax * hcur));

  /* enforce minimum bound time step reduction */
  h_acc = int_dir * SUNMAX(SUNRabs(h_acc), SUNRabs(hadapt_mem->etamin * hcur));

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_INFO, "ARKODE::arkAdapt",
                     "new-step-after-max-min-bounds",
                     "h_acc = %" RSYM ", h_cfl = %" RSYM, h_acc, h_cfl);
#endif

  /* increment the relevant step counter, set desired step */
  if (SUNRabs(h_acc) < SUNRabs(h_cfl)) { hadapt_mem->nst_acc++; }
  else { hadapt_mem->nst_exp++; }
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(h_cfl));

  /* enforce adaptivity bounds to retain Jacobian/preconditioner accuracy */
  if (dsm <= ONE)
  {
    if ((SUNRabs(h_acc) > SUNRabs(hcur * hadapt_mem->lbound * ONEMSM)) &&
        (SUNRabs(h_acc) < SUNRabs(hcur * hadapt_mem->ubound * ONEPSM)))
    {
      h_acc = hcur;
    }
  }

  /* set basic value of ark_eta */
  ark_mem->eta = h_acc / hcur;

  /* enforce minimum time step size */
  ark_mem->eta = SUNMAX(ark_mem->eta, ark_mem->hmin / SUNRabs(hcur));

  /* enforce maximum time step size */
  ark_mem->eta /= SUNMAX(ONE, SUNRabs(hcur) * ark_mem->hmax_inv * ark_mem->eta);

#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_DEBUG
  SUNLogger_QueueMsg(ARK_LOGGER, SUN_LOGLEVEL_DEBUG, "ARKODE::arkAdapt",
                     "new-step-eta", "eta = %" RSYM, ark_mem->eta);
#endif

  return (retval);
}

/*===============================================================
  EOF
  ===============================================================*/
