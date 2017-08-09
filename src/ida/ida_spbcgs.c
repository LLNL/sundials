/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the IDA scaled preconditioned
 * Bi-CGSTAB linear solver module, IDASPBCG.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <ida/ida_spbcgs.h>
#include "ida_spils_impl.h"
#include "ida_impl.h"

#include <sundials/sundials_spbcgs.h>
#include <sundials/sundials_math.h>

/* Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define PT9  RCONST(0.9)
#define PT05 RCONST(0.05)

/* IDASPBCG linit, lsetup, lsolve, lperf, and lfree routines */

static int IDASpbcgInit(IDAMem IDA_mem);

static int IDASpbcgSetup(IDAMem IDA_mem, 
                         N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int IDASpbcgSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector yy_now, N_Vector yp_now, N_Vector rr_now);

static int IDASpbcgPerf(IDAMem IDA_mem, int perftask);

static int IDASpbcgFree(IDAMem IDA_mem);


/*
 * -----------------------------------------------------------------
 * Function : IDASpbcg
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDASPBCG linear solver module.
 *
 * IDASpbcg first calls the existing lfree routine if this is not NULL.
 * It then sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDASpbcgInit, IDASpbcgSetup,
 * IDASpbcgSolve, IDASpbcgPerf, and IDASpbcgFree, respectively.
 * It allocates memory for a structure of type IDASpilsMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem). It then sets various fields
 * in the IDASpilsMemRec structure. Finally, IDASpbcg allocates memory
 * for ytemp, yptemp, and xx, and calls SpbcgMalloc to allocate memory
 * for the Spbcg solver.
 *
 * The return value of IDASpbcg is:
 *   IDASPILS_SUCCESS   =  0 if successful
 *   IDASPILS_MEM_FAIL  = -1 if IDA_mem is NULL or a memory
 *                           allocation failed
 *   IDASPILS_ILL_INPUT = -2 if a required vector operation is not
 *                           implemented.
 * -----------------------------------------------------------------
 */

int IDASpbcg(void *ida_mem, int maxl)
{
  IDAMem IDA_mem;
  IDASpilsMem idaspils_mem;
  SpbcgMem spbcg_mem;
  int maxl1;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_NULL, "IDASPBCG", "IDASpbcg", MSGS_IDAMEM_NULL);
    return(IDASPILS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Check if N_VDotProd is present */
  if (IDA_mem->ida_tempv1->ops->nvdotprod == NULL) {
    IDAProcessError(NULL, IDASPILS_ILL_INPUT, "IDASPBCG", "IDASpbcg", MSGS_BAD_NVECTOR);
    return(IDASPILS_ILL_INPUT);
  }

  if (IDA_mem->ida_lfree != NULL) IDA_mem->ida_lfree((IDAMem) ida_mem);

  /* Set five main function fields in ida_mem */
  IDA_mem->ida_linit  = IDASpbcgInit;
  IDA_mem->ida_lsetup = IDASpbcgSetup;
  IDA_mem->ida_lsolve = IDASpbcgSolve;
  IDA_mem->ida_lperf  = IDASpbcgPerf;
  IDA_mem->ida_lfree  = IDASpbcgFree;

  /* Get memory for IDASpilsMemRec */
  idaspils_mem = NULL;
  idaspils_mem = (IDASpilsMem) malloc(sizeof(struct IDASpilsMemRec));
  if (idaspils_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPBCG", "IDASpbcg", MSGS_MEM_FAIL);
    return(IDASPILS_MEM_FAIL);
  }

  /* Set ILS type */
  idaspils_mem->s_type = SPILS_SPBCG;

  /* Set SPBCG parameters that were passed in call sequence */
  maxl1 = (maxl <= 0) ? IDA_SPILS_MAXL : maxl;
  idaspils_mem->s_maxl = maxl1;

  /* Set defaults for Jacobian-related fileds */
  idaspils_mem->s_jtimesDQ = TRUE;
  idaspils_mem->s_jtimes   = NULL;
  idaspils_mem->s_jdata    = NULL;

  /* Set defaults for preconditioner-related fields */
  idaspils_mem->s_pset   = NULL;
  idaspils_mem->s_psolve = NULL;
  idaspils_mem->s_pfree  = NULL;
  idaspils_mem->s_pdata  = IDA_mem->ida_user_data;

  /* Set default values for the rest of the Spbcg parameters */
  idaspils_mem->s_eplifac   = PT05;
  idaspils_mem->s_dqincfac  = ONE;

  idaspils_mem->s_last_flag = IDASPILS_SUCCESS;

  idaSpilsInitializeCounters(idaspils_mem);

  /* Set setupNonNull to FALSE */
  IDA_mem->ida_setupNonNull = FALSE;

  /* Allocate memory for ytemp, yptemp, and xx */

  idaspils_mem->s_ytemp = N_VClone(IDA_mem->ida_tempv1);
  if (idaspils_mem->s_ytemp == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPBCG", "IDASpbcg", MSGS_MEM_FAIL);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  idaspils_mem->s_yptemp = N_VClone(IDA_mem->ida_tempv1);
  if (idaspils_mem->s_yptemp == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPBCG", "IDASpbcg", MSGS_MEM_FAIL);
    N_VDestroy(idaspils_mem->s_ytemp);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  idaspils_mem->s_xx = N_VClone(IDA_mem->ida_tempv1);
  if (idaspils_mem->s_xx == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPBCG", "IDASpbcg", MSGS_MEM_FAIL);
    N_VDestroy(idaspils_mem->s_ytemp);
    N_VDestroy(idaspils_mem->s_yptemp);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, idaspils_mem->s_ytemp);
  idaspils_mem->s_sqrtN = SUNRsqrt(N_VDotProd(idaspils_mem->s_ytemp, idaspils_mem->s_ytemp));

  /* Call SpbcgMalloc to allocate workspace for Spbcg */
  spbcg_mem = NULL;
  spbcg_mem = SpbcgMalloc(maxl1, IDA_mem->ida_tempv1);
  if (spbcg_mem == NULL) {
    IDAProcessError(NULL, IDASPILS_MEM_FAIL, "IDASPBCG", "IDASpbcg", MSGS_MEM_FAIL);
    N_VDestroy(idaspils_mem->s_ytemp);
    N_VDestroy(idaspils_mem->s_yptemp);
    N_VDestroy(idaspils_mem->s_xx);
    free(idaspils_mem); idaspils_mem = NULL;
    return(IDASPILS_MEM_FAIL);
  }

  /* Attach SPBCG memory to spils memory structure */
  idaspils_mem->s_spils_mem = (void *)spbcg_mem;

  /* Attach linear solver memory to the integrator memory */
  IDA_mem->ida_lmem = idaspils_mem;

  return(IDASPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDASPBCG interface routines
 * -----------------------------------------------------------------
 */

static int IDASpbcgInit(IDAMem IDA_mem)
{
  IDASpilsMem idaspils_mem;
  SpbcgMem spbcg_mem;

  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;
  spbcg_mem = (SpbcgMem) idaspils_mem->s_spils_mem;

  idaSpilsInitializeCounters(idaspils_mem);

  /* Set setupNonNull to TRUE iff there is preconditioning with setup */
  IDA_mem->ida_setupNonNull = (idaspils_mem->s_psolve != NULL) && (idaspils_mem->s_pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (idaspils_mem->s_jtimesDQ) {
    idaspils_mem->s_jtimes = IDASpilsDQJtimes;
    idaspils_mem->s_jdata = IDA_mem;
  } else {
    idaspils_mem->s_jdata = IDA_mem->ida_user_data;
  }

  /*  Set maxl in the SPBCG memory in case it was changed by the user */
  spbcg_mem->l_max  = idaspils_mem->s_maxl;

  idaspils_mem->s_last_flag = IDASPILS_SUCCESS;

  return(0);
}

static int IDASpbcgSetup(IDAMem IDA_mem, 
                         N_Vector yy_p, N_Vector yp_p, N_Vector rr_p, 
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int retval;
  IDASpilsMem idaspils_mem;

  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  /* Call user setup routine pset and update counter npe */
  retval = idaspils_mem->s_pset(IDA_mem->ida_tn, yy_p, yp_p, rr_p, IDA_mem->ida_cj, idaspils_mem->s_pdata,
                tmp1, tmp2, tmp3);
  idaspils_mem->s_npe++;

  if (retval < 0) {
    IDAProcessError(IDA_mem, SPBCG_PSET_FAIL_UNREC, "IDASPBCG", "IDASpbcgSetup", MSGS_PSET_FAILED);
    idaspils_mem->s_last_flag = SPBCG_PSET_FAIL_UNREC;
    return(-1);
  }
  if (retval > 0) {
    idaspils_mem->s_last_flag = SPBCG_PSET_FAIL_REC;
    return(+1);
  }

  idaspils_mem->s_last_flag = SPBCG_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcgSolve
 * -----------------------------------------------------------------
 * Note: The x-scaling and b-scaling arrays are both equal to weight.
 *
 * We set the initial guess, x = 0, then call SpbcgSolve.
 * We copy the solution x into b, and update the counters nli, nps,
 * and ncfl. If SpbcgSolve returned nli_inc = 0 (hence x = 0), we
 * take the SPBCG vtemp vector (= P_inverse F) as the correction
 * vector instead. Finally, we set the return value according to the
 * success of SpbcgSolve.
 * -----------------------------------------------------------------
 */

static int IDASpbcgSolve(IDAMem IDA_mem, N_Vector bb, N_Vector weight,
                         N_Vector yy_now, N_Vector yp_now, N_Vector rr_now)
{
  IDASpilsMem idaspils_mem;
  SpbcgMem spbcg_mem;
  int pretype, nli_inc, nps_inc, retval;
  realtype res_norm;

  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  spbcg_mem = (SpbcgMem) idaspils_mem->s_spils_mem;

  /* Set SpbcgSolve convergence test constant epslin, in terms of the
     Newton convergence test constant epsNewt and safety factors. The factor
     sqrt(Neq) assures that the Bi-CGSTAB convergence test is applied to the
     WRMS norm of the residual vector, rather than the weighted L2 norm. */
  idaspils_mem->s_epslin = idaspils_mem->s_sqrtN * idaspils_mem->s_eplifac * IDA_mem->ida_epsNewt;

  /* Set vectors ycur, ypcur, and rcur for use by the Atimes and Psolve */
  idaspils_mem->s_ycur = yy_now;
  idaspils_mem->s_ypcur = yp_now;
  idaspils_mem->s_rcur = rr_now;

  /* Set SpbcgSolve inputs pretype and initial guess xx = 0 */  
  pretype = (idaspils_mem->s_psolve == NULL) ? PREC_NONE : PREC_LEFT;
  N_VConst(ZERO, idaspils_mem->s_xx);
  
  /* Call SpbcgSolve and copy xx to bb */
  retval = SpbcgSolve(spbcg_mem, IDA_mem, idaspils_mem->s_xx, bb, pretype, idaspils_mem->s_epslin,
                      IDA_mem, weight, weight, IDASpilsAtimes,
                      IDASpilsPSolve, &res_norm, &nli_inc, &nps_inc);

  if (nli_inc == 0) N_VScale(ONE, SPBCG_VTEMP(spbcg_mem), bb);
  else N_VScale(ONE, idaspils_mem->s_xx, bb);
  
  /* Increment counters nli, nps, and return if successful */
  idaspils_mem->s_nli += nli_inc;
  idaspils_mem->s_nps += nps_inc;
  if (retval != SPBCG_SUCCESS) idaspils_mem->s_ncfl++;

  /* Interpret return value from SpbcgSolve */

  idaspils_mem->s_last_flag = retval;

  switch(retval) {

  case SPBCG_SUCCESS:
    return(0);
    break;
  case SPBCG_RES_REDUCED:
    return(1);
    break;
  case SPBCG_CONV_FAIL:
    return(1);
    break;
  case SPBCG_PSOLVE_FAIL_REC:
    return(1);
    break;
  case SPBCG_ATIMES_FAIL_REC:
    return(1);
    break;
  case SPBCG_MEM_NULL:
    return(-1);
    break;
  case SPBCG_ATIMES_FAIL_UNREC:
    IDAProcessError(IDA_mem, SPBCG_ATIMES_FAIL_UNREC, "IDaSPBCG", "IDASpbcgSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SPBCG_PSOLVE_FAIL_UNREC:
    IDAProcessError(IDA_mem, SPBCG_PSOLVE_FAIL_UNREC, "IDASPBCG", "IDASpbcgSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : IDASpbcgPerf
 * -----------------------------------------------------------------
 * This routine handles performance monitoring specific to the
 * IDASPBCG linear solver. When perftask = 0, it saves values of
 * various counters. When perftask = 1, it examines difference
 * quotients in these counters, and depending on their values, it
 * prints up to three warning messages. Messages are printed up to
 * a maximum of 10 times.
 * -----------------------------------------------------------------
 */

static int IDASpbcgPerf(IDAMem IDA_mem, int perftask)
{
  IDASpilsMem idaspils_mem;
  realtype avdim, rcfn, rcfl;
  long int nstd, nnid;
  booleantype lavd, lcfn, lcfl;

  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  if (perftask == 0) {
    idaspils_mem->s_nst0 = IDA_mem->ida_nst;  idaspils_mem->s_nni0 = IDA_mem->ida_nni;  idaspils_mem->s_nli0 = idaspils_mem->s_nli;
    idaspils_mem->s_ncfn0 = IDA_mem->ida_ncfn;  idaspils_mem->s_ncfl0 = idaspils_mem->s_ncfl;  
    idaspils_mem->s_nwarn = 0;
    return(0);
  }

  nstd = IDA_mem->ida_nst - idaspils_mem->s_nst0;  nnid = IDA_mem->ida_nni - idaspils_mem->s_nni0;
  if (nstd == 0 || nnid == 0) return(0);
  avdim = (realtype) ((idaspils_mem->s_nli - idaspils_mem->s_nli0)/((realtype) nnid));
  rcfn = (realtype) ((IDA_mem->ida_ncfn - idaspils_mem->s_ncfn0)/((realtype) nstd));
  rcfl = (realtype) ((idaspils_mem->s_ncfl - idaspils_mem->s_ncfl0)/((realtype) nnid));
  lavd = (avdim > ((realtype) idaspils_mem->s_maxl));
  lcfn = (rcfn > PT9);
  lcfl = (rcfl > PT9);
  if (!(lavd || lcfn || lcfl)) return(0);
  idaspils_mem->s_nwarn++;
  if (idaspils_mem->s_nwarn > 10) return(1);
  if (lavd) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPBCG", "IDASpbcgPerf", MSGS_AVD_WARN, IDA_mem->ida_tn, avdim);
  if (lcfn) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPBCG", "IDASpbcgPerf", MSGS_CFN_WARN, IDA_mem->ida_tn, rcfn);
  if (lcfl) 
    IDAProcessError(IDA_mem, IDA_WARNING, "IDASPBCG", "IDASpbcgPerf", MSGS_CFL_WARN, IDA_mem->ida_tn, rcfl);

  return(0);
}

static int IDASpbcgFree(IDAMem IDA_mem)
{
  IDASpilsMem idaspils_mem;
  SpbcgMem spbcg_mem;

  idaspils_mem = (IDASpilsMem) IDA_mem->ida_lmem;

  N_VDestroy(idaspils_mem->s_ytemp);
  N_VDestroy(idaspils_mem->s_yptemp);
  N_VDestroy(idaspils_mem->s_xx);

  spbcg_mem = (SpbcgMem) idaspils_mem->s_spils_mem;
  SpbcgFree(spbcg_mem);

  if (idaspils_mem->s_pfree != NULL) (idaspils_mem->s_pfree)(IDA_mem);

  free(idaspils_mem); idaspils_mem = NULL;

  return(0);
}
