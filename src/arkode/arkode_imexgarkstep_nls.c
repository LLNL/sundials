/* ----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------------------
 * Based on arkode_arkstep.c written by Daniel R. Reynolds @ SMU
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * This is the interface between IMEXGARKStep and the
 * SUNNonlinearSolver object
 * --------------------------------------------------------------------------*/

/*---------------------------------------------------------------
 imexgarkStep_Nls

 This routine attempts to solve the nonlinear system associated
 with a single implicit step of the linear multistep method.
 It calls one of imexgarkStep_NlsAccelFP, imexgarkStep_Ls or
 imexgarkStep_NlsNewton to do the work.

 Upon entry, the predicted solution is held in step_mem->zpred;
 this array is never changed throughout this routine.  If an
 initial attempt at solving the nonlinear system fails (e.g. due
 to a stale Jacobian), this allows for new attempts at the
 solution.

 Upon a successful solve, the solution is held in arkode_mem->ycur.
---------------------------------------------------------------*/
int imexgarkStep_Nls(ARKodeMem arkode_mem, int nflag)
{
  ARKodeIMEXGARKStepMem step_mem;

  /* access step memory structure */
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::IMEXGARKStep",
                    "imexgarkStep_Nls", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* call the appropriate solver */

  /*   fixed point */
  if (step_mem->use_fp)
    return(imexgarkStep_NlsAccelFP(arkode_mem, nflag));

  /*   linearly implicit (one Newton iteration) */
  if (step_mem->linear)
    return(imexgarkStep_Ls(arkode_mem, nflag));

  /*   Newton */
  return(imexgarkStep_NlsNewton(arkode_mem, nflag));

}


/*---------------------------------------------------------------
 imexgarkStep_NlsResid

 This routine evaluates the negative nonlinear residual for the
 additive Runge-Kutta method.  It assumes that any data from
 previous time steps/stages is contained in step_mem, and
 merely combines this old data with the current implicit ODE RHS
 vector, fy = fi(t,y), to compute the nonlinear residual r.

 At the ith stage, we compute the residual vector:
    r = -M*zi + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                     + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = -M*zp - M*zc + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                            + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = (-M*zc + gamma*Fi(zi)) + (M*yn - M*zp + data)
 where zi = zp + zc, and where
    zc is stored in the input z
    Fi(zi) is stored in fz
    (M*yn-M*zp+data) is stored in step_mem->sdata,
 so we really just compute
    r = -M*z + gamma*fz + step_mem->sdata

 Possible return values:  ARK_SUCCESS  or  ARK_RHSFUNC_FAIL
---------------------------------------------------------------*/
int imexgarkStep_NlsResid(ARKodeMem arkode_mem, N_Vector z, N_Vector fz,
                     N_Vector r)
{

  /* temporary variables */
  int retval;
  realtype c[3];
  N_Vector X[3];
  ARKodeIMEXGARKStepMem step_mem;

  /* access step memory structure */
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::IMEXGARKStep",
                    "imexgarkStep_NlsResid", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  /* put M*z in r */
  if (step_mem->mass_mem != NULL) {
    retval = step_mem->mmult((void *) arkode_mem, z, r);
    if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
    X[0] = r;
  } else {
    X[0] = z;
  }

  /* update with My, sdata and gamma*fy */
  c[0] = -ONE;
  c[1] = ONE;
  X[1] = step_mem->sdata;
  c[2] = step_mem->gamma;
  X[2] = fz;
  retval = N_VLinearCombination(3, c, X, r);
  if (retval != 0)  return(ARK_VECTOROP_ERR);
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
 imexgarkStep_Ls

 This routine attempts to solve the linear system associated
 with a single implicit step of the ARK method.  This should
 only be called if the user has specified that the implicit
 problem is linear.  In this routine, we assume that the problem
 depends linearly on the solution.  Additionally, if the Jacobian
 is not time dependent we only call lsetup on changes to gamma;
 otherwise we call lsetup at every call.  In all cases, we then
 call the user-specified linear solver (with a tight tolerance)
 to compute the time-evolved solution.

 Upon entry, the predicted solution is held in step_mem->zpred.

 Upon a successful solve, the solution is held in arkode_mem->ycur.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  -+
   ARK_LSETUP_FAIL    |-> halt the integration
   ARK_LSOLVE_FAIL   -+

   RHSFUNC_RECVR     --> predict again or stop if too many
---------------------------------------------------------------*/
int imexgarkStep_Ls(ARKodeMem arkode_mem, int nflag)
{
  N_Vector b, z, zcor;
  int convfail, retval, ier, is;
  booleantype callSetup;
  realtype del;
  ARKodeIMEXGARKStepMem step_mem;

  /* access step memory structure */
  if (arkode_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE::IMEXGARKStep",
                    "imexgarkStep_Ls", MSG_IMEXGARK_NO_STEP_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeIMEXGARKStepMem) arkode_mem->step_mem;

  b    = arkode_mem->tempv1;   /* rename tempv1 as b for readability */
  z    = arkode_mem->ycur;     /* rename ycur as z for readability */
  zcor = arkode_mem->tempv2;   /* rename tempv2 as zcor for readability */
  is   = step_mem->istage;

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = (nflag == FIRST_CALL) ? ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (step_mem->lsetup) {
    callSetup = (arkode_mem->firststage) ||
      (step_mem->linear_timedep) || (step_mem->msbp < 0) ||
      (SUNRabs(step_mem->gamrat-ONE) > step_mem->dgmax);
  } else {
    callSetup = SUNFALSE;
  }

  /* update implicit RHS, store in arkode_mem->Fi[is] */
  retval = step_mem->fi(arkode_mem->tcur, step_mem->zpred,
                           step_mem->Fi[is], arkode_mem->user_data);
  step_mem->nfi++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  /* update system matrix/factorization if necessary */
  if (callSetup) {

    /* Solver diagnostics reporting */
    if (arkode_mem->report)  fprintf(arkode_mem->diagfp, "  lsetup\n");

    /* call lsetup, using zcor, z and b as temporary vectors */
    ier = step_mem->lsetup(arkode_mem, convfail, arkode_mem->tcur,
                              step_mem->zpred, step_mem->Fi[is],
                              &step_mem->jcur, zcor, z, b);
    step_mem->nsetups++;
    callSetup = SUNFALSE;
    arkode_mem->firststage = SUNFALSE;
    step_mem->gamrat = step_mem->crate = ONE;
    step_mem->gammap = step_mem->gamma;
    step_mem->nstlp  = arkode_mem->nst;

    /* Return if lsetup failed */
    if (ier < 0) return(ARK_LSETUP_FAIL);
    if (ier > 0) return(CONV_FAIL);
  }

  /* Set zcor to zero and load prediction into y vector */
  N_VConst(ZERO, zcor);
  N_VScale(ONE, step_mem->zpred, z);


  /* Do a single Newton iteration */

  /*   Initialize temporary variables for use in iteration */
  step_mem->mnewt = 0;
  del = ZERO;

  /*   Set the stored residual norm to force an "accurate" initial linear solve */
  step_mem->eRNrm = RCONST(0.1) * step_mem->nlscoef;

  /*   Evaluate the nonlinear system residual, put result into b */
  retval = imexgarkStep_NlsResid(arkode_mem, zcor, step_mem->Fi[is], b);
  if (retval != ARK_SUCCESS)  return (ARK_RHSFUNC_FAIL);

  /*   Call the lsolve function */
  retval = step_mem->lsolve(arkode_mem, b, arkode_mem->tcur, z, step_mem->Fi[is],
                               step_mem->eRNrm, step_mem->mnewt);
  step_mem->nni++;
  if (retval != 0)  return (ARK_LSOLVE_FAIL);

  /*   Get WRMS norm of correction; add correction to zcor and z */
  del = N_VWrmsNorm(b, arkode_mem->ewt);
  N_VLinearSum(ONE, zcor, ONE, b, zcor);
  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, z);

  /*   Solver diagnostics reporting */
  if (arkode_mem->report)
    fprintf(arkode_mem->diagfp, "    newt  %i  %"RSYM"  %g\n", 0, del, 0.0);

  /* clean up and return */
  step_mem->jcur = SUNFALSE;
  return (ARK_SUCCESS);
}
