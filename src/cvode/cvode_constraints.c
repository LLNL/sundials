/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * Inequality constraint handling for CVODE
 * ---------------------------------------------------------------------------*/

/*
 * cvCheckConstraints
 *
 * This routine determines if the constraints of the problem
 * are satisfied by the proposed step
 *
 * Possible return values are:
 *
 *   CV_SUCCESS     ---> allows stepping forward
 *
 *   CONSTR_RECVR   ---> values failed to satisfy constraints
 *
 *   CV_CONSTR_FAIL ---> values failed to satisfy constraints with hmin
 */

static int cvCheckConstraints(CVodeMem cv_mem)
{
  sunbooleantype constraintsPassed;
  sunrealtype vnorm;
  N_Vector mm  = cv_mem->cv_ftemp;
  N_Vector tmp = cv_mem->cv_tempv;

  /* Get mask vector mm, set where constraints failed */
  constraintsPassed = N_VConstrMask(cv_mem->cv_constraints, cv_mem->cv_y, mm);
  if (constraintsPassed) { return (CV_SUCCESS); }

  printf("tn = %Lg\n", cv_mem->cv_tn);
  printf("hn = %Lg\n", cv_mem->cv_h);
  printf("constraints:\n");
  N_VPrint(cv_mem->cv_constraints);
  printf("y:\n");
  N_VPrint(cv_mem->cv_y);
  printf("mm:\n");
  N_VPrint(mm);

  /* Constraints not met */

  /* Compute correction to satisfy constraints */
#ifdef SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS
  if (cv_mem->cv_usefused)
  {
    cvCheckConstraints_fused(cv_mem->cv_constraints, cv_mem->cv_ewt,
                             cv_mem->cv_y, mm, tmp);
  }
  else
#endif
  {
    N_VCompare(ONEPT5, cv_mem->cv_constraints, tmp); /* a[i]=1 when |c[i]|=2  */
    N_VProd(tmp, cv_mem->cv_constraints, tmp);       /* a * c                 */
    N_VDiv(tmp, cv_mem->cv_ewt, tmp);                /* a * c * wt            */
    N_VLinearSum(ONE, cv_mem->cv_y, -PT1, tmp, tmp); /* y - 0.1 * a * c * wt  */
    N_VProd(tmp, mm, tmp);                           /* v = mm*(y-0.1*a*c*wt) */
  }

  vnorm = N_VWrmsNorm(tmp, cv_mem->cv_ewt); /* ||v|| */

  printf("correction:\n");
  N_VPrint(tmp);
  printf("correction norm = %Lg:\n", vnorm);
  printf("threshold = %Lg:\n", cv_mem->cv_tq[4]);

  /* If vector v of constraint corrections is small in norm, correct and
     accept this step */
  if (vnorm <= cv_mem->cv_tq[4])
  {
    printf("<<<<< Correct\n\n");
    N_VLinearSum(ONE, cv_mem->cv_acor, -ONE, tmp,
                 cv_mem->cv_acor); /* acor <- acor - v */
    return (CV_SUCCESS);
  }

  /* Return with error if |h| == hmin */
  if (SUNRabs(cv_mem->cv_h) <= cv_mem->cv_hmin * ONEPSM)
  {
    return (CV_CONSTR_FAIL);
  }

  /* Constraint correction is too large, reduce h by computing eta = h'/h */
  N_VLinearSum(ONE, cv_mem->cv_zn[0], -ONE, cv_mem->cv_y, tmp);
  N_VProd(mm, tmp, tmp);

  printf("tmp to compute eta:\n");
  N_VPrint(tmp);
  printf("z[0]:\n");
  N_VPrint(cv_mem->cv_zn[0]);

  cv_mem->cv_eta = PT9 * N_VMinQuotient(cv_mem->cv_zn[0], tmp);
  printf("eta1 = %Lg\n", cv_mem->cv_eta);
  cv_mem->cv_eta = SUNMAX(cv_mem->cv_eta, PT1);
  printf("eta2 = %Lg\n", cv_mem->cv_eta);
  cv_mem->cv_eta = SUNMAX(cv_mem->cv_eta,
                          cv_mem->cv_hmin / SUNRabs(cv_mem->cv_h));
  printf("eta3 = %Lg\n", cv_mem->cv_eta);

  printf(">>>>> Constraint limited step, eta = %Lg\n\n", cv_mem->cv_eta);

  /* Reattempt step with new step size */
  return (CONSTR_RECVR);
}
