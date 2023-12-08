/* -----------------------------------------------------------------
 * Programmer(s): Alan Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for the optional inputs and
 * outputs for the IDA solver.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "ida_ls_impl.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_types.h"

#define ZERO   SUN_RCONST(0.0)
#define HALF   SUN_RCONST(0.5)
#define ONE    SUN_RCONST(1.0)
#define TWOPT5 SUN_RCONST(2.5)

/*
 * =================================================================
 * IDA optional input functions
 * =================================================================
 */

int IDASetDeltaCjLSetup(void* ida_mem, sunrealtype dcj)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetDeltaCjLSetup", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (dcj < ZERO || dcj >= ONE) { IDA_mem->ida_dcj = DCJ_DEFAULT; }
  else { IDA_mem->ida_dcj = dcj; }

  return (IDA_SUCCESS);
}

int IDASetErrHandlerFn(void* ida_mem, IDAErrHandlerFn ehfun, void* eh_data)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetErrHandlerFn", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_ehfun   = ehfun;
  IDA_mem->ida_eh_data = eh_data;

  return (IDA_SUCCESS);
}

int IDASetErrFile(void* ida_mem, FILE* errfp)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetErrFile", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_errfp = errfp;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetUserData(void* ida_mem, void* user_data)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetUserData", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_user_data = user_data;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetEtaFixedStepBounds(void* ida_mem, sunrealtype eta_min_fx,
                             sunrealtype eta_max_fx)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetEtaFixedStepBounds",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* set allowed value or use default */
  if (eta_min_fx >= ZERO && eta_min_fx <= ONE)
  {
    IDA_mem->ida_eta_min_fx = eta_min_fx;
  }
  else { IDA_mem->ida_eta_min_fx = ETA_MIN_FX_DEFAULT; }

  if (eta_max_fx >= ONE) { IDA_mem->ida_eta_max_fx = eta_max_fx; }
  else { IDA_mem->ida_eta_max_fx = ETA_MAX_FX_DEFAULT; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetEtaMax(void* ida_mem, sunrealtype eta_max)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetEtaMax", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* set allowed value or use default */
  if (eta_max <= ONE) { IDA_mem->ida_eta_max = ETA_MAX_DEFAULT; }
  else { IDA_mem->ida_eta_max = eta_max; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetEtaMin(void* ida_mem, sunrealtype eta_min)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetEtaMin", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* set allowed value or use default */
  if (eta_min <= ZERO || eta_min >= ONE)
  {
    IDA_mem->ida_eta_min = ETA_MIN_DEFAULT;
  }
  else { IDA_mem->ida_eta_min = eta_min; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetEtaLow(void* ida_mem, sunrealtype eta_low)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetEtaLow", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* set allowed value or use default */
  if (eta_low <= ZERO || eta_low >= ONE)
  {
    IDA_mem->ida_eta_low = ETA_LOW_DEFAULT;
  }
  else { IDA_mem->ida_eta_low = eta_low; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetEtaMinErrFail(void* ida_mem, sunrealtype eta_min_ef)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetEtaMinErrFail", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* set allowed value or use default */
  if (eta_min_ef <= ZERO || eta_min_ef >= ONE)
  {
    IDA_mem->ida_eta_min_ef = ETA_MIN_EF_DEFAULT;
  }
  else { IDA_mem->ida_eta_min_ef = eta_min_ef; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetEtaConvFail(void* ida_mem, sunrealtype eta_cf)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetEtaConvFail", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  /* set allowed value or use default */
  if (eta_cf <= ZERO || eta_cf >= ONE) { IDA_mem->ida_eta_cf = ETA_CF_DEFAULT; }
  else { IDA_mem->ida_eta_cf = eta_cf; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxOrd(void* ida_mem, int maxord)
{
  IDAMem IDA_mem;
  int maxord_alloc;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxOrd", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (maxord <= 0)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxOrd",
                    MSG_NEG_MAXORD);
    return (IDA_ILL_INPUT);
  }

  /* Cannot increase maximum order beyond the value that
     was used when allocating memory */
  maxord_alloc = IDA_mem->ida_maxord_alloc;

  if (maxord > maxord_alloc)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxOrd",
                    MSG_BAD_MAXORD);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_maxord = SUNMIN(maxord, MAXORD_DEFAULT);

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumSteps(void* ida_mem, long int mxsteps)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumSteps", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */

  if (mxsteps == 0) { IDA_mem->ida_mxstep = MXSTEP_DEFAULT; }
  else { IDA_mem->ida_mxstep = mxsteps; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetInitStep(void* ida_mem, sunrealtype hin)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetInitStep", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_hin = hin;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxStep(void* ida_mem, sunrealtype hmax)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxStep", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (hmax < ZERO)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxStep", MSG_NEG_HMAX);
    return (IDA_ILL_INPUT);
  }

  /* Passing 0 sets hmax = infinity */
  if (hmax == ZERO)
  {
    IDA_mem->ida_hmax_inv = HMAX_INV_DEFAULT;
    return (IDA_SUCCESS);
  }

  IDA_mem->ida_hmax_inv = ONE / hmax;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMinStep(void* ida_mem, sunrealtype hmin)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMinStep", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (hmin < ZERO)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMinStep", MSG_NEG_HMIN);
    return (IDA_ILL_INPUT);
  }

  /* Passing 0 sets hmin = zero */
  if (hmin == ZERO)
  {
    IDA_mem->ida_hmin = HMIN_DEFAULT;
    return (IDA_SUCCESS);
  }

  IDA_mem->ida_hmin = hmin;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetStopTime(void* ida_mem, sunrealtype tstop)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetStopTime", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  /* If IDASolve was called at least once, test if tstop is legal
   * (i.e. if it was not already passed).
   * If IDASetStopTime is called before the first call to IDASolve,
   * tstop will be checked in IDASolve. */
  if (IDA_mem->ida_nst > 0)
  {
    if ((tstop - IDA_mem->ida_tn) * IDA_mem->ida_hh < ZERO)
    {
      IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetStopTime",
                      MSG_BAD_TSTOP, tstop, IDA_mem->ida_tn);
      return (IDA_ILL_INPUT);
    }
  }

  IDA_mem->ida_tstop    = tstop;
  IDA_mem->ida_tstopset = SUNTRUE;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAClearStopTime(void* ida_mem)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAClearStopTime", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_tstopset = SUNFALSE;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetNonlinConvCoef(void* ida_mem, sunrealtype epcon)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetNonlinConvCoef",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (epcon <= ZERO)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetNonlinConvCoef",
                    MSG_NEG_EPCON);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_epcon = epcon;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxErrTestFails(void* ida_mem, int maxnef)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxErrTestFails",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_maxnef = maxnef;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxConvFails(void* ida_mem, int maxncf)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxConvFails", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_maxncf = maxncf;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNonlinIters(void* ida_mem, int maxcor)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNonlinIters",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  /* check that the NLS is non-NULL */
  if (IDA_mem->NLS == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_FAIL, "IDA", "IDASetMaxNonlinIters",
                    MSG_MEM_FAIL);
    return (IDA_MEM_FAIL);
  }

  return (SUNNonlinSolSetMaxIters(IDA_mem->NLS, maxcor));
}

/*-----------------------------------------------------------------*/

int IDASetSuppressAlg(void* ida_mem, sunbooleantype suppressalg)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetSuppressAlg", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_suppressalg = suppressalg;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetId(void* ida_mem, N_Vector id)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetId", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (id == NULL)
  {
    if (IDA_mem->ida_idMallocDone)
    {
      N_VDestroy(IDA_mem->ida_id);
      IDA_mem->ida_lrw -= IDA_mem->ida_lrw1;
      IDA_mem->ida_liw -= IDA_mem->ida_liw1;
    }
    IDA_mem->ida_idMallocDone = SUNFALSE;
    return (IDA_SUCCESS);
  }

  if (!(IDA_mem->ida_idMallocDone))
  {
    IDA_mem->ida_id = N_VClone(id);
    IDA_mem->ida_lrw += IDA_mem->ida_lrw1;
    IDA_mem->ida_liw += IDA_mem->ida_liw1;
    IDA_mem->ida_idMallocDone = SUNTRUE;
  }

  /* Load the id vector */

  N_VScale(ONE, id, IDA_mem->ida_id);

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetConstraints(void* ida_mem, N_Vector constraints)
{
  IDAMem IDA_mem;
  sunrealtype temptest;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetConstraints", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (constraints == NULL)
  {
    if (IDA_mem->ida_constraintsMallocDone)
    {
      N_VDestroy(IDA_mem->ida_constraints);
      IDA_mem->ida_lrw -= IDA_mem->ida_lrw1;
      IDA_mem->ida_liw -= IDA_mem->ida_liw1;
    }
    IDA_mem->ida_constraintsMallocDone = SUNFALSE;
    IDA_mem->ida_constraintsSet        = SUNFALSE;
    return (IDA_SUCCESS);
  }

  /* Test if required vector ops. are defined */

  if (constraints->ops->nvdiv == NULL || constraints->ops->nvmaxnorm == NULL ||
      constraints->ops->nvcompare == NULL ||
      constraints->ops->nvconstrmask == NULL ||
      constraints->ops->nvminquotient == NULL)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetConstraints",
                    MSG_BAD_NVECTOR);
    return (IDA_ILL_INPUT);
  }

  /*  Check the constraints vector */

  temptest = N_VMaxNorm(constraints);
  if ((temptest > TWOPT5) || (temptest < HALF))
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetConstraints",
                    MSG_BAD_CONSTR);
    return (IDA_ILL_INPUT);
  }

  if (!(IDA_mem->ida_constraintsMallocDone))
  {
    IDA_mem->ida_constraints = N_VClone(constraints);
    IDA_mem->ida_lrw += IDA_mem->ida_lrw1;
    IDA_mem->ida_liw += IDA_mem->ida_liw1;
    IDA_mem->ida_constraintsMallocDone = SUNTRUE;
  }

  /* Load the constraints vector */

  N_VScale(ONE, constraints, IDA_mem->ida_constraints);

  IDA_mem->ida_constraintsSet = SUNTRUE;

  return (IDA_SUCCESS);
}

/*
 * IDASetRootDirection
 *
 * Specifies the direction of zero-crossings to be monitored.
 * The default is to monitor both crossings.
 */

int IDASetRootDirection(void* ida_mem, int* rootdir)
{
  IDAMem IDA_mem;
  int i, nrt;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetRootDirection", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  nrt = IDA_mem->ida_nrtfn;
  if (nrt == 0)
  {
    IDAProcessError(NULL, IDA_ILL_INPUT, "IDA", "IDASetRootDirection",
                    MSG_NO_ROOT);
    return (IDA_ILL_INPUT);
  }

  for (i = 0; i < nrt; i++) { IDA_mem->ida_rootdir[i] = rootdir[i]; }

  return (IDA_SUCCESS);
}

/*
 * IDASetNoInactiveRootWarn
 *
 * Disables issuing a warning if some root function appears
 * to be identically zero at the beginning of the integration
 */

int IDASetNoInactiveRootWarn(void* ida_mem)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetNoInactiveRootWarn",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_mxgnull = 0;

  return (IDA_SUCCESS);
}

/*
 * =================================================================
 * IDA IC optional input functions
 * =================================================================
 */

int IDASetNonlinConvCoefIC(void* ida_mem, sunrealtype epiccon)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetNonlinConvCoefIC",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (epiccon <= ZERO)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetNonlinConvCoefIC",
                    MSG_BAD_EPICCON);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_epiccon = epiccon;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumStepsIC(void* ida_mem, int maxnh)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumStepsIC", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (maxnh <= 0)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxNumStepsIC",
                    MSG_BAD_MAXNH);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_maxnh = maxnh;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumJacsIC(void* ida_mem, int maxnj)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumJacsIC", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (maxnj <= 0)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxNumJacsIC",
                    MSG_BAD_MAXNJ);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_maxnj = maxnj;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxNumItersIC(void* ida_mem, int maxnit)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxNumItersIC", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (maxnit <= 0)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxNumItersIC",
                    MSG_BAD_MAXNIT);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_maxnit = maxnit;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetMaxBacksIC(void* ida_mem, int maxbacks)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetMaxBacksIC", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (maxbacks <= 0)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetMaxBacksIC",
                    MSG_IC_BAD_MAXBACKS);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_maxbacks = maxbacks;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetLineSearchOffIC(void* ida_mem, sunbooleantype lsoff)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetLineSearchOffIC",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  IDA_mem->ida_lsoff = lsoff;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDASetStepToleranceIC(void* ida_mem, sunrealtype steptol)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDASetStepToleranceIC",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (steptol <= ZERO)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDASetStepToleranceIC",
                    MSG_BAD_STEPTOL);
    return (IDA_ILL_INPUT);
  }

  IDA_mem->ida_steptol = steptol;

  return (IDA_SUCCESS);
}

/*
 * =================================================================
 * IDA optional input functions
 * =================================================================
 */

int IDAGetNumSteps(void* ida_mem, long int* nsteps)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumSteps", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nsteps = IDA_mem->ida_nst;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumResEvals(void* ida_mem, long int* nrevals)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumResEvals", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nrevals = IDA_mem->ida_nre;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumLinSolvSetups(void* ida_mem, long int* nlinsetups)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumLinSolvSetups",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nlinsetups = IDA_mem->ida_nsetups;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumErrTestFails(void* ida_mem, long int* netfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumErrTestFails",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *netfails = IDA_mem->ida_netf;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumBacktrackOps(void* ida_mem, long int* nbacktracks)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumBacktrackOps",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nbacktracks = IDA_mem->ida_nbacktr;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetConsistentIC(void* ida_mem, N_Vector yy0, N_Vector yp0)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetConsistentIC", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  if (IDA_mem->ida_kused != 0)
  {
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAGetConsistentIC",
                    MSG_TOO_LATE);
    return (IDA_ILL_INPUT);
  }

  if (yy0 != NULL) { N_VScale(ONE, IDA_mem->ida_phi[0], yy0); }
  if (yp0 != NULL) { N_VScale(ONE, IDA_mem->ida_phi[1], yp0); }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetLastOrder(void* ida_mem, int* klast)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetLastOrder", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *klast = IDA_mem->ida_kused;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentOrder(void* ida_mem, int* kcur)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentOrder", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *kcur = IDA_mem->ida_kk;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentCj(void* ida_mem, sunrealtype* cj)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentCjRatio",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *cj = IDA_mem->ida_cj;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentY(void* ida_mem, N_Vector* ycur)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentY", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *ycur = IDA_mem->ida_yy;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentYp(void* ida_mem, N_Vector* ypcur)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentYp", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *ypcur = IDA_mem->ida_yp;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetActualInitStep(void* ida_mem, sunrealtype* hinused)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetActualInitStep",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *hinused = IDA_mem->ida_h0u;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetLastStep(void* ida_mem, sunrealtype* hlast)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetLastStep", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *hlast = IDA_mem->ida_hused;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentStep(void* ida_mem, sunrealtype* hcur)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentStep", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *hcur = IDA_mem->ida_hh;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetCurrentTime(void* ida_mem, sunrealtype* tcur)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetCurrentTime", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *tcur = IDA_mem->ida_tn;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetTolScaleFactor(void* ida_mem, sunrealtype* tolsfact)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetTolScaleFactor",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *tolsfact = IDA_mem->ida_tolsf;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetErrWeights(void* ida_mem, N_Vector eweight)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetErrWeights", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  N_VScale(ONE, IDA_mem->ida_ewt, eweight);

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetEstLocalErrors(void* ida_mem, N_Vector ele)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetEstLocalErrors",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }
  IDA_mem = (IDAMem)ida_mem;

  N_VScale(ONE, IDA_mem->ida_ee, ele);

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetWorkSpace(void* ida_mem, long int* lenrw, long int* leniw)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetWorkSpace", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *leniw = IDA_mem->ida_liw;
  *lenrw = IDA_mem->ida_lrw;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetIntegratorStats(void* ida_mem, long int* nsteps, long int* nrevals,
                          long int* nlinsetups, long int* netfails, int* klast,
                          int* kcur, sunrealtype* hinused, sunrealtype* hlast,
                          sunrealtype* hcur, sunrealtype* tcur)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetIntegratorStats",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nsteps     = IDA_mem->ida_nst;
  *nrevals    = IDA_mem->ida_nre;
  *nlinsetups = IDA_mem->ida_nsetups;
  *netfails   = IDA_mem->ida_netf;
  *klast      = IDA_mem->ida_kused;
  *kcur       = IDA_mem->ida_kk;
  *hinused    = IDA_mem->ida_h0u;
  *hlast      = IDA_mem->ida_hused;
  *hcur       = IDA_mem->ida_hh;
  *tcur       = IDA_mem->ida_tn;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumGEvals(void* ida_mem, long int* ngevals)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumGEvals", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *ngevals = IDA_mem->ida_nge;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetRootInfo(void* ida_mem, int* rootsfound)
{
  IDAMem IDA_mem;
  int i, nrt;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetRootInfo", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  nrt = IDA_mem->ida_nrtfn;

  for (i = 0; i < nrt; i++) { rootsfound[i] = IDA_mem->ida_iroots[i]; }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumNonlinSolvIters(void* ida_mem, long int* nniters)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumNonlinSolvIters",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nniters = IDA_mem->ida_nni;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumNonlinSolvConvFails(void* ida_mem, long int* nnfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumNonlinSolvConvFails",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nnfails = IDA_mem->ida_nnf;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNonlinSolvStats(void* ida_mem, long int* nniters, long int* nnfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNonlinSolvStats",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nniters = IDA_mem->ida_nni;
  *nnfails = IDA_mem->ida_nnf;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetNumStepSolveFails(void* ida_mem, long int* nncfails)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetNumStepSolveFails",
                    MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *nncfails = IDA_mem->ida_ncfn;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAGetUserData(void* ida_mem, void** user_data)
{
  IDAMem IDA_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAGetUserData", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  *user_data = IDA_mem->ida_user_data;

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

int IDAPrintAllStats(void* ida_mem, FILE* outfile, SUNOutputFormat fmt)
{
  IDAMem IDA_mem;
  IDALsMem idals_mem;

  if (ida_mem == NULL)
  {
    IDAProcessError(NULL, IDA_MEM_NULL, "IDA", "IDAPrintAllStats", MSG_NO_MEM);
    return (IDA_MEM_NULL);
  }

  IDA_mem = (IDAMem)ida_mem;

  switch (fmt)
  {
  case SUN_OUTPUTFORMAT_TABLE:
    /* step and method stats */
    fprintf(outfile, "Current time                 = %" RSYM "\n",
            IDA_mem->ida_tn);
    fprintf(outfile, "Steps                        = %ld\n", IDA_mem->ida_nst);
    fprintf(outfile, "Error test fails             = %ld\n", IDA_mem->ida_netf);
    fprintf(outfile, "NLS step fails               = %ld\n", IDA_mem->ida_ncfn);
    fprintf(outfile, "Initial step size            = %" RSYM "\n",
            IDA_mem->ida_h0u);
    fprintf(outfile, "Last step size               = %" RSYM "\n",
            IDA_mem->ida_hused);
    fprintf(outfile, "Current step size            = %" RSYM "\n",
            IDA_mem->ida_hh);
    fprintf(outfile, "Last method order            = %d\n", IDA_mem->ida_kused);
    fprintf(outfile, "Current method order         = %d\n", IDA_mem->ida_kk);

    /* function evaluations */
    fprintf(outfile, "Residual fn evals            = %ld\n", IDA_mem->ida_nre);

    /* IC calculation stats */
    fprintf(outfile, "IC linesearch backtrack ops  = %d\n", IDA_mem->ida_nbacktr);

    /* nonlinear solver stats */
    fprintf(outfile, "NLS iters                    = %ld\n", IDA_mem->ida_nni);
    fprintf(outfile, "NLS fails                    = %ld\n", IDA_mem->ida_nnf);
    if (IDA_mem->ida_nst > 0)
    {
      fprintf(outfile, "NLS iters per step           = %" RSYM "\n",
              (sunrealtype)IDA_mem->ida_nre / (sunrealtype)IDA_mem->ida_nst);
    }

    /* linear solver stats */
    fprintf(outfile, "LS setups                    = %ld\n",
            IDA_mem->ida_nsetups);
    if (IDA_mem->ida_lmem)
    {
      idals_mem = (IDALsMem)(IDA_mem->ida_lmem);
      fprintf(outfile, "Jac fn evals                 = %ld\n", idals_mem->nje);
      fprintf(outfile, "LS residual fn evals         = %ld\n", idals_mem->nreDQ);
      fprintf(outfile, "Prec setup evals             = %ld\n", idals_mem->npe);
      fprintf(outfile, "Prec solves                  = %ld\n", idals_mem->nps);
      fprintf(outfile, "LS iters                     = %ld\n", idals_mem->nli);
      fprintf(outfile, "LS fails                     = %ld\n", idals_mem->ncfl);
      fprintf(outfile, "Jac-times setups             = %ld\n",
              idals_mem->njtsetup);
      fprintf(outfile, "Jac-times evals              = %ld\n",
              idals_mem->njtimes);
      if (IDA_mem->ida_nni > 0)
      {
        fprintf(outfile, "LS iters per NLS iter        = %" RSYM "\n",
                (sunrealtype)idals_mem->nli / (sunrealtype)IDA_mem->ida_nni);
        fprintf(outfile, "Jac evals per NLS iter       = %" RSYM "\n",
                (sunrealtype)idals_mem->nje / (sunrealtype)IDA_mem->ida_nni);
        fprintf(outfile, "Prec evals per NLS iter      = %" RSYM "\n",
                (sunrealtype)idals_mem->npe / (sunrealtype)IDA_mem->ida_nni);
      }
    }

    /* rootfinding stats */
    fprintf(outfile, "Root fn evals                = %ld\n", IDA_mem->ida_nge);
    break;

  case SUN_OUTPUTFORMAT_CSV:
    /* step and method stats */
    fprintf(outfile, "Time,%" RSYM, IDA_mem->ida_tn);
    fprintf(outfile, ",Steps,%ld", IDA_mem->ida_nst);
    fprintf(outfile, ",Error test fails,%ld", IDA_mem->ida_netf);
    fprintf(outfile, ",NLS step fails,%ld", IDA_mem->ida_ncfn);
    fprintf(outfile, ",Initial step size,%" RSYM, IDA_mem->ida_h0u);
    fprintf(outfile, ",Last step size,%" RSYM, IDA_mem->ida_hused);
    fprintf(outfile, ",Current step size,%" RSYM, IDA_mem->ida_hh);
    fprintf(outfile, ",Last method order,%d", IDA_mem->ida_kused);
    fprintf(outfile, ",Current method order,%d", IDA_mem->ida_kk);

    /* function evaluations */
    fprintf(outfile, ",Residual fn evals,%ld", IDA_mem->ida_nre);

    /* IC calculation stats */
    fprintf(outfile, ",IC linesearch backtrack ops,%d", IDA_mem->ida_nbacktr);

    /* nonlinear solver stats */
    fprintf(outfile, ",NLS iters,%ld", IDA_mem->ida_nni);
    fprintf(outfile, ",NLS fails,%ld", IDA_mem->ida_nnf);
    if (IDA_mem->ida_nst > 0)
    {
      fprintf(outfile, ",NLS iters per step,%" RSYM,
              (sunrealtype)IDA_mem->ida_nre / (sunrealtype)IDA_mem->ida_nst);
    }
    else { fprintf(outfile, ",NLS iters per step,0"); }

    /* linear solver stats */
    fprintf(outfile, ",LS setups,%ld", IDA_mem->ida_nsetups);
    if (IDA_mem->ida_lmem)
    {
      idals_mem = (IDALsMem)(IDA_mem->ida_lmem);
      fprintf(outfile, ",Jac fn evals,%ld", idals_mem->nje);
      fprintf(outfile, ",LS residual evals,%ld", idals_mem->nreDQ);
      fprintf(outfile, ",Prec setup evals,%ld", idals_mem->npe);
      fprintf(outfile, ",Prec solves,%ld", idals_mem->nps);
      fprintf(outfile, ",LS iters,%ld", idals_mem->nli);
      fprintf(outfile, ",LS fails,%ld", idals_mem->ncfl);
      fprintf(outfile, ",Jac-times setups,%ld", idals_mem->njtsetup);
      fprintf(outfile, ",Jac-times evals,%ld", idals_mem->njtimes);
      if (IDA_mem->ida_nni > 0)
      {
        fprintf(outfile, ",LS iters per NLS iter,%" RSYM,
                (sunrealtype)idals_mem->nli / (sunrealtype)IDA_mem->ida_nni);
        fprintf(outfile, ",Jac evals per NLS iter,%" RSYM,
                (sunrealtype)idals_mem->nje / (sunrealtype)IDA_mem->ida_nni);
        fprintf(outfile, ",Prec evals per NLS iter,%" RSYM,
                (sunrealtype)idals_mem->npe / (sunrealtype)IDA_mem->ida_nni);
      }
      else
      {
        fprintf(outfile, ",LS iters per NLS iter,0");
        fprintf(outfile, ",Jac evals per NLS iter,0");
        fprintf(outfile, ",Prec evals per NLS iter,0");
      }
    }

    /* rootfinding stats */
    fprintf(outfile, ",Root fn evals,%ld", IDA_mem->ida_nge);
    fprintf(outfile, "\n");
    break;

  default:
    IDAProcessError(IDA_mem, IDA_ILL_INPUT, "IDA", "IDAPrintAllStats",
                    "Invalid formatting option.");
    return (IDA_ILL_INPUT);
  }

  return (IDA_SUCCESS);
}

/*-----------------------------------------------------------------*/

char* IDAGetReturnFlagName(long int flag)
{
  char* name;

  name = (char*)malloc(24 * sizeof(char));

  switch (flag)
  {
  case IDA_SUCCESS: sprintf(name, "IDA_SUCCESS"); break;
  case IDA_TSTOP_RETURN: sprintf(name, "IDA_TSTOP_RETURN"); break;
  case IDA_ROOT_RETURN: sprintf(name, "IDA_ROOT_RETURN"); break;
  case IDA_TOO_MUCH_WORK: sprintf(name, "IDA_TOO_MUCH_WORK"); break;
  case IDA_TOO_MUCH_ACC: sprintf(name, "IDA_TOO_MUCH_ACC"); break;
  case IDA_ERR_FAIL: sprintf(name, "IDA_ERR_FAIL"); break;
  case IDA_CONV_FAIL: sprintf(name, "IDA_CONV_FAIL"); break;
  case IDA_LINIT_FAIL: sprintf(name, "IDA_LINIT_FAIL"); break;
  case IDA_LSETUP_FAIL: sprintf(name, "IDA_LSETUP_FAIL"); break;
  case IDA_LSOLVE_FAIL: sprintf(name, "IDA_LSOLVE_FAIL"); break;
  case IDA_CONSTR_FAIL: sprintf(name, "IDA_CONSTR_FAIL"); break;
  case IDA_RES_FAIL: sprintf(name, "IDA_RES_FAIL"); break;
  case IDA_FIRST_RES_FAIL: sprintf(name, "IDA_FIRST_RES_FAIL"); break;
  case IDA_REP_RES_ERR: sprintf(name, "IDA_REP_RES_ERR"); break;
  case IDA_RTFUNC_FAIL: sprintf(name, "IDA_RTFUNC_FAIL"); break;
  case IDA_MEM_FAIL: sprintf(name, "IDA_MEM_FAIL"); break;
  case IDA_MEM_NULL: sprintf(name, "IDA_MEM_NULL"); break;
  case IDA_ILL_INPUT: sprintf(name, "IDA_ILL_INPUT"); break;
  case IDA_NO_MALLOC: sprintf(name, "IDA_NO_MALLOC"); break;
  case IDA_BAD_T: sprintf(name, "IDA_BAD_T"); break;
  case IDA_BAD_EWT: sprintf(name, "IDA_BAD_EWT"); break;
  case IDA_NO_RECOVERY: sprintf(name, "IDA_NO_RECOVERY"); break;
  case IDA_LINESEARCH_FAIL: sprintf(name, "IDA_LINESEARCH_FAIL"); break;
  case IDA_NLS_SETUP_FAIL: sprintf(name, "IDA_NLS_SETUP_FAIL"); break;
  case IDA_NLS_FAIL: sprintf(name, "IDA_NLS_FAIL"); break;
  default: sprintf(name, "NONE");
  }

  return (name);
}
