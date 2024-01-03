/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKode MRIStep module.
 * -----------------------------------------------------------------*/

#ifndef _MRISTEP_H
#define _MRISTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_ls.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * MRIStep Constants
 * ----------------- */

/* MRIStep method types */
typedef enum
{
  MRISTEP_EXPLICIT,
  MRISTEP_IMPLICIT,
  MRISTEP_IMEX
} MRISTEP_METHOD_TYPE;

typedef enum
{
  ARKODE_MRI_NONE    = -1, /* ensure enum is signed int */
  ARKODE_MIN_MRI_NUM = 200,
  ARKODE_MIS_KW3     = ARKODE_MIN_MRI_NUM,
  ARKODE_MRI_GARK_ERK33a,
  ARKODE_MRI_GARK_ERK45a,
  ARKODE_MRI_GARK_IRK21a,
  ARKODE_MRI_GARK_ESDIRK34a,
  ARKODE_MRI_GARK_ESDIRK46a,
  ARKODE_IMEX_MRI_GARK3a,
  ARKODE_IMEX_MRI_GARK3b,
  ARKODE_IMEX_MRI_GARK4,
  ARKODE_MAX_MRI_NUM = ARKODE_IMEX_MRI_GARK4
} ARKODE_MRITableID;

/* Default MRI coupling tables for each order */

static const int MRISTEP_DEFAULT_3         = ARKODE_MIS_KW3;
static const int MRISTEP_DEFAULT_EXPL_3    = ARKODE_MIS_KW3;
static const int MRISTEP_DEFAULT_EXPL_4    = ARKODE_MRI_GARK_ERK45a;
static const int MRISTEP_DEFAULT_IMPL_SD_2 = ARKODE_MRI_GARK_IRK21a;
static const int MRISTEP_DEFAULT_IMPL_SD_3 = ARKODE_MRI_GARK_ESDIRK34a;
static const int MRISTEP_DEFAULT_IMPL_SD_4 = ARKODE_MRI_GARK_ESDIRK46a;
static const int MRISTEP_DEFAULT_IMEX_SD_3 = ARKODE_IMEX_MRI_GARK3b;
static const int MRISTEP_DEFAULT_IMEX_SD_4 = ARKODE_IMEX_MRI_GARK4;

/* ------------------------------------
 * MRIStep Inner Stepper Function Types
 * ------------------------------------ */

typedef int (*MRIStepInnerEvolveFn)(MRIStepInnerStepper stepper, sunrealtype t0,
                                    sunrealtype tout, N_Vector y);

typedef int (*MRIStepInnerFullRhsFn)(MRIStepInnerStepper stepper, sunrealtype t,
                                     N_Vector y, N_Vector f, int mode);

typedef int (*MRIStepInnerResetFn)(MRIStepInnerStepper stepper, sunrealtype tR,
                                   N_Vector yR);

/*---------------------------------------------------------------
  MRI coupling data structure and associated utility routines
  ---------------------------------------------------------------*/
struct MRIStepCouplingMem
{
  int nmat;         /* number of MRI coupling matrices                   */
  int stages;       /* size of coupling matrices (stages * stages)       */
  int q;            /* method order of accuracy                          */
  int p;            /* embedding order of accuracy                       */
  sunrealtype* c;   /* stage abscissae                                   */
  sunrealtype*** W; /* explicit coupling matrices [nmat][stages][stages] */
  sunrealtype*** G; /* implicit coupling matrices [nmat][stages][stages] */
};

typedef _SUNDIALS_STRUCT_ MRIStepCouplingMem* MRIStepCoupling;

/* Accessor routine to load built-in MRI table */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID method);

/* Accessor routine to load built-in MRI table from string */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_LoadTableByName(const char* method);

/* Utility routines to allocate/free/output coupling table structures */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages,
                                                      MRISTEP_METHOD_TYPE type);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages,
                                                       int q, int p,
                                                       sunrealtype* W,
                                                       sunrealtype* G,
                                                       sunrealtype* c);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B,
                                                         int q, int p);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling MRIC);
SUNDIALS_EXPORT void MRIStepCoupling_Space(MRIStepCoupling MRIC,
                                           sunindextype* liw, sunindextype* lrw);
SUNDIALS_EXPORT void MRIStepCoupling_Free(MRIStepCoupling MRIC);
SUNDIALS_EXPORT void MRIStepCoupling_Write(MRIStepCoupling MRIC, FILE* outfile);

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*MRIStepPreInnerFn)(sunrealtype t, N_Vector* f, int nvecs,
                                 void* user_data);

typedef int (*MRIStepPostInnerFn)(sunrealtype t, N_Vector y, void* user_data);

/* -------------------
 * Exported Functions
 * ------------------- */

/* Create, Resize, and Reinitialization functions */
SUNDIALS_EXPORT void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0,
                                    N_Vector y0, MRIStepInnerStepper stepper,
                                    SUNContext sunctx);

SUNDIALS_EXPORT int MRIStepResize(void* arkode_mem, N_Vector ynew, sunrealtype t0,
                                  ARKVecResizeFn resize, void* resize_data);

SUNDIALS_EXPORT int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi,
                                  sunrealtype t0, N_Vector y0);

SUNDIALS_EXPORT int MRIStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);

/* Tolerance input functions */
SUNDIALS_EXPORT int MRIStepSStolerances(void* arkode_mem, sunrealtype reltol,
                                        sunrealtype abstol);
SUNDIALS_EXPORT int MRIStepSVtolerances(void* arkode_mem, sunrealtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int MRIStepWFtolerances(void* arkode_mem, ARKEwtFn efun);

/* Linear solver set function */
SUNDIALS_EXPORT int MRIStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS,
                                           SUNMatrix A);

/* Rootfinding initialization */
SUNDIALS_EXPORT int MRIStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);

/* Optional input functions -- must be called AFTER MRIStepCreate */
SUNDIALS_EXPORT int MRIStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int MRIStepSetOrder(void* arkode_mem, int ord);
SUNDIALS_EXPORT int MRIStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_EXPORT int MRIStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_EXPORT int MRIStepSetDenseOrder(void* arkode_mem, int dord);
SUNDIALS_EXPORT int MRIStepSetNonlinearSolver(void* arkode_mem,
                                              SUNNonlinearSolver NLS);
SUNDIALS_EXPORT int MRIStepSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fs);
SUNDIALS_EXPORT int MRIStepSetLinear(void* arkode_mem, int timedepend);
SUNDIALS_EXPORT int MRIStepSetNonlinear(void* arkode_mem);
SUNDIALS_EXPORT int MRIStepSetCoupling(void* arkode_mem, MRIStepCoupling MRIC);
SUNDIALS_EXPORT int MRIStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_EXPORT int MRIStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown);
SUNDIALS_EXPORT int MRIStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv);
SUNDIALS_EXPORT int MRIStepSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax);
SUNDIALS_EXPORT int MRIStepSetLSetupFrequency(void* arkode_mem, int msbp);
SUNDIALS_EXPORT int MRIStepSetPredictorMethod(void* arkode_mem, int method);
SUNDIALS_EXPORT int MRIStepSetMaxNonlinIters(void* arkode_mem, int maxcor);
SUNDIALS_EXPORT int MRIStepSetNonlinConvCoef(void* arkode_mem,
                                             sunrealtype nlscoef);
SUNDIALS_EXPORT int MRIStepSetMaxHnilWarns(void* arkode_mem, int mxhnil);
SUNDIALS_EXPORT int MRIStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_EXPORT int MRIStepSetInterpolateStopTime(void* arkode_mem,
                                                  sunbooleantype interp);
SUNDIALS_EXPORT int MRIStepClearStopTime(void* arkode_mem);
SUNDIALS_EXPORT int MRIStepSetFixedStep(void* arkode_mem, sunrealtype hsfixed);
SUNDIALS_EXPORT int MRIStepSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_EXPORT int MRIStepSetNoInactiveRootWarn(void* arkode_mem);
SUNDIALS_EXPORT int MRIStepSetUserData(void* arkode_mem, void* user_data);
SUNDIALS_EXPORT int MRIStepSetPostprocessStepFn(void* arkode_mem,
                                                ARKPostProcessFn ProcessStep);
SUNDIALS_EXPORT int MRIStepSetPostprocessStageFn(void* arkode_mem,
                                                 ARKPostProcessFn ProcessStage);
SUNDIALS_EXPORT int MRIStepSetPreInnerFn(void* arkode_mem,
                                         MRIStepPreInnerFn prefn);
SUNDIALS_EXPORT int MRIStepSetPostInnerFn(void* arkode_mem,
                                          MRIStepPostInnerFn postfn);
SUNDIALS_EXPORT int MRIStepSetStagePredictFn(void* arkode_mem,
                                             ARKStagePredictFn PredictStage);
SUNDIALS_EXPORT int MRIStepSetDeduceImplicitRhs(void* arkode_mem,
                                                sunbooleantype deduce);

/* Linear solver interface optional input functions -- must be called
   AFTER MRIStepSetLinearSolver */
SUNDIALS_EXPORT int MRIStepSetJacFn(void* arkode_mem, ARKLsJacFn jac);
SUNDIALS_EXPORT int MRIStepSetJacEvalFrequency(void* arkode_mem, long int msbj);
SUNDIALS_EXPORT int MRIStepSetLinearSolutionScaling(void* arkode_mem,
                                                    sunbooleantype onoff);
SUNDIALS_EXPORT int MRIStepSetEpsLin(void* arkode_mem, sunrealtype eplifac);
SUNDIALS_EXPORT int MRIStepSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac);
SUNDIALS_EXPORT int MRIStepSetPreconditioner(void* arkode_mem,
                                             ARKLsPrecSetupFn psetup,
                                             ARKLsPrecSolveFn psolve);
SUNDIALS_EXPORT int MRIStepSetJacTimes(void* arkode_mem,
                                       ARKLsJacTimesSetupFn jtsetup,
                                       ARKLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int MRIStepSetJacTimesRhsFn(void* arkode_mem,
                                            ARKRhsFn jtimesRhsFn);
SUNDIALS_EXPORT int MRIStepSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys);

/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int MRIStepEvolve(void* arkode_mem, sunrealtype tout,
                                  N_Vector yout, sunrealtype* tret, int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int MRIStepGetDky(void* arkode_mem, sunrealtype t, int k,
                                  N_Vector dky);

/* Utility function to update/compute y based on zcor */
SUNDIALS_EXPORT int MRIStepComputeState(void* arkode_mem, N_Vector zcor,
                                        N_Vector z);

/* Optional output functions */
SUNDIALS_EXPORT int MRIStepGetNumRhsEvals(void* arkode_mem, long int* nfse_evals,
                                          long int* nfsi_evals);
SUNDIALS_EXPORT int MRIStepGetNumLinSolvSetups(void* arkode_mem,
                                               long int* nlinsetups);
SUNDIALS_EXPORT int MRIStepGetCurrentCoupling(void* arkode_mem,
                                              MRIStepCoupling* MRIC);
SUNDIALS_EXPORT int MRIStepGetWorkSpace(void* arkode_mem, long int* lenrw,
                                        long int* leniw);
SUNDIALS_EXPORT int MRIStepGetNumSteps(void* arkode_mem, long int* nssteps);
SUNDIALS_EXPORT int MRIStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_EXPORT int MRIStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_EXPORT int MRIStepGetCurrentState(void* arkode_mem, N_Vector* state);
SUNDIALS_EXPORT int MRIStepGetCurrentGamma(void* arkode_mem, sunrealtype* gamma);
SUNDIALS_EXPORT int MRIStepGetTolScaleFactor(void* arkode_mem,
                                             sunrealtype* tolsfac);
SUNDIALS_EXPORT int MRIStepGetErrWeights(void* arkode_mem, N_Vector eweight);
SUNDIALS_EXPORT int MRIStepGetNumGEvals(void* arkode_mem, long int* ngevals);
SUNDIALS_EXPORT int MRIStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_EXPORT int MRIStepGetLastInnerStepFlag(void* arkode_mem, int* flag);
SUNDIALS_EXPORT int MRIStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_EXPORT int MRIStepPrintAllStats(void* arkode_mem, FILE* outfile,
                                         SUNOutputFormat fmt);
SUNDIALS_EXPORT char* MRIStepGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int MRIStepWriteParameters(void* arkode_mem, FILE* fp);

SUNDIALS_EXPORT int MRIStepWriteCoupling(void* arkode_mem, FILE* fp);

/* Nonlinear solver optional output functions */
SUNDIALS_EXPORT int MRIStepGetNonlinearSystemData(
  void* arkode_mem, sunrealtype* tcur, N_Vector* zpred, N_Vector* z,
  N_Vector* F, sunrealtype* gamma, N_Vector* sdata, void** user_data);
SUNDIALS_EXPORT int MRIStepGetNumNonlinSolvIters(void* arkode_mem,
                                                 long int* nniters);
SUNDIALS_EXPORT int MRIStepGetNumNonlinSolvConvFails(void* arkode_mem,
                                                     long int* nnfails);
SUNDIALS_EXPORT int MRIStepGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                                              long int* nnfails);
SUNDIALS_EXPORT int MRIStepGetNumStepSolveFails(void* arkode_mem,
                                                long int* nncfails);

/* Linear solver optional output functions */
SUNDIALS_EXPORT int MRIStepGetJac(void* arkode_mem, SUNMatrix* J);
SUNDIALS_EXPORT int MRIStepGetJacTime(void* arkode_mem, sunrealtype* t_J);
SUNDIALS_EXPORT int MRIStepGetJacNumSteps(void* arkode_mem, long* nst_J);
SUNDIALS_EXPORT int MRIStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS,
                                           long int* leniwLS);
SUNDIALS_EXPORT int MRIStepGetNumJacEvals(void* arkode_mem, long int* njevals);
SUNDIALS_EXPORT int MRIStepGetNumPrecEvals(void* arkode_mem, long int* npevals);
SUNDIALS_EXPORT int MRIStepGetNumPrecSolves(void* arkode_mem, long int* npsolves);
SUNDIALS_EXPORT int MRIStepGetNumLinIters(void* arkode_mem, long int* nliters);
SUNDIALS_EXPORT int MRIStepGetNumLinConvFails(void* arkode_mem,
                                              long int* nlcfails);
SUNDIALS_EXPORT int MRIStepGetNumJTSetupEvals(void* arkode_mem,
                                              long int* njtsetups);
SUNDIALS_EXPORT int MRIStepGetNumJtimesEvals(void* arkode_mem,
                                             long int* njvevals);
SUNDIALS_EXPORT int MRIStepGetNumLinRhsEvals(void* arkode_mem,
                                             long int* nfevalsLS);
SUNDIALS_EXPORT int MRIStepGetLastLinFlag(void* arkode_mem, long int* flag);

SUNDIALS_EXPORT char* MRIStepGetLinReturnFlagName(long int flag);

/* Free function */
SUNDIALS_EXPORT void MRIStepFree(void** arkode_mem);

/* Output the MRIStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void MRIStepPrintMem(void* arkode_mem, FILE* outfile);

/* Custom inner stepper functions */
SUNDIALS_EXPORT int MRIStepInnerStepper_Create(SUNContext sunctx,
                                               MRIStepInnerStepper* stepper);

SUNDIALS_EXPORT int MRIStepInnerStepper_Free(MRIStepInnerStepper* stepper);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetContent(MRIStepInnerStepper stepper,
                                                   void* content);

SUNDIALS_EXPORT int MRIStepInnerStepper_GetContent(MRIStepInnerStepper stepper,
                                                   void** content);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetEvolveFn(MRIStepInnerStepper stepper,
                                                    MRIStepInnerEvolveFn fn);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetFullRhsFn(MRIStepInnerStepper stepper,
                                                     MRIStepInnerFullRhsFn fn);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetResetFn(MRIStepInnerStepper stepper,
                                                   MRIStepInnerResetFn fn);

SUNDIALS_EXPORT int MRIStepInnerStepper_AddForcing(MRIStepInnerStepper stepper,
                                                   sunrealtype t, N_Vector f);

SUNDIALS_EXPORT int MRIStepInnerStepper_GetForcingData(
  MRIStepInnerStepper stepper, sunrealtype* tshift, sunrealtype* tscale,
  N_Vector** forcing, int* nforcing);

#ifdef __cplusplus
}
#endif

#endif
