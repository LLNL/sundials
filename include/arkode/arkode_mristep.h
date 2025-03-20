/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKODE MRIStep module.
 * -----------------------------------------------------------------*/

#ifndef _MRISTEP_H
#define _MRISTEP_H

#include <arkode/arkode.h>
#include <arkode/arkode_butcher_dirk.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_ls.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_stepper.h>

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
  MRISTEP_IMEX,
  MRISTEP_MERK,
  MRISTEP_SR
} MRISTEP_METHOD_TYPE;

/* MRI coupling table IDs */
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
  ARKODE_MRI_GARK_FORWARD_EULER,
  ARKODE_MRI_GARK_RALSTON2,
  ARKODE_MRI_GARK_ERK22a,
  ARKODE_MRI_GARK_ERK22b,
  ARKODE_MRI_GARK_RALSTON3,
  ARKODE_MRI_GARK_BACKWARD_EULER,
  ARKODE_MRI_GARK_IMPLICIT_MIDPOINT,
  ARKODE_IMEX_MRI_GARK_EULER,
  ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL,
  ARKODE_IMEX_MRI_GARK_MIDPOINT,
  ARKODE_MERK21,
  ARKODE_MERK32,
  ARKODE_MERK43,
  ARKODE_MERK54,
  ARKODE_IMEX_MRI_SR21,
  ARKODE_IMEX_MRI_SR32,
  ARKODE_IMEX_MRI_SR43,
  ARKODE_MAX_MRI_NUM = ARKODE_IMEX_MRI_SR43
} ARKODE_MRITableID;

/* Default MRI coupling tables for each order and type */
static const int MRISTEP_DEFAULT_EXPL_1 = ARKODE_MRI_GARK_FORWARD_EULER;
static const int MRISTEP_DEFAULT_EXPL_2 = ARKODE_MRI_GARK_ERK22b;
static const int MRISTEP_DEFAULT_EXPL_3 = ARKODE_MIS_KW3;
static const int MRISTEP_DEFAULT_EXPL_4 = ARKODE_MRI_GARK_ERK45a;

static const int MRISTEP_DEFAULT_EXPL_2_AD = ARKODE_MRI_GARK_ERK22b;
static const int MRISTEP_DEFAULT_EXPL_3_AD = ARKODE_MRI_GARK_ERK33a;
static const int MRISTEP_DEFAULT_EXPL_4_AD = ARKODE_MRI_GARK_ERK45a;
static const int MRISTEP_DEFAULT_EXPL_5_AD = ARKODE_MERK54;

static const int MRISTEP_DEFAULT_IMPL_SD_1 = ARKODE_MRI_GARK_BACKWARD_EULER;
static const int MRISTEP_DEFAULT_IMPL_SD_2 = ARKODE_MRI_GARK_IRK21a;
static const int MRISTEP_DEFAULT_IMPL_SD_3 = ARKODE_MRI_GARK_ESDIRK34a;
static const int MRISTEP_DEFAULT_IMPL_SD_4 = ARKODE_MRI_GARK_ESDIRK46a;

static const int MRISTEP_DEFAULT_IMEX_SD_1 = ARKODE_IMEX_MRI_GARK_EULER;
static const int MRISTEP_DEFAULT_IMEX_SD_2 = ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL;
static const int MRISTEP_DEFAULT_IMEX_SD_3 = ARKODE_IMEX_MRI_GARK3b;
static const int MRISTEP_DEFAULT_IMEX_SD_4 = ARKODE_IMEX_MRI_GARK4;

static const int MRISTEP_DEFAULT_IMEX_SD_2_AD = ARKODE_IMEX_MRI_SR21;
static const int MRISTEP_DEFAULT_IMEX_SD_3_AD = ARKODE_IMEX_MRI_SR32;
static const int MRISTEP_DEFAULT_IMEX_SD_4_AD = ARKODE_IMEX_MRI_SR43;

/* ------------------------------------
 * MRIStep Inner Stepper Function Types
 * ------------------------------------ */

typedef int (*MRIStepInnerEvolveFn)(MRIStepInnerStepper stepper, sunrealtype t0,
                                    sunrealtype tout, N_Vector y);

typedef int (*MRIStepInnerFullRhsFn)(MRIStepInnerStepper stepper, sunrealtype t,
                                     N_Vector y, N_Vector f, int mode);

typedef int (*MRIStepInnerResetFn)(MRIStepInnerStepper stepper, sunrealtype tR,
                                   N_Vector yR);

typedef int (*MRIStepInnerGetAccumulatedError)(MRIStepInnerStepper stepper,
                                               sunrealtype* accum_error);

typedef int (*MRIStepInnerResetAccumulatedError)(MRIStepInnerStepper stepper);

typedef int (*MRIStepInnerSetRTol)(MRIStepInnerStepper stepper, sunrealtype rtol);

/*---------------------------------------------------------------
  MRI coupling data structure and associated utility routines
  ---------------------------------------------------------------*/
struct MRIStepCouplingMem
{
  MRISTEP_METHOD_TYPE type; /* flag to encode the MRI method type                  */
  int nmat;         /* number of MRI coupling matrices                     */
  int stages;       /* size of coupling matrices ((stages+1) * stages)     */
  int q;            /* method order of accuracy                            */
  int p;            /* embedding order of accuracy                         */
  sunrealtype* c;   /* stage abscissae                                     */
  sunrealtype*** W; /* explicit coupling matrices [nmat][stages+1][stages] */
  sunrealtype*** G; /* implicit coupling matrices [nmat][stages+1][stages] */

  int ngroup;  /* number of stage groups (MERK-specific)              */
  int** group; /* stages to integrate together (MERK-specific)        */
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
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
void MRIStepCoupling_Space(MRIStepCoupling MRIC, sunindextype* liw,
                           sunindextype* lrw);
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

/* Creation and Reinitialization functions */
SUNDIALS_EXPORT void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, sunrealtype t0,
                                    N_Vector y0, MRIStepInnerStepper stepper,
                                    SUNContext sunctx);
SUNDIALS_EXPORT int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi,
                                  sunrealtype t0, N_Vector y0);

/* Optional input functions -- must be called AFTER MRIStepCreate */
SUNDIALS_EXPORT int MRIStepSetCoupling(void* arkode_mem, MRIStepCoupling MRIC);
SUNDIALS_EXPORT int MRIStepSetPreInnerFn(void* arkode_mem,
                                         MRIStepPreInnerFn prefn);
SUNDIALS_EXPORT int MRIStepSetPostInnerFn(void* arkode_mem,
                                          MRIStepPostInnerFn postfn);

/* Optional output functions */
SUNDIALS_EXPORT int MRIStepGetCurrentCoupling(void* arkode_mem,
                                              MRIStepCoupling* MRIC);
SUNDIALS_EXPORT int MRIStepGetLastInnerStepFlag(void* arkode_mem, int* flag);
SUNDIALS_EXPORT int MRIStepGetNumInnerStepperFails(void* arkode_mem,
                                                   long int* inner_fails);

/* Custom inner stepper functions */
SUNDIALS_EXPORT int MRIStepInnerStepper_Create(SUNContext sunctx,
                                               MRIStepInnerStepper* stepper);

SUNDIALS_EXPORT int MRIStepInnerStepper_CreateFromSUNStepper(
  SUNStepper sunstepper, MRIStepInnerStepper* stepper);

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
SUNDIALS_EXPORT int MRIStepInnerStepper_SetAccumulatedErrorGetFn(
  MRIStepInnerStepper stepper, MRIStepInnerGetAccumulatedError fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetAccumulatedErrorResetFn(
  MRIStepInnerStepper stepper, MRIStepInnerResetAccumulatedError fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_SetRTolFn(MRIStepInnerStepper stepper,
                                                  MRIStepInnerSetRTol fn);
SUNDIALS_EXPORT int MRIStepInnerStepper_AddForcing(MRIStepInnerStepper stepper,
                                                   sunrealtype t, N_Vector f);
SUNDIALS_EXPORT int MRIStepInnerStepper_GetForcingData(
  MRIStepInnerStepper stepper, sunrealtype* tshift, sunrealtype* tscale,
  N_Vector** forcing, int* nforcing);

/* --------------------------------------------------------------------------
 * Deprecated Functions -- all are superseded by shared ARKODE-level routines
 * -------------------------------------------------------------------------- */

SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeResize instead")
int MRIStepResize(void* arkode_mem, N_Vector ynew, sunrealtype t0,
                  ARKVecResizeFn resize, void* resize_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeReset instead")
int MRIStepReset(void* arkode_mem, sunrealtype tR, N_Vector yR);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSStolerances instead")
int MRIStepSStolerances(void* arkode_mem, sunrealtype reltol, sunrealtype abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSVtolerances instead")
int MRIStepSVtolerances(void* arkode_mem, sunrealtype reltol, N_Vector abstol);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWFtolerances instead")
int MRIStepWFtolerances(void* arkode_mem, ARKEwtFn efun);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinearSolver instead")
int MRIStepSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix A);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeRootInit instead")
int MRIStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDefaults instead")
int MRIStepSetDefaults(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetOrder instead")
int MRIStepSetOrder(void* arkode_mem, int ord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantType instead")
int MRIStepSetInterpolantType(void* arkode_mem, int itype);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int MRIStepSetInterpolantDegree(void* arkode_mem, int degree);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolantDegree instead")
int MRIStepSetDenseOrder(void* arkode_mem, int dord);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinearSolver instead")
int MRIStepSetNonlinearSolver(void* arkode_mem, SUNNonlinearSolver NLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNlsRhsFn instead")
int MRIStepSetNlsRhsFn(void* arkode_mem, ARKRhsFn nls_fs);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinear instead")
int MRIStepSetLinear(void* arkode_mem, int timedepend);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinear instead")
int MRIStepSetNonlinear(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNumSteps instead")
int MRIStepSetMaxNumSteps(void* arkode_mem, long int mxsteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinCRDown instead")
int MRIStepSetNonlinCRDown(void* arkode_mem, sunrealtype crdown);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinRDiv instead")
int MRIStepSetNonlinRDiv(void* arkode_mem, sunrealtype rdiv);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDeltaGammaMax instead")
int MRIStepSetDeltaGammaMax(void* arkode_mem, sunrealtype dgmax);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLSetupFrequency instead")
int MRIStepSetLSetupFrequency(void* arkode_mem, int msbp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPredictorMethod instead")
int MRIStepSetPredictorMethod(void* arkode_mem, int method);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxNonlinIters instead")
int MRIStepSetMaxNonlinIters(void* arkode_mem, int maxcor);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNonlinConvCoef instead")
int MRIStepSetNonlinConvCoef(void* arkode_mem, sunrealtype nlscoef);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetMaxHnilWarns instead")
int MRIStepSetMaxHnilWarns(void* arkode_mem, int mxhnil);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetInterpolateStopTime instead")
int MRIStepSetInterpolateStopTime(void* arkode_mem, sunbooleantype interp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStopTime instead")
int MRIStepSetStopTime(void* arkode_mem, sunrealtype tstop);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeClearStopTime instead")
int MRIStepClearStopTime(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetFixedStep instead")
int MRIStepSetFixedStep(void* arkode_mem, sunrealtype hsfixed);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetRootDirection instead")
int MRIStepSetRootDirection(void* arkode_mem, int* rootdir);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetNoInactiveRootWarn instead")
int MRIStepSetNoInactiveRootWarn(void* arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetUserData instead")
int MRIStepSetUserData(void* arkode_mem, void* user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStepFn instead")
int MRIStepSetPostprocessStepFn(void* arkode_mem, ARKPostProcessFn ProcessStep);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPostprocessStageFn instead")
int MRIStepSetPostprocessStageFn(void* arkode_mem, ARKPostProcessFn ProcessStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetStagePredictFn instead")
int MRIStepSetStagePredictFn(void* arkode_mem, ARKStagePredictFn PredictStage);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetDeduceImplicitRhs instead")
int MRIStepSetDeduceImplicitRhs(void* arkode_mem, sunbooleantype deduce);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacFn instead")
int MRIStepSetJacFn(void* arkode_mem, ARKLsJacFn jac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacEvalFrequency instead")
int MRIStepSetJacEvalFrequency(void* arkode_mem, long int msbj);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinearSolutionScaling instead")
int MRIStepSetLinearSolutionScaling(void* arkode_mem, sunbooleantype onoff);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetEpsLin instead")
int MRIStepSetEpsLin(void* arkode_mem, sunrealtype eplifac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLSNormFactor instead")
int MRIStepSetLSNormFactor(void* arkode_mem, sunrealtype nrmfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetPreconditioner instead")
int MRIStepSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacTimes instead")
int MRIStepSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetJacTimesRhsFn instead")
int MRIStepSetJacTimesRhsFn(void* arkode_mem, ARKRhsFn jtimesRhsFn);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeSetLinSysFn instead")
int MRIStepSetLinSysFn(void* arkode_mem, ARKLsLinSysFn linsys);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeEvolve instead")
int MRIStepEvolve(void* arkode_mem, sunrealtype tout, N_Vector yout,
                  sunrealtype* tret, int itask);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetDky instead")
int MRIStepGetDky(void* arkode_mem, sunrealtype t, int k, N_Vector dky);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeComputeState instead")
int MRIStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinSolvSetups instead")
int MRIStepGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int MRIStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumSteps instead")
int MRIStepGetNumSteps(void* arkode_mem, long int* nssteps);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastStep instead")
int MRIStepGetLastStep(void* arkode_mem, sunrealtype* hlast);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentTime instead")
int MRIStepGetCurrentTime(void* arkode_mem, sunrealtype* tcur);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentState instead")
int MRIStepGetCurrentState(void* arkode_mem, N_Vector* state);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetCurrentGamma instead")
int MRIStepGetCurrentGamma(void* arkode_mem, sunrealtype* gamma);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetTolScaleFactor instead")
int MRIStepGetTolScaleFactor(void* arkode_mem, sunrealtype* tolsfac);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetErrWeights instead")
int MRIStepGetErrWeights(void* arkode_mem, N_Vector eweight);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumGEvals instead")
int MRIStepGetNumGEvals(void* arkode_mem, long int* ngevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetRootInfo instead")
int MRIStepGetRootInfo(void* arkode_mem, int* rootsfound);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetUserData instead")
int MRIStepGetUserData(void* arkode_mem, void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintAllStats instead")
int MRIStepPrintAllStats(void* arkode_mem, FILE* outfile, SUNOutputFormat fmt);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetReturnFlagName instead")
char* MRIStepGetReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeWriteParameters instead")
int MRIStepWriteParameters(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "use MRIStepGetCurrentCoupling and MRIStepCoupling_Write instead")
int MRIStepWriteCoupling(void* arkode_mem, FILE* fp);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNonlinearSystemData instead")
int MRIStepGetNonlinearSystemData(void* arkode_mem, sunrealtype* tcur,
                                  N_Vector* zpred, N_Vector* z, N_Vector* F,
                                  sunrealtype* gamma, N_Vector* sdata,
                                  void** user_data);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumNonlinSolvIters instead")
int MRIStepGetNumNonlinSolvIters(void* arkode_mem, long int* nniters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumNonlinSolvConvFails instead")
int MRIStepGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNonlinSolvStats instead")
int MRIStepGetNonlinSolvStats(void* arkode_mem, long int* nniters,
                              long int* nnfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumStepSolveFails instead")
int MRIStepGetNumStepSolveFails(void* arkode_mem, long int* nncfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJac instead")
int MRIStepGetJac(void* arkode_mem, SUNMatrix* J);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJacTime instead")
int MRIStepGetJacTime(void* arkode_mem, sunrealtype* t_J);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetJacNumSteps instead")
int MRIStepGetJacNumSteps(void* arkode_mem, long* nst_J);
SUNDIALS_DEPRECATED_EXPORT_MSG(
  "Work space functions will be removed in version 8.0.0")
int MRIStepGetLinWorkSpace(void* arkode_mem, long int* lenrwLS,
                           long int* leniwLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJacEvals instead")
int MRIStepGetNumJacEvals(void* arkode_mem, long int* njevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumPrecEvals instead")
int MRIStepGetNumPrecEvals(void* arkode_mem, long int* npevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumPrecSolves instead")
int MRIStepGetNumPrecSolves(void* arkode_mem, long int* npsolves);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinIters instead")
int MRIStepGetNumLinIters(void* arkode_mem, long int* nliters);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinConvFails instead")
int MRIStepGetNumLinConvFails(void* arkode_mem, long int* nlcfails);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJTSetupEvals instead")
int MRIStepGetNumJTSetupEvals(void* arkode_mem, long int* njtsetups);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumJtimesEvals instead")
int MRIStepGetNumJtimesEvals(void* arkode_mem, long int* njvevals);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumLinRhsEvals instead")
int MRIStepGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLastLinFlag instead")
int MRIStepGetLastLinFlag(void* arkode_mem, long int* flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetLinReturnFlagName instead")
char* MRIStepGetLinReturnFlagName(long int flag);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeFree instead")
void MRIStepFree(void** arkode_mem);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodePrintMem instead")
void MRIStepPrintMem(void* arkode_mem, FILE* outfile);
SUNDIALS_DEPRECATED_EXPORT_MSG("use ARKodeGetNumRhsEvals instead")
int MRIStepGetNumRhsEvals(void* arkode_mem, long int* nfse_evals,
                          long int* nfsi_evals);

#ifdef __cplusplus
}
#endif

#endif
