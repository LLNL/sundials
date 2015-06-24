/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for the optional input and 
 * output functions for the ARKODE solver.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


/*===============================================================
 ARKODE optional input functions
===============================================================*/

/*---------------------------------------------------------------
 ARKodeSetDefaults:

 Resets all optional inputs to ARKode default values.  Does not 
 change problem-defining function pointers fe and fi or 
 user_data pointer.  Also leaves alone any data 
 structures/options related to root-finding (those can be reset 
 using ARKodeRootInit).
---------------------------------------------------------------*/
int ARKodeSetDefaults(void *arkode_mem)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set default values for integrator optional inputs */
  ark_mem->ark_q                = Q_DEFAULT;      /* method order */
  ark_mem->ark_p                = 0;              /* embedding order */
  ark_mem->ark_dense_q          = QDENSE_DEF;     /* dense output order */
  ark_mem->ark_expstab          = arkExpStab;     /* explicit stability fn */
  ark_mem->ark_estab_data       = ark_mem;        /* explicit stability data */
  ark_mem->ark_fixedstep        = FALSE;          /* default to use adaptive steps */
  ark_mem->ark_hadapt           = NULL;           /* step adaptivity fn */
  ark_mem->ark_hadapt_data      = NULL;           /* step adaptivity data */
  ark_mem->ark_hadapt_cfl       = CFLFAC;         /* explicit stability factor */
  ark_mem->ark_hadapt_safety    = SAFETY;         /* step adaptivity safety factor  */
  ark_mem->ark_hadapt_bias      = BIAS;           /* step adaptivity error bias */
  ark_mem->ark_hadapt_growth    = GROWTH;         /* step adaptivity growth factor */
  ark_mem->ark_hadapt_lbound    = HFIXED_LB;      /* step adaptivity no-change lower bound */
  ark_mem->ark_hadapt_ubound    = HFIXED_UB;      /* step adaptivity no-change upper bound */
  ark_mem->ark_hadapt_pq        = FALSE;          /* use embedding order */
  ark_mem->ark_hadapt_imethod   = 0;              /* PID controller */
  ark_mem->ark_hadapt_k1        = AD0_K1;         /* step adaptivity parameter */
  ark_mem->ark_hadapt_k2        = AD0_K2;         /* step adaptivity parameter */
  ark_mem->ark_hadapt_k3        = AD0_K3;         /* step adaptivity parameter */
  ark_mem->ark_predictor        = 3;              /* max order close, first order far */
  ark_mem->ark_reltol           = 1.e-4;          /* relative tolerance */
  ark_mem->ark_itol             = ARK_SS;         /* scalar-scalar solution tolerances */
  ark_mem->ark_ritol            = ARK_SS;         /* scalar-scalar residual tolerances */
  ark_mem->ark_Sabstol          = 1.e-9;          /* solution absolute tolerance */
  ark_mem->ark_SRabstol         = 1.e-9;          /* residual absolute tolerance */
  ark_mem->ark_user_efun        = FALSE;          /* no user-supplied ewt function */
  ark_mem->ark_efun             = arkEwtSet;      /* built-in ewt function */
  ark_mem->ark_e_data           = NULL;           /* ewt function data */
  ark_mem->ark_user_rfun        = FALSE;          /* no user-supplied rwt function */
  ark_mem->ark_rfun             = arkRwtSet;      /* built-in rwt function */
  ark_mem->ark_e_data           = NULL;           /* rwt function data */
  ark_mem->ark_linear           = FALSE;          /* nonlinear problem */
  ark_mem->ark_linear_timedep   = TRUE;           /* dfi/dy depends on t */
  ark_mem->ark_explicit         = FALSE;          /* fi(t,y) will be used */
  ark_mem->ark_implicit         = FALSE;          /* fe(t,y) will be used */
  ark_mem->ark_ehfun            = arkErrHandler;  /* default error handler fn */
  ark_mem->ark_eh_data          = ark_mem;        /* error handler data */
  ark_mem->ark_errfp            = stderr;         /* output stream for errors */
  ark_mem->ark_mxstep           = MXSTEP_DEFAULT; /* max number of steps */
  ark_mem->ark_mxhnil           = MXHNIL;         /* max warns of t+h==t */
  ark_mem->ark_hin              = ZERO;           /* determine initial step on-the-fly */
  ark_mem->ark_hmin             = ZERO;           /* no minimum step size */
  ark_mem->ark_hmax_inv         = ZERO;           /* no maximum step size */
  ark_mem->ark_tstopset         = FALSE;          /* no stop time set */
  ark_mem->ark_tstop            = ZERO;           /* no fixed stop time */
  ark_mem->ark_maxcor           = MAXCOR;         /* max nonlinear iters/stage */
  ark_mem->ark_maxnef           = MAXNEF;         /* max error test fails */
  ark_mem->ark_maxncf           = MAXNCF;         /* max convergence fails */
  ark_mem->ark_nlscoef          = NLSCOEF;        /* nonlinear tolerance coefficient */
  ark_mem->ark_etamx1           = ETAMX1;         /* max change on first step */
  ark_mem->ark_etamxf           = ETAMXF;         /* max change on error-failed step */
  ark_mem->ark_small_nef        = SMALL_NEF;      /* num error fails before ETAMXF enforced */
  ark_mem->ark_etacf            = ETACF;          /* max change on convergence failure */
  ark_mem->ark_crdown           = CRDOWN;         /* nonlinear convergence estimate coeff. */
  ark_mem->ark_rdiv             = RDIV;           /* nonlinear divergence tolerance */
  ark_mem->ark_dgmax            = DGMAX;          /* max step change before recomputing J or P */
  ark_mem->ark_msbp             = MSBP;           /* max steps between updates to J or P */
  ark_mem->ark_use_fp           = FALSE;          /* use Newton solver */
  ark_mem->ark_fp_m             = FP_ACCEL_M;     /* num Anderson acceleration vectors */
  ark_mem->ark_diagfp           = NULL;           /* no solver diagnostics file */
  ark_mem->ark_report           = FALSE;          /* don't report solver diagnostics */
  ark_mem->ark_stages           = 0;              /* no stages */
  ark_mem->ark_istage           = 0;              /* current stage */
  for (i=0; i<ARK_S_MAX; i++) {                   /* no Butcher table */
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(ark_mem->ark_Ae,i,j) = ZERO;
      ARK_A(ark_mem->ark_Ai,i,j) = ZERO;
    }
    ark_mem->ark_c[i]   = ZERO;
    ark_mem->ark_b[i]   = ZERO;
    ark_mem->ark_b2[i]  = ZERO;
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetOptimalParams:

 Sets all adaptivity and solver parameters to our 'best guess' 
 values, for a given integration method (ERK, DIRK, ARK), a 
 given method order, and a given nonlinear solver type.  Should
 only be called after the method order, solver, and integration
 method have been set.
---------------------------------------------------------------*/
int ARKodeSetOptimalParams(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetOptimalParams", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Choose values based on method, order */

  /*    explicit */
  if (ark_mem->ark_explicit) {
    ark_mem->ark_hadapt_imethod = 1;
    ark_mem->ark_hadapt_safety  = RCONST(0.99);
    ark_mem->ark_hadapt_bias    = RCONST(1.2);
    ark_mem->ark_hadapt_growth  = RCONST(25.0);
    ark_mem->ark_hadapt_k1      = RCONST(0.8);
    ark_mem->ark_hadapt_k2      = RCONST(0.31);
    ark_mem->ark_etamxf         = RCONST(0.3);

  /*    implicit */
  } else if (ark_mem->ark_implicit) {
    switch (ark_mem->ark_q) {
    case 2:   /* just use standard defaults since better ones unknown */
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = SAFETY;
      ark_mem->ark_hadapt_bias      = BIAS;
      ark_mem->ark_hadapt_growth    = GROWTH;
      ark_mem->ark_etamxf           = ETAMXF;
      ark_mem->ark_nlscoef          = RCONST(0.001);
      ark_mem->ark_maxcor           = 5;
      ark_mem->ark_crdown           = CRDOWN;
      ark_mem->ark_rdiv             = RDIV;
      ark_mem->ark_dgmax            = DGMAX;
      ark_mem->ark_msbp             = MSBP;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 3:
      ark_mem->ark_hadapt_imethod   = 2;
      ark_mem->ark_hadapt_safety    = RCONST(0.957);
      ark_mem->ark_hadapt_bias      = RCONST(1.9);
      ark_mem->ark_hadapt_growth    = RCONST(17.6);
      ark_mem->ark_etamxf           = RCONST(0.45);
      ark_mem->ark_nlscoef          = RCONST(0.22);
      ark_mem->ark_crdown           = RCONST(0.17);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.19);
      ark_mem->ark_msbp             = 60;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 4:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.988);
      ark_mem->ark_hadapt_bias      = RCONST(1.2);
      ark_mem->ark_hadapt_growth    = RCONST(31.5);
      ark_mem->ark_hadapt_k1        = RCONST(0.535);
      ark_mem->ark_hadapt_k2        = RCONST(0.209);
      ark_mem->ark_hadapt_k3        = RCONST(0.148);
      ark_mem->ark_etamxf           = RCONST(0.33);
      ark_mem->ark_nlscoef          = RCONST(0.24);
      ark_mem->ark_crdown           = RCONST(0.26);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.16);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 5:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.937);
      ark_mem->ark_hadapt_bias      = RCONST(3.3);
      ark_mem->ark_hadapt_growth    = RCONST(22.0);
      ark_mem->ark_hadapt_k1        = RCONST(0.56);
      ark_mem->ark_hadapt_k2        = RCONST(0.338);
      ark_mem->ark_hadapt_k3        = RCONST(0.14);
      ark_mem->ark_etamxf           = RCONST(0.44);
      ark_mem->ark_nlscoef          = RCONST(0.25);
      ark_mem->ark_crdown           = RCONST(0.4);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.32);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    }

    /* newton vs fixed-point */
    if (ark_mem->ark_use_fp)  ark_mem->ark_maxcor = 10;

  /*    imex */
  } else {
    switch (ark_mem->ark_q) {
    case 3:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.965);
      ark_mem->ark_hadapt_bias      = RCONST(1.42);
      ark_mem->ark_hadapt_growth    = RCONST(28.7);
      ark_mem->ark_hadapt_k1        = RCONST(0.54);
      ark_mem->ark_hadapt_k2        = RCONST(0.36);
      ark_mem->ark_hadapt_k3        = RCONST(0.14);
      ark_mem->ark_etamxf           = RCONST(0.46);
      ark_mem->ark_nlscoef          = RCONST(0.22);
      ark_mem->ark_crdown           = RCONST(0.17);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.19);
      ark_mem->ark_msbp             = 60;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 4:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.97);
      ark_mem->ark_hadapt_bias      = RCONST(1.35);
      ark_mem->ark_hadapt_growth    = RCONST(25.0);
      ark_mem->ark_hadapt_k1        = RCONST(0.543);
      ark_mem->ark_hadapt_k2        = RCONST(0.297);
      ark_mem->ark_hadapt_k3        = RCONST(0.14);
      ark_mem->ark_etamxf           = RCONST(0.47);
      ark_mem->ark_nlscoef          = RCONST(0.24);
      ark_mem->ark_crdown           = RCONST(0.26);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.16);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 5:
      ark_mem->ark_hadapt_imethod   = 1;
      ark_mem->ark_hadapt_safety    = RCONST(0.993);
      ark_mem->ark_hadapt_bias      = RCONST(1.15);
      ark_mem->ark_hadapt_growth    = RCONST(28.5);
      ark_mem->ark_hadapt_k1        = RCONST(0.8);
      ark_mem->ark_hadapt_k2        = RCONST(0.35);
      ark_mem->ark_etamxf           = RCONST(0.3);
      ark_mem->ark_nlscoef          = RCONST(0.25);
      ark_mem->ark_crdown           = RCONST(0.4);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.32);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    }

    /* newton vs fixed-point */
    if (ark_mem->ark_use_fp)  ark_mem->ark_maxcor = 10;

  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetErrHandlerFn:

 Specifies the error handler function
---------------------------------------------------------------*/
int ARKodeSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun, 
			  void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetErrHandlerFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set user-provided values, or defaults, depending on argument */
  if (ehfun == NULL) {
    ark_mem->ark_ehfun   = arkErrHandler;
    ark_mem->ark_eh_data = ark_mem;
  } else {
    ark_mem->ark_ehfun   = ehfun;
    ark_mem->ark_eh_data = eh_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetErrFile:

 Specifies the FILE pointer for output (NULL means no messages)
---------------------------------------------------------------*/
int ARKodeSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetErrFile", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_errfp = errfp;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetUserData:

 Specifies the user data pointer for f
---------------------------------------------------------------*/
int ARKodeSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetUserData", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_user_data = user_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetDiagnostics:

 Specifies to enable solver diagnostics, and specifies the FILE
 pointer for output (diagfp==NULL disables output)
---------------------------------------------------------------*/
int ARKodeSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDiagnostics", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_diagfp = diagfp;
  if (diagfp != NULL) {
    ark_mem->ark_report = TRUE;
  } else {
    ark_mem->ark_report = FALSE;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetOrder:

 Specifies the method order

 ** Note in documentation that this should not be called along 
    with ARKodeSetERKTable, ARKodeSetIRKTable, ARKodeSetARKTable, 
    ARKodeSetERKTableNum, ARKodeSetIRKTableNum or 
    ARKodeSetARKTableNum.  This routine is used to specify a 
    desired method order using default Butcher tables, whereas 
    any user-supplied table will have their own order associated 
    with them.
---------------------------------------------------------------*/
int ARKodeSetOrder(void *arkode_mem, int ord)
{
  int i, j;

  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetOrder", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) {
    ark_mem->ark_q = Q_DEFAULT;
  } else {
    ark_mem->ark_q = ord;
  }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  ark_mem->ark_stages = 0;
  ark_mem->ark_istage = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(ark_mem->ark_Ae,i,j) = ZERO;
      ARK_A(ark_mem->ark_Ai,i,j) = ZERO;
    }
    ark_mem->ark_c[i]   = ZERO;
    ark_mem->ark_b[i]   = ZERO;
    ark_mem->ark_b2[i]  = ZERO;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetDenseOrder:

 Specifies the polynomial order for dense output.  Allowed values
 range from 0 to min(q,5), where q is the order of the time 
 integration method.  Illegal values imply to use the default.
---------------------------------------------------------------*/
int ARKodeSetDenseOrder(void *arkode_mem, int dord)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDenseOrder", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check input */
  if (dord > 5) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetDenseOrder", "Dense output order must be <= 5");
    return(ARK_ILL_INPUT);
  }
  /* NOTE: we check that dord < q internally, to allow for subsequent 
     changes via ARKodeSetOrder */

  /* set user-provided value, or default, depending on argument */
  /* if ((dord < 0) || (dord > 5)) { */
  if ((dord < 0) || (dord > 3)) {
    ark_mem->ark_dense_q = QDENSE_DEF;
  } else {
    ark_mem->ark_dense_q = dord;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetLinear:

 Specifies that the implicit portion of the problem is linear, 
 and to tighten the linear solver tolerances while taking only 
 one Newton iteration.  Not useful when used in combination with
 the fixed-point solver.  Automatically tightens DeltaGammaMax 
 to ensure that step size changes cause Jacobian recomputation.

 The argument should be 1 or 0, where 1 indicates that the 
 Jacobian of fi with respect to y depends on time, and
 0 indicates that it is not time dependent.  Alternately, when 
 using an iterative linear solver this flag denotes time 
 dependence of the preconditioner. 
---------------------------------------------------------------*/
int ARKodeSetLinear(void *arkode_mem, int timedepend)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetLinear", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_linear = TRUE;
  ark_mem->ark_linear_timedep = (timedepend == 1);
  ark_mem->ark_dgmax = RCONST(100.0)*UNIT_ROUNDOFF;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNonlinear:

 Specifies that the implicit portion of the problem is nonlinear.
 Used to undo a previous call to ARKodeSetLinear.  Automatically
 loosens DeltaGammaMax back to default value.
---------------------------------------------------------------*/
int ARKodeSetNonlinear(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNonlinear", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_linear = FALSE;
  ark_mem->ark_linear_timedep = TRUE;
  ark_mem->ark_dgmax = DGMAX;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetFixedPoint:

 Specifies that the implicit portion of the problem should be 
 solved using the accelerated fixed-point solver instead of the 
 Newton iteration.  Allowed values for the dimension of the 
 acceleration space, fp_m, must be non-negative.  Illegal 
 values imply to use the default.
---------------------------------------------------------------*/
int ARKodeSetFixedPoint(void *arkode_mem, long int fp_m)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetFixedPoint", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_use_fp = TRUE;
  if (fp_m < 0) {
    ark_mem->ark_fp_m = FP_ACCEL_M;
  } else {
    ark_mem->ark_fp_m = fp_m;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNewton:

 Specifies that the implicit portion of the problem should be 
 solved using the modified Newton solver.  Used to undo a 
 previous call to ARKodeSetFixedPoint.
---------------------------------------------------------------*/
int ARKodeSetNewton(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNewton", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_use_fp = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetExplicit:

 Specifies that the implicit portion of the problem is disabled, 
 and to use an explicit RK method.
---------------------------------------------------------------*/
int ARKodeSetExplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetExplicit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fe is defined */
  if (ark_mem->ark_fe == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetExplicit", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_explicit = TRUE;
  ark_mem->ark_implicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetImplicit:

 Specifies that the explicit portion of the problem is disabled, 
 and to use an implicit RK method.
---------------------------------------------------------------*/
int ARKodeSetImplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetImplicit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fi is defined */
  if (ark_mem->ark_fi == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetImplicit", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_implicit = TRUE;
  ark_mem->ark_explicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetImEx:

 Specifies that the specifies that problem has both implicit and
 explicit parts, and to use an ARK method (this is the default).
---------------------------------------------------------------*/
int ARKodeSetImEx(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetImEx", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fe and fi are defined */
  if (ark_mem->ark_fe == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetImEx", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }
  if (ark_mem->ark_fi == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetImEx", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_explicit = FALSE;
  ark_mem->ark_implicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetERKTable:

 Specifies to use a customized Butcher table for the explicit 
 portion of the system (automatically calls ARKodeSetExplicit).
---------------------------------------------------------------*/
int ARKodeSetERKTable(void *arkode_mem, int s, int q, int p,
		      realtype *c, realtype *A, realtype *b, 
		      realtype *bembed)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > ARK_S_MAX) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERKTable", "s exceeds ARK_S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (A == NULL) || (b == NULL) || (bembed == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters */
  ark_mem->ark_stages = 0;
  ark_mem->ark_q = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    ark_mem->ark_c[i]  = 0.0;
    ark_mem->ark_b[i]  = 0.0;
    ark_mem->ark_b2[i] = 0.0;
    for (j=0; j<ARK_S_MAX; j++) 
      ARK_A(ark_mem->ark_Ae,i,j) = 0.0;
  }

  /* set the relevant parameters */
  ark_mem->ark_stages = s;
  ark_mem->ark_q = q;
  ark_mem->ark_p = p;
  for (i=0; i<s; i++) {
    ark_mem->ark_c[i]  = c[i];
    ark_mem->ark_b[i]  = b[i];
    ark_mem->ark_b2[i] = bembed[i];
    for (j=0; j<s; j++) {
      ARK_A(ark_mem->ark_Ae,i,j) = A[i*s + j];
    }
  }
  
  /* set method as purely explicit */
  if (ARKodeSetExplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetIRKTable:

 Specifies to use a customized Butcher table for the implicit 
 portion of the system (automatically calls ARKodeSetImplicit).
---------------------------------------------------------------*/
int ARKodeSetIRKTable(void *arkode_mem, int s, int q, int p,
		      realtype *c, realtype *A, realtype *b, 
		      realtype *bembed)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > ARK_S_MAX) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRKTable", "s exceeds ARK_S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (A == NULL) || (b == NULL) || (bembed == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters */
  ark_mem->ark_stages = 0;
  ark_mem->ark_q = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    ark_mem->ark_c[i]  = 0.0;
    ark_mem->ark_b[i]  = 0.0;
    ark_mem->ark_b2[i] = 0.0;
    for (j=0; j<ARK_S_MAX; j++) 
      ARK_A(ark_mem->ark_Ai,i,j) = 0.0;
  }

  /* set the relevant parameters */
  ark_mem->ark_stages = s;
  ark_mem->ark_q = q;
  ark_mem->ark_p = p;
  for (i=0; i<s; i++) {
    ark_mem->ark_c[i]  = c[i];
    ark_mem->ark_b[i]  = b[i];
    ark_mem->ark_b2[i] = bembed[i];
    for (j=0; j<s; j++) {
      ARK_A(ark_mem->ark_Ai,i,j) = A[i*s + j];
    }
  }

  /* set method as purely implicit */
  if (ARKodeSetImplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetARKTables:

 Specifies to use customized Butcher tables for the ImEx system
 (automatically calls ARKodeSetImEx).
---------------------------------------------------------------*/
int ARKodeSetARKTables(void *arkode_mem, int s, int q, int p,
		       realtype *c, realtype *Ai, realtype *Ae, 
		       realtype *b, realtype *bembed)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > ARK_S_MAX) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTables", "s exceeds ARK_S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (Ai == NULL) || (Ae == NULL) || 
      (b == NULL) || (bembed == NULL)) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters */
  ark_mem->ark_stages = 0;
  ark_mem->ark_q = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    ark_mem->ark_c[i]  = 0.0;
    ark_mem->ark_b[i]  = 0.0;
    ark_mem->ark_b2[i] = 0.0;
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(ark_mem->ark_Ai,i,j) = 0.0;
      ARK_A(ark_mem->ark_Ae,i,j) = 0.0;
    }
  }

  /* set the relevant parameters */
  ark_mem->ark_stages = s;
  ark_mem->ark_q = q;
  ark_mem->ark_p = p;
  for (i=0; i<s; i++) {
    ark_mem->ark_c[i]  = c[i];
    ark_mem->ark_b[i]  = b[i];
    ark_mem->ark_b2[i] = bembed[i];
    for (j=0; j<s; j++) {
      ARK_A(ark_mem->ark_Ai,i,j) = Ai[i*s + j];
      ARK_A(ark_mem->ark_Ae,i,j) = Ae[i*s + j];
    }
  }

  /* set method as ImEx */
  if (ARKodeSetImEx(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTables", MSGARK_MISSING_F);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetERKTableNum:

 Specifies to use a pre-existing Butcher table for the explicit 
 portion of the problem, based on the integer flag held in 
 ARKodeLoadButcherTable() within the file arkode_butcher.c 
 (automatically calls ARKodeSetExplicit).
---------------------------------------------------------------*/
int ARKodeSetERKTableNum(void *arkode_mem, int itable)
{
  int iflag;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTableNum", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check that argument specifies an explicit table (0-12) */
  if (itable<0 || itable>10) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTableNum", 
		    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* fill in table based on argument */
  iflag = ARKodeLoadButcherTable(itable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ae, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);
  if (iflag != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTableNum", 
		    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }

  /* set method as purely explicit */
  if (ARKodeSetExplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERKTableNum", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetIRKTableNum:

 Specifies to use a pre-existing Butcher table for the implicit 
 portion of the problem, based on the integer flag held in 
 ARKodeLoadButcherTable() within the file arkode_butcher.c
 (automatically calls ARKodeSetImplicit).
---------------------------------------------------------------*/
int ARKodeSetIRKTableNum(void *arkode_mem, int itable)
{
  int iflag;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTableNum", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check that argument specifies an implicit table (13-27) */
  if (itable<11) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTableNum", 
		    "Illegal IRK table number");
    return(ARK_ILL_INPUT);
  }

  /* fill in table based on argument */
  iflag = ARKodeLoadButcherTable(itable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ai, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);
  if (iflag != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTableNum", 
		    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }

  /* set method as purely implicit */
  if (ARKodeSetImplicit(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRKTableNum", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetARKTableNum:

 Specifies to use pre-existing Butcher tables for the ImEx system,
 based on the integer flags held in ARKodeLoadButcherTable() 
 within the file arkode_butcher.c (automatically calls ARKodeSetImEx).
---------------------------------------------------------------*/
int ARKodeSetARKTableNum(void *arkode_mem, int itable, int etable)
{
  int iflag, eflag;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTableNum", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that tables match */
  iflag = 1;
  if ((etable == 2) && (itable == 15))  iflag = 0;
  if ((etable == 4) && (itable == 20))  iflag = 0;
  if ((etable == 9) && (itable == 22))  iflag = 0;
  if (iflag) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTableNum", 
		    "Incompatible Butcher tables for ARK method");
    return(ARK_ILL_INPUT);
  }

  /* fill in tables based on arguments */
  iflag = ARKodeLoadButcherTable(itable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ai, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);
  eflag = ARKodeLoadButcherTable(etable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ae, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);

  /* check that requested tables are legal */
  if (iflag != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTableNum", 
		    "Illegal IRK table number");
    return(ARK_ILL_INPUT);
  }
  if (eflag != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTableNum", 
		    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* set method as ImEx */
  if (ARKodeSetImEx(arkode_mem) != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTableNum", MSGARK_MISSING_F);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxNumSteps:

 Specifies the maximum number of integration steps
---------------------------------------------------------------*/
int ARKodeSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */
  if (mxsteps == 0)
    ark_mem->ark_mxstep = MXSTEP_DEFAULT;
  else
    ark_mem->ark_mxstep = mxsteps;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxHnilWarns:

 Specifies the maximum number of warnings for small h
---------------------------------------------------------------*/
int ARKodeSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxHnilWarns", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing mxhnil=0 sets the default, otherwise use input. */
  if (mxhnil == 0) {
    ark_mem->ark_mxhnil = 10;
  } else {
    ark_mem->ark_mxhnil = mxhnil;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetInitStep:

 Specifies the initial step size
---------------------------------------------------------------*/
int ARKodeSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetInitStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing hin=0 sets the default, otherwise use input. */
  if (hin == ZERO) {
    ark_mem->ark_hin = ZERO;
  } else {
    ark_mem->ark_hin = hin;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMinStep:

 Specifies the minimum step size
---------------------------------------------------------------*/
int ARKodeSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMinStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing a value <= 0 sets hmax = infinity */
  if (hmin <= ZERO) {
    ark_mem->ark_hmin = ZERO;
    return(ARK_SUCCESS);
  }

  /* check that hmin and hmax are agreeable */
  if (hmin * ark_mem->ark_hmax_inv > ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMinStep", MSGARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->ark_hmin = hmin;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxStep:

 Specifies the maximum step size
---------------------------------------------------------------*/
int ARKodeSetMaxStep(void *arkode_mem, realtype hmax)
{
  realtype hmax_inv;
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxStep", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing a value <= 0 sets hmax = infinity */
  if (hmax <= ZERO) {
    ark_mem->ark_hmax_inv = ZERO;
    return(ARK_SUCCESS);
  }

  /* check that hmax and hmin are agreeable */
  hmax_inv = ONE/hmax;
  if (hmax_inv * ark_mem->ark_hmin > ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMaxStep", MSGARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->ark_hmax_inv = hmax_inv;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetStopTime:

 Specifies the time beyond which the integration is not to proceed.
---------------------------------------------------------------*/
int ARKodeSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetStopTime", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* If ARKode was called at least once, test if tstop is legal
     (i.e. if it was not already passed).
     If ARKodeSetStopTime is called before the first call to ARKode,
     tstop will be checked in ARKode. */
  if (ark_mem->ark_nst > 0) {
    if ( (tstop - ark_mem->ark_tn) * ark_mem->ark_h < ZERO ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeSetStopTime", MSGARK_BAD_TSTOP, 
		      tstop, ark_mem->ark_tn);
      return(ARK_ILL_INPUT);
    }
  }

  ark_mem->ark_tstop    = tstop;
  ark_mem->ark_tstopset = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetFixedStep:

 Specifies to use a fixed time step size instead of performing 
 any form of temporal adaptivity.  ARKode will use this step size 
 for all steps (unless tstop is set, in which case it may need to 
 modify that last step approaching tstop.  If any (non)linear
 solver failure occurs, ARKode will immediately return with an 
 error message since the time step size cannot be modified.  

 Any nonzero argument will result in the use of that fixed step 
 size; an argument of 0 will re-enable temporal adaptivity.
---------------------------------------------------------------*/
int ARKodeSetFixedStep(void *arkode_mem, realtype hfixed)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetFixedStep", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set ark_mem entry */
  if (hfixed != ZERO) {
    ark_mem->ark_fixedstep = TRUE;
    ark_mem->ark_hin = hfixed;
  } else {
    ark_mem->ark_fixedstep = FALSE;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetCFLFraction:

 Specifies the safety factor to use on the maximum explicitly-
 stable step size.  Allowable values must be within the open 
 interval (0,1).  A non-positive input implies a reset to
 the default value.
---------------------------------------------------------------*/
int ARKodeSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetCFLFraction", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for allowable parameters */
  if (cfl_frac >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetCFLFraction", "Illegal CFL fraction");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters into ark_mem, 
     otherwise set default */
  if (cfl_frac <= ZERO) {
    ark_mem->ark_hadapt_cfl = CFLFAC;
  } else {
    ark_mem->ark_hadapt_cfl = cfl_frac;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetSafetyFactor:

 Specifies the safety factor to use on the error-based predicted 
 time step size.  Allowable values must be within the open 
 interval (0,1).  A non-positive input implies a reset to the 
 default value.
---------------------------------------------------------------*/
int ARKodeSetSafetyFactor(void *arkode_mem, realtype safety)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetSafetyFactoy", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for allowable parameters */
  if (safety >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetSafetyFactor", "Illegal safety factor");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters into ark_mem, 
     otherwise set default */
  if (safety <= ZERO) {
    ark_mem->ark_hadapt_safety = SAFETY;
  } else {
    ark_mem->ark_hadapt_safety = safety;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetErrorBias:

 Specifies the error bias to use when performing adaptive-step
 error control.  Allowable values must be >= 1.0.  Any illegal
 value implies a reset to the default value.
---------------------------------------------------------------*/
int ARKodeSetErrorBias(void *arkode_mem, realtype bias)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetErrorBias", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set allowed value into ark_mem, otherwise set default */
  if (bias < 1.0) {
    ark_mem->ark_hadapt_bias = BIAS;
  } else {
    ark_mem->ark_hadapt_bias = bias;
  }
    
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxGrowth:

 Specifies the maximum step size growth factor to be allowed
 between successive integration steps.  Note: the first step uses 
 a separate maximum growth factor.  Allowable values must be 
 > 1.0.  Any illegal value implies a reset to the default.
---------------------------------------------------------------*/
int ARKodeSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxGrowth", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set allowed value into ark_mem, otherwise set default */
  if (mx_growth == ZERO) {
    ark_mem->ark_hadapt_growth = GROWTH;
  } else {
    ark_mem->ark_hadapt_growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetFixedStepBounds:

 Specifies the step size growth interval within which the step 
 size will remain unchanged.  Allowable values must enclose the 
 value 1.0.  Any illegal interval implies a reset to the default.
---------------------------------------------------------------*/
int ARKodeSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetFixedStepBounds", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set allowable interval into ark_mem, otherwise set defaults */
  if ((lb <= 1.0) && (ub >= 1.0)) {
    ark_mem->ark_hadapt_lbound = lb;
    ark_mem->ark_hadapt_ubound = ub;
  } else {
    ark_mem->ark_hadapt_lbound = HFIXED_LB;
    ark_mem->ark_hadapt_ubound = HFIXED_UB;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetAdaptivityMethod:

 Specifies the built-in time step adaptivity algorithm (and 
 optionally, its associated parameters) to use.  All parameters 
 will be checked for validity when used by the solver.
---------------------------------------------------------------*/
int ARKodeSetAdaptivityMethod(void *arkode_mem, int imethod, 
			      int idefault, int pq, 
			      realtype *adapt_params)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetAdaptivityMethod", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for allowable parameters */
  if ((imethod > 5) || (imethod < 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* set adaptivity method */
  ark_mem->ark_hadapt_imethod = imethod;

  /* set flag whether to use p or q */
  ark_mem->ark_hadapt_pq = (pq != 0);

  /* set method parameters */
  if (idefault == 1) {
    switch (ark_mem->ark_hadapt_imethod) {
    case (0):
      ark_mem->ark_hadapt_k1 = AD0_K1; 
      ark_mem->ark_hadapt_k2 = AD0_K2;
      ark_mem->ark_hadapt_k3 = AD0_K3; break;
    case (1):
      ark_mem->ark_hadapt_k1 = AD1_K1;
      ark_mem->ark_hadapt_k2 = AD1_K2; break;
    case (2):
      ark_mem->ark_hadapt_k1 = AD2_K1; break;
    case (3):
      ark_mem->ark_hadapt_k1 = AD3_K1;
      ark_mem->ark_hadapt_k2 = AD3_K2; break;
    case (4):
      ark_mem->ark_hadapt_k1 = AD4_K1;
      ark_mem->ark_hadapt_k2 = AD4_K2; break;
    case (5):
      ark_mem->ark_hadapt_k1 = AD5_K1;
      ark_mem->ark_hadapt_k2 = AD5_K2; 
      ark_mem->ark_hadapt_k3 = AD5_K3; break;
    }
  } else {
    ark_mem->ark_hadapt_k1 = adapt_params[0];
    ark_mem->ark_hadapt_k2 = adapt_params[1];
    ark_mem->ark_hadapt_k3 = adapt_params[2];
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetAdaptivityFn:

 Specifies the user-provided time step adaptivity function to use.
---------------------------------------------------------------*/
int ARKodeSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun,
			  void *h_data)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetAdaptivityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  if (hfun == NULL) {
    ark_mem->ark_hadapt         = NULL;
    ark_mem->ark_hadapt_data    = NULL;
    ark_mem->ark_hadapt_imethod = 0;
  } else {
    ark_mem->ark_hadapt         = hfun;
    ark_mem->ark_hadapt_data    = h_data;
    ark_mem->ark_hadapt_imethod = -1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxFirstGrowth:

 Specifies the user-provided time step adaptivity constant 
 etamx1.  Legal values are greater than 1.0.  Illegal values 
 imply a reset to the default value. 
---------------------------------------------------------------*/
int ARKodeSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxFirstGrowth", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    ark_mem->ark_etamx1 = ETAMX1;
  } else {
    ark_mem->ark_etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxEFailGrowth:

 Specifies the user-provided time step adaptivity constant 
 etamxf. Legal values are in the interval (0,1].  Illegal values 
 imply a reset to the default value. 
---------------------------------------------------------------*/
int ARKodeSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxEFailGrowth", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    ark_mem->ark_etamxf = ETAMXF;
  } else {
    ark_mem->ark_etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetSmallNumEFails:

 Specifies the user-provided time step adaptivity constant
 small_nef.  Legal values are > 0.  Illegal values 
 imply a reset to the default value. 
---------------------------------------------------------------*/
int ARKodeSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetSmallNumEFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    ark_mem->ark_small_nef = SMALL_NEF;
  } else {
    ark_mem->ark_small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxCFailGrowth:

 Specifies the user-provided time step adaptivity constant
 etacf. Legal values are in the interval (0,1].  Illegal values 
 imply a reset to the default value. 
---------------------------------------------------------------*/
int ARKodeSetMaxCFailGrowth(void *arkode_mem, realtype etacf)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxCFailGrowth", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) {
    ark_mem->ark_etacf = ETACF;
  } else {
    ark_mem->ark_etacf = etacf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNonlinCRDown:

 Specifies the user-provided nonlinear convergence constant
 crdown.  Legal values are strictly positive; illegal values 
 imply a reset to the default.
---------------------------------------------------------------*/
int ARKodeSetNonlinCRDown(void *arkode_mem, realtype crdown)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNonlinCRDown", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if (crdown <= ZERO) {
    ark_mem->ark_crdown = CRDOWN;
  } else {
    ark_mem->ark_crdown = crdown;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNonlinRDiv:

 Specifies the user-provided nonlinear convergence constant
 rdiv.  Legal values are strictly positive; illegal values 
 imply a reset to the default.
---------------------------------------------------------------*/
int ARKodeSetNonlinRDiv(void *arkode_mem, realtype rdiv)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNonlinRDiv", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if (rdiv <= ZERO) {
    ark_mem->ark_rdiv = RDIV;
  } else {
    ark_mem->ark_rdiv = rdiv;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetDeltaGammaMax:

 Specifies the user-provided linear setup decision constant
 dgmax.  Legal values are strictly positive; illegal values imply 
 a reset to the default. 
---------------------------------------------------------------*/
int ARKodeSetDeltaGammaMax(void *arkode_mem, realtype dgmax)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDeltaGammaMax", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if (dgmax <= ZERO) {
    ark_mem->ark_dgmax = DGMAX;
  } else {
    ark_mem->ark_dgmax = dgmax;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxStepsBetweenLSet:

 Specifies the user-provided linear setup decision constant
 msbp.  Positive values give the number of time steps to wait 
 before calling lsetup; negative values imply recomputation of 
 lsetup at each Newton iteration; a zero value implies a reset 
 to the default. 
---------------------------------------------------------------*/
int ARKodeSetMaxStepsBetweenLSet(void *arkode_mem, int msbp)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxStepsBetweenLSet", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) {
    ark_mem->ark_msbp = MSBP;
  } else {
    ark_mem->ark_msbp = msbp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetPredictorMethod:

 Specifies the method to use for predicting implicit solutions.  
 Non-default choices are {1,2,3,4}, all others will use default 
 (trivial) predictor.
---------------------------------------------------------------*/
int ARKodeSetPredictorMethod(void *arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetPredictorMethod", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetStabilityFn:

 Specifies the user-provided explicit time step stability 
 function to use.  A NULL input function implies a reset to
 the default function (empty).
---------------------------------------------------------------*/
int ARKodeSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab,
			 void *estab_data)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetStabilityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  if (EStab == NULL) {
    ark_mem->ark_expstab    = arkExpStab;
    ark_mem->ark_estab_data = ark_mem;
  } else {
    ark_mem->ark_expstab    = EStab;
    ark_mem->ark_estab_data = estab_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxErrTestFails:

 Specifies the maximum number of error test failures during one
 step try.  A non-positive input implies a reset to
 the default value.
---------------------------------------------------------------*/
int ARKodeSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxErrTestFails", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    ark_mem->ark_maxnef = MAXNEF;
  } else {
    ark_mem->ark_maxnef = maxnef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxConvFails:

 Specifies the maximum number of nonlinear convergence failures 
 during one step try.  A non-positive input implies a reset to
 the default value.
---------------------------------------------------------------*/
int ARKodeSetMaxConvFails(void *arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxConvFails", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    ark_mem->ark_maxncf = MAXNCF;
  } else {
    ark_mem->ark_maxncf = maxncf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxNonlinIters:

 Specifies the maximum number of nonlinear iterations during
 one solve.  A non-positive input implies a reset to the 
 default value.
---------------------------------------------------------------*/
int ARKodeSetMaxNonlinIters(void *arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxNonlinIters", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) {
    ark_mem->ark_maxcor = MAXCOR;
  } else {
    ark_mem->ark_maxcor = maxcor;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNonlinConvCoef:

 Specifies the coefficient in the nonlinear solver convergence
 test.  A non-positive input implies a reset to the default value.
---------------------------------------------------------------*/
int ARKodeSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNonlinConvCoef", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) {
    ark_mem->ark_nlscoef = NLSCOEF;
  } else {
    ark_mem->ark_nlscoef = nlscoef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetRootDirection:

 Specifies the direction of zero-crossings to be monitored.
 The default is to monitor both crossings.
---------------------------------------------------------------*/
int ARKodeSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  int i, nrt;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetRootDirection", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  nrt = ark_mem->ark_nrtfn;
  if (nrt==0) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetRootDirection", MSGARK_NO_ROOT);
    return(ARK_ILL_INPUT);    
  }

  for(i=0; i<nrt; i++) ark_mem->ark_rootdir[i] = rootdir[i];

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNoInactiveRootWarn:

 Disables issuing a warning if some root function appears
 to be identically zero at the beginning of the integration
---------------------------------------------------------------*/
int ARKodeSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNoInactiveRootWarn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_mxgnull = 0;
  
  return(ARK_SUCCESS);
}


/*===============================================================
 ARKODE optional output functions
===============================================================*/

/*---------------------------------------------------------------
 ARKodeGetNumSteps:

 Returns the current number of integration steps
---------------------------------------------------------------*/
int ARKodeGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumExpSteps:

 Returns the current number of stability-limited steps
---------------------------------------------------------------*/
int ARKodeGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumExpSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst_exp;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumAccSteps:

 Returns the current number of accuracy-limited steps
---------------------------------------------------------------*/
int ARKodeGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumAccSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumStepAttempts:

 Returns the current number of steps attempted by the solver
---------------------------------------------------------------*/
int ARKodeGetNumStepAttempts(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumStepAttempts", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst_attempts;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumRhsEvals:

 Returns the current number of calls to fe and fi
---------------------------------------------------------------*/
int ARKodeGetNumRhsEvals(void *arkode_mem, long int *fe_evals,
			 long int *fi_evals)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumRhsEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *fe_evals = ark_mem->ark_nfe;
  *fi_evals = ark_mem->ark_nfi;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumLinSolvSetups:

 Returns the current number of calls to the lsetup routine
---------------------------------------------------------------*/
int ARKodeGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumLinSolvSetups", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nlinsetups = ark_mem->ark_nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumMassSolves:

 Returns the current number of calls to the mass matrix solver.
---------------------------------------------------------------*/
int ARKodeGetNumMassSolves(void *arkode_mem, long int *nMassSolves)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumMassSolves", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nMassSolves = ark_mem->ark_mass_solves;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumMassMultiplies:

 Returns the current number of calls to the mass matrix product.
---------------------------------------------------------------*/
int ARKodeGetNumMassMultiplies(void *arkode_mem, long int *nMassMult)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumMassMult", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nMassMult = ark_mem->ark_mass_mult;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumErrTestFails:

 Returns the current number of error test failures
---------------------------------------------------------------*/
int ARKodeGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumErrTestFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *netfails = ark_mem->ark_netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetActualInitStep:

 Returns the step size used on the first step
---------------------------------------------------------------*/
int ARKodeGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetActualInitStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *hinused = ark_mem->ark_h0u;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetLastStep:

 Returns the step size used on the last successful step
---------------------------------------------------------------*/
int ARKodeGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetLastStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hlast = ark_mem->ark_hold;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentStep:

 Returns the step size to be attempted on the next step
---------------------------------------------------------------*/
int ARKodeGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetCurrentStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  
  *hcur = ark_mem->ark_next_h;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentTime:

 Returns the current value of the independent variable
---------------------------------------------------------------*/
int ARKodeGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetCurrentTime", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tcur = ark_mem->ark_tn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentButcherTables:

 Returns the explicit and implicit Butcher tables currently in 
 use.  The tables should be declared statically as 
 Ae[ARK_S_MAX*ARK_S_MAX] and Ai[ARK_S_MAX*ARK_S_MAX], and the 
 arrays c, b and b2 should all have length ARK_S_MAX.
---------------------------------------------------------------*/
int ARKodeGetCurrentButcherTables(void *arkode_mem, 
				  int *s, int *q, int *p,
				  realtype *Ai, realtype *Ae, 
				  realtype *c, realtype *b,
				  realtype *b2)
{
  int i,j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetCurrentButcherTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *s = ark_mem->ark_stages;
  *q = ark_mem->ark_q;
  *p = ark_mem->ark_p;
  for (i=0; i<ARK_S_MAX; i++) {
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(Ae,i,j) = ARK_A(ark_mem->ark_Ae,i,j);
      ARK_A(Ai,i,j) = ARK_A(ark_mem->ark_Ai,i,j);
    }
    c[i]  = ark_mem->ark_c[i];
    b[i]  = ark_mem->ark_b[i];
    b2[i] = ark_mem->ark_b2[i];
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetTolScaleFactor:

 Returns a suggested factor for scaling tolerances
---------------------------------------------------------------*/
int ARKodeGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetTolScaleFactor", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tolsfact = ark_mem->ark_tolsf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetErrWeights:

 This routine returns the current weight vector.
---------------------------------------------------------------*/
int ARKodeGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetErrWeights", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->ark_ewt, eweight);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetEstLocalErrors: (updated to the correct vector, but 
   need to verify that it is unchanged between filling the 
   estimated error and the end of the time step)

 Returns an estimate of the local error
---------------------------------------------------------------*/
int ARKodeGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetEstLocalErrors", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->ark_tempv, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetWorkSpace:

 Returns integrator work space requirements
---------------------------------------------------------------*/
int ARKodeGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetWorkSpace", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *leniw = ark_mem->ark_liw;
  *lenrw = ark_mem->ark_lrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetIntegratorStats:

 Returns integrator statistics
---------------------------------------------------------------*/
int ARKodeGetIntegratorStats(void *arkode_mem, long int *nsteps, 
			     long int *expsteps, long int *accsteps, 
			     long int *step_attempts, long int *fe_evals, 
			     long int *fi_evals, long int *nlinsetups, 
			     long int *netfails, realtype *hinused, 
			     realtype *hlast, realtype *hcur, 
			     realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetIntegratorStats", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps        = ark_mem->ark_nst;
  *expsteps      = ark_mem->ark_nst_exp;
  *accsteps      = ark_mem->ark_nst_acc;
  *step_attempts = ark_mem->ark_nst_attempts;
  *fe_evals      = ark_mem->ark_nfe;
  *fi_evals      = ark_mem->ark_nfi;
  *nlinsetups    = ark_mem->ark_nsetups;
  *netfails      = ark_mem->ark_netf;
  *hinused       = ark_mem->ark_h0u;
  *hlast         = ark_mem->ark_hold;
  *hcur          = ark_mem->ark_next_h;
  *tcur          = ark_mem->ark_tn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumGEvals:

 Returns the current number of calls to g (for rootfinding)
---------------------------------------------------------------*/
int ARKodeGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumGEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *ngevals = ark_mem->ark_nge;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetRootInfo:

 Returns pointer to array rootsfound showing roots found
---------------------------------------------------------------*/
int ARKodeGetRootInfo(void *arkode_mem, int *rootsfound)
{
  int i, nrt;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetRootInfo", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  nrt = ark_mem->ark_nrtfn;

  for (i=0; i<nrt; i++) rootsfound[i] = ark_mem->ark_iroots[i];

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumNonlinSolvIters:

 Returns the current number of nonlinear solver iterations 
---------------------------------------------------------------*/
int ARKodeGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumNonlinSolvIters", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nniters = ark_mem->ark_nni;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumNonlinSolvConvFails:

 Returns the current number of nonlinear solver convergence fails
---------------------------------------------------------------*/
int ARKodeGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumNonlinSolvConvFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nncfails = ark_mem->ark_ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNonlinSolvStats:

 Returns nonlinear solver statistics
---------------------------------------------------------------*/
int ARKodeGetNonlinSolvStats(void *arkode_mem, long int *nniters, 
			     long int *nncfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNonlinSolvStats", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nniters = ark_mem->ark_nni;
  *nncfails = ark_mem->ark_ncfn;

  return(ARK_SUCCESS);
}

/*-----------------------------------------------------------------*/

char *ARKodeGetReturnFlagName(long int flag)
{
  char *name;
  name = (char *)malloc(24*sizeof(char));

  switch(flag) {
  case ARK_SUCCESS:
    sprintf(name,"ARK_SUCCESS");
    break;
  case ARK_TSTOP_RETURN:
    sprintf(name,"ARK_TSTOP_RETURN");
    break;
  case ARK_ROOT_RETURN:
    sprintf(name,"ARK_ROOT_RETURN");
    break;
  case ARK_TOO_MUCH_WORK:
    sprintf(name,"ARK_TOO_MUCH_WORK");
    break;
  case ARK_TOO_MUCH_ACC:
    sprintf(name,"ARK_TOO_MUCH_ACC");
    break;
  case ARK_ERR_FAILURE:
    sprintf(name,"ARK_ERR_FAILURE");
    break;
  case ARK_CONV_FAILURE:
    sprintf(name,"ARK_CONV_FAILURE");
    break;
  case ARK_LINIT_FAIL:
    sprintf(name,"ARK_LINIT_FAIL");
    break;
  case ARK_LSETUP_FAIL:
    sprintf(name,"ARK_LSETUP_FAIL");
    break;
  case ARK_LSOLVE_FAIL:
    sprintf(name,"ARK_LSOLVE_FAIL");
    break;
  case ARK_RHSFUNC_FAIL:
    sprintf(name,"ARK_RHSFUNC_FAIL");
    break;
  case ARK_FIRST_RHSFUNC_ERR:
    sprintf(name,"ARK_FIRST_RHSFUNC_ERR");
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    sprintf(name,"ARK_REPTD_RHSFUNC_ERR");
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    sprintf(name,"ARK_UNREC_RHSFUNC_ERR");
    break;
  case ARK_RTFUNC_FAIL:
    sprintf(name,"ARK_RTFUNC_FAIL");
    break;
  case ARK_LFREE_FAIL:
    sprintf(name,"ARK_LFREE_FAIL");
    break;
  case ARK_MASSINIT_FAIL:
    sprintf(name,"ARK_MASSINIT_FAIL");
    break;
  case ARK_MASSSETUP_FAIL:
    sprintf(name,"ARK_MASSSETUP_FAIL");
    break;
  case ARK_MASSSOLVE_FAIL:
    sprintf(name,"ARK_MASSSOLVE_FAIL");
    break;
  case ARK_MASSFREE_FAIL:
    sprintf(name,"ARK_MASSFREE_FAIL");
    break;
  case ARK_MASSMULT_FAIL:
    sprintf(name,"ARK_MASSMULT_FAIL");
    break;
  case ARK_MEM_FAIL:
    sprintf(name,"ARK_MEM_FAIL");
    break;
  case ARK_MEM_NULL:
    sprintf(name,"ARK_MEM_NULL");
    break;
  case ARK_ILL_INPUT:
    sprintf(name,"ARK_ILL_INPUT");
    break;
  case ARK_NO_MALLOC:
    sprintf(name,"ARK_NO_MALLOC");
    break;
  case ARK_BAD_K:
    sprintf(name,"ARK_BAD_K");
    break;
  case ARK_BAD_T:
    sprintf(name,"ARK_BAD_T");
    break;
  case ARK_BAD_DKY:
    sprintf(name,"ARK_BAD_DKY");
    break;
  case ARK_TOO_CLOSE:
    sprintf(name,"ARK_TOO_CLOSE");
    break;    
  default:
    sprintf(name,"NONE");
  }

  return(name);
}



/*===============================================================
  ARKODE parameter output
===============================================================*/

/*---------------------------------------------------------------
 ARKodeWriteParameters:

 Outputs all solver parameters to the provided file pointer.
---------------------------------------------------------------*/
int ARKodeWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeWriteParameters", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* print integrator parameters to file */
  fprintf(fp, "ARKode solver parameters:\n");
  fprintf(fp, "  Method order %i\n",ark_mem->ark_q);
  fprintf(fp, "  Dense output order %i\n",ark_mem->ark_dense_q);
  if (ark_mem->ark_linear) {
    fprintf(fp, "  Linear implicit problem");
    if (ark_mem->ark_linear_timedep) {
      fprintf(fp, " (time-dependent Jacobian)\n");
    } else {
      fprintf(fp, " (time-independent Jacobian)\n");
    }
  }
  if (ark_mem->ark_explicit) {
    fprintf(fp, "  Explicit integrator\n");
  } else if (ark_mem->ark_implicit) {
    fprintf(fp, "  Implicit integrator\n");
  } else {
    fprintf(fp, "  ImEx integrator\n");
  }
  if (ark_mem->ark_fixedstep) {
    fprintf(fp, "  Fixed time-stepping enabled\n");
  } else {
    if (ark_mem->ark_hadapt == NULL) {
      fprintf(fp, "  Time step adaptivity method %i\n", ark_mem->ark_hadapt_imethod);
      fprintf(fp, "     Safety factor = %g\n", ark_mem->ark_hadapt_safety);
      fprintf(fp, "     Bias factor = %g\n", ark_mem->ark_hadapt_bias);
      fprintf(fp, "     Growth factor = %g\n", ark_mem->ark_hadapt_growth);
      fprintf(fp, "     Step growth lower bound = %g\n", ark_mem->ark_hadapt_lbound);
      fprintf(fp, "     Step growth upper bound = %g\n", ark_mem->ark_hadapt_ubound);
      fprintf(fp, "     k1 = %g\n", ark_mem->ark_hadapt_k1);
      fprintf(fp, "     k2 = %g\n", ark_mem->ark_hadapt_k2);
      fprintf(fp, "     k3 = %g\n", ark_mem->ark_hadapt_k3);
    } else {
      fprintf(fp, "  User provided time step adaptivity function\n");
    }
  }
  if (ark_mem->ark_itol == ARK_WF) {
    fprintf(fp, "  User provided error weight function\n");
  } else {
    fprintf(fp, "  Solver relative tolerance = %g\n", ark_mem->ark_reltol);
    if (ark_mem->ark_itol == ARK_SS) {
      fprintf(fp, "  Solver absolute tolerance = %g\n", ark_mem->ark_Sabstol);
    } else {
      fprintf(fp, "  Vector-valued solver absolute tolerance\n");
    }
  }
  if (!ark_mem->ark_rwt_is_ewt) {
    if (ark_mem->ark_ritol == ARK_WF) {
      fprintf(fp, "  User provided residual weight function\n");
    } else {
      if (ark_mem->ark_ritol == ARK_SS) {
	fprintf(fp, "  Absolute residual tolerance = %g\n", ark_mem->ark_SRabstol);
      } else {
	fprintf(fp, "  Vector-valued residual absolute tolerance\n");
      }
    }
  }
  if (ark_mem->ark_hin != ZERO)  
    fprintf(fp, "  Initial step size = %g\n",ark_mem->ark_hin);
  if (ark_mem->ark_hmin != ZERO)  
    fprintf(fp, "  Minimum step size = %g\n",ark_mem->ark_hmin);
  if (ark_mem->ark_hmax_inv != ZERO)  
    fprintf(fp, "  Maximum step size = %g\n",ONE/ark_mem->ark_hmax_inv);
  fprintf(fp, "  Maximum number of error test failures = %i\n",ark_mem->ark_maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",ark_mem->ark_maxncf);
  fprintf(fp, "  Maximum step increase (first step) = %g\n",ark_mem->ark_etamx1);
  fprintf(fp, "  Step reduction factor on multiple error fails = %g\n",ark_mem->ark_etamxf);
  fprintf(fp, "  Minimum error fails before above factor is used = %i\n",ark_mem->ark_small_nef);
  fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %g\n",ark_mem->ark_etacf);

  if (!ark_mem->ark_implicit) {
    if (ark_mem->ark_expstab == arkExpStab) {
      fprintf(fp, "  Default explicit stability function\n");
    } else {
      fprintf(fp, "  User provided explicit stability function\n");
    }
    fprintf(fp, "  Explicit safety factor = %g\n",ark_mem->ark_hadapt_cfl);
  }
  if (!ark_mem->ark_explicit) {
    fprintf(fp, "  Implicit predictor method = %i\n",ark_mem->ark_predictor);
    fprintf(fp, "  Implicit solver tolerance coefficient = %g\n",ark_mem->ark_nlscoef);
    fprintf(fp, "  Maximum number of nonlinear corrections = %i\n",ark_mem->ark_maxcor);
    fprintf(fp, "  Nonlinear convergence rate constant = %g\n",ark_mem->ark_crdown);
    fprintf(fp, "  Nonlinear divergence tolerance = %g\n",ark_mem->ark_rdiv);
    fprintf(fp, "  Gamma factor LSetup tolerance = %g\n",ark_mem->ark_dgmax);
    fprintf(fp, "  Number of steps between LSetup calls = %i\n",ark_mem->ark_msbp);
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeWriteButcher:

 Outputs Butcher tables to the provided file pointer.
---------------------------------------------------------------*/
int ARKodeWriteButcher(void *arkode_mem, FILE *fp)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeWriteButcher", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* print Butcher tables to file */
  fprintf(fp, "\nARKode Butcher tables (stages = %i):\n", ark_mem->ark_stages);
  if (!ark_mem->ark_implicit) {
    fprintf(fp, "  Explicit Butcher table:\n");
    for (i=0; i<ark_mem->ark_stages; i++) {
      fprintf(fp, "     %.5f",ark_mem->ark_c[i]);
      for (j=0; j<ark_mem->ark_stages; j++) 
	fprintf(fp, " %.5f",ARK_A(ark_mem->ark_Ae,i,j));
      fprintf(fp,"\n");
    }
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b[j]);
    fprintf(fp,"\n");
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b2[j]);
    fprintf(fp,"\n");
  }
  if (!ark_mem->ark_explicit) {
    fprintf(fp, "  Implicit Butcher table:\n");
    for (i=0; i<ark_mem->ark_stages; i++) {
      fprintf(fp, "     %.5f",ark_mem->ark_c[i]);
      for (j=0; j<ark_mem->ark_stages; j++) 
	fprintf(fp, " %.5f",ARK_A(ark_mem->ark_Ai,i,j));
      fprintf(fp,"\n");
    }
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b[j]);
    fprintf(fp,"\n");
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b2[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}




/*---------------------------------------------------------------
      EOF
---------------------------------------------------------------*/
