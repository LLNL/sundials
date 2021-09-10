/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2012-03-07 21:50:32 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Header file for the KINSOL Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _KIM_H
#define _KIM_H

#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_bbdpre.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ---------------------------------------------------------------------------------
 * Constants
 * ---------------------------------------------------------------------------------
 */


/* Linear solver types */

enum {LS_NONE, LS_DENSE, LS_BAND, LS_SPGMR, LS_SPBCG, LS_SPTFQMR};

/* Preconditioner modules */

enum {PM_NONE, PM_BBDPRE};

/*
 * ---------------------------------------------------------------------------------
 * Types for global data structures
 * ---------------------------------------------------------------------------------
 */

typedef struct kimInterfaceData_ {

  void *kin_mem;        /* KINSOL solver memory */

  sunindextype n;           /* problem dimension */

  N_Vector Y;           /* solution vector */

  int LS;               /* linear solver type */
  int PM;               /* preconditioner module */

  booleantype errMsg;   /* post error/warning messages? */

  int fig_handle;       /* figure for posting info */

  /* Matlab functions and data associated with this problem */

  mxArray *SYSfct;

  mxArray *JACfct;

  mxArray *PSETfct;
  mxArray *PSOLfct;

  mxArray *GLOCfct;
  mxArray *GCOMfct;

  mxArray *mtlb_data;

} *kimInterfaceData;

/*
 * ---------------------------------------------------------------------------------
 * Error and info handler functions
 * ---------------------------------------------------------------------------------
 */

void kimErrHandler(int error_code, 
                   const char *module, const char *function, 
                   char *msg, void *eh_data); 

void kimInfoHandler(const char *module, const char *function, 
                    char *msg, void *ih_data); 

/*
 * ---------------------------------------------------------------------------------
 * Wrapper functions
 * ---------------------------------------------------------------------------------
 */


int mxW_KINSys(N_Vector y, N_Vector fy, void *user_data );

/* Dense direct linear solver */

int mxW_KINDenseJac(sunindextype N,
                    N_Vector y, N_Vector fy, 
                    DlsMat J, void *user_data,
                    N_Vector tmp1, N_Vector tmp2);

/* Band direct linear solver */

int mxW_KINBandJac(sunindextype N, sunindextype mupper, sunindextype mlower,
                   N_Vector u, N_Vector fu, 
                   DlsMat J, void *user_data,
                   N_Vector tmp1, N_Vector tmp2);

/* Scaled Preconditioned Iterative Linear Solver (SPGMR or SPBCG) */

int mxW_KINSpilsJac(N_Vector v, N_Vector Jv,
                    N_Vector y, booleantype *new_y, 
                    void *user_data);
int mxW_KINSpilsPset(N_Vector y, N_Vector yscale,
                     N_Vector fy, N_Vector fscale,
                     void *user_data, N_Vector vtemp1,
                     N_Vector vtemp2);
int mxW_KINSpilsPsol(N_Vector y, N_Vector yscale, 
                     N_Vector fy, N_Vector fscale, 
                     N_Vector v, void *user_data,
                     N_Vector vtemp);

/* BBD Preconditioner */

int mxW_KINGloc(sunindextype Nlocal, N_Vector y, N_Vector gval, void *user_data);
int mxW_KINGcom(sunindextype Nlocal, N_Vector y, void *user_data);

/*
 * ---------------------------------------------------------------------------------
 * Option handling functions
 * ---------------------------------------------------------------------------------
 */

int get_SolverOptions(const mxArray *options,
                      booleantype *verbose, booleantype *errmsg,
                      int *mxiter, int *msbset, int *msbsetsub,
                      int *etachoice, int *mxnbcf,
                      double *eta, double *egamma, double *ealpha, double *mxnewtstep, 
                      double *relfunc, double *fnormtol, double *scsteptol,
                      double **constraints,
                      booleantype *noInitSetup, booleantype *noMinEps);

int get_LinSolvOptions(const mxArray *options,
                       sunindextype *mupper, sunindextype *mlower,
                       sunindextype *mudq, sunindextype *mldq, double *dqrely,
                       int *ptype, int *maxrs, int *maxl);

#ifdef __cplusplus
}
#endif

#endif
