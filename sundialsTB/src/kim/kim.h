/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006-02-02 00:39:06 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/kinsol/LICENSE.
 * -----------------------------------------------------------------
 * Header file for the KINSOL Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _KIM_H
#define _KIM_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "kinsol.h"
#include "kinsol_dense.h"
#include "kinsol_spgmr.h"
#include "kinsol_spbcgs.h"
#include "kinsol_bbdpre.h"
#include "sundials_nvector.h"

  /*
   * ---------------------------------------------------------------------------------
   * Constants
   * ---------------------------------------------------------------------------------
   */


  /* Linear solver types */

  enum {LS_NONE, LS_DENSE, LS_SPGMR, LS_SPBCG, LS_SPTFQMR};

  /* Preconditioner modules */

  enum {PM_NONE, PM_BBDPRE};

  /*
   * ---------------------------------------------------------------------------------
   * Declarations for global variables
   * ---------------------------------------------------------------------------------
   */

  /* KINSOL data */

  extern void *kin_mem;     /* KINSOL solver memory */
  extern void *bbd_data;    /* BBD preconditioner data */
  extern N_Vector y;        /* solution vector */
  extern int N;             /* problem dimension */
  extern int ls;            /* linear solver type */
  extern int pm;            /* preconditioner module */

  /* Matlab data */

  extern mxArray *mx_mtlb_SYSfct;
  extern mxArray *mx_mtlb_JACfct;
  extern mxArray *mx_mtlb_PSETfct;
  extern mxArray *mx_mtlb_PSOLfct;
  extern mxArray *mx_mtlb_GLOCfct;
  extern mxArray *mx_mtlb_GCOMfct;

  extern mxArray *mx_mtlb_data;

  extern int fig_handle;
  
  /*
   * ---------------------------------------------------------------------------------
   * Wrapper functions
   * ---------------------------------------------------------------------------------
   */

  void mtlb_KINErrHandler(int error_code, 
                          const char *module, const char *function, 
                          char *msg, void *eh_data); 
  
  void mtlb_KINInfoHandler(const char *module, const char *function, 
                           char *msg, void *ih_data); 

  int mtlb_KINSys(N_Vector y, N_Vector fy, void *f_data );

  /* Dense direct linear solver */

  int mtlb_KINDenseJac(long int N, DenseMat J, 
                       N_Vector y, N_Vector fy, void *jac_data,
                       N_Vector tmp1, N_Vector tmp2);

  /* Scaled Preconditioned Iterative Linear Solver (SPGMR or SPBCG) */

  int mtlb_KINSpilsJac(N_Vector v, N_Vector Jv,
                       N_Vector y, booleantype *new_y, 
                       void *J_data);
  int mtlb_KINSpilsPset(N_Vector y, N_Vector yscale,
                        N_Vector fy, N_Vector fscale,
                        void *P_data, N_Vector vtemp1,
                        N_Vector vtemp2);
  int mtlb_KINSpilsPsol(N_Vector y, N_Vector yscale, 
                        N_Vector fy, N_Vector fscale, 
                        N_Vector v, void *P_data,
                        N_Vector vtemp);

  /* BBD Preconditioner */
  
  int mtlb_KINGloc(long int Nlocal, N_Vector y, N_Vector gval, void *f_data);
  int mtlb_KINGcom(long int Nlocal, N_Vector y, void *f_data);

  /*
   * ---------------------------------------------------------------------------------
   * Option handling functions
   * ---------------------------------------------------------------------------------
   */

  int get_SolverOptions(const mxArray *options,
                        booleantype *verbose,
                        int *mxiter, int *msbset, int *etachoice, int *mxnbcf,
                        double *eta, double *egamma, double *ealpha, double *mxnewtstep, 
                        double *relfunc, double *fnormtol, double *scsteptol,
                        double **constraints,
                        booleantype *noInitSetup, booleantype *noMinEps);

  int get_LinSolvOptions(const mxArray *options,
                         int *ls_tmp,
                         int *maxl, int *maxrs, int *ptype, int *pm_tmp,
                         int *mudq, int *mldq, int *mukeep, int *mlkeep, 
                         double *dq,
                         mxArray **mx_tmp_JACfct,
                         mxArray **mx_tmp_PSETfct, mxArray **mx_tmp_PSOLfct,
                         mxArray **mx_tmp_GLOCfct, mxArray **mx_tmp_GCOMfct);
  
#ifdef __cplusplus
}
#endif

#endif
