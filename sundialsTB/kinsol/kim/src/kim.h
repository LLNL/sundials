/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 16:00:44 $
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
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_bbdpre.h>

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

  /* KINSOL data */

  typedef struct kim_KINSOLdataStruct {

    void *kin_mem;     /* KINSOL solver memory */
    void *bbd_data;    /* BBD preconditioner data */
    N_Vector y;        /* solution vector */
    int N;             /* problem dimension */
    int ls;            /* linear solver type */
    int pm;            /* preconditioner module */

  } *kim_KINSOLdata;

  /* Matlab data */

  typedef struct kim_MATLABdataStruct {

    mxArray *mx_SYSfct;
    mxArray *mx_JACfct;
    mxArray *mx_PSETfct;
    mxArray *mx_PSOLfct;
    mxArray *mx_GLOCfct;
    mxArray *mx_GCOMfct;
    mxArray *mx_data;
    int fig_handle;

  } *kim_MATLABdata;


  /*
   * ---------------------------------------------------------------------------------
   * Declarations for global variables (defined in kim.c)
   * ---------------------------------------------------------------------------------
   */

  extern kim_KINSOLdata kim_Kdata;  /* KINSOL data */
  extern kim_MATLABdata kim_Mdata;  /* MATLAB data */

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

  /* Band direct linear solver */

  int mtlb_KINBandJac(long int N, long int mupper, long int mlower,
                      BandMat J, N_Vector u, N_Vector fu, void *jac_data,
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
                        int *mxiter, int *msbset, int *msbsetsub,
                        int *etachoice, int *mxnbcf,
                        double *eta, double *egamma, double *ealpha, double *mxnewtstep, 
                        double *relfunc, double *fnormtol, double *scsteptol,
                        double **constraints,
                        booleantype *noInitSetup, booleantype *noMinEps);

  int get_LinSolvOptions(const mxArray *options,
                         int *mupper, int *mlower,
                         int *mudq, int *mldq, double *dqrely,
                         int *ptype, int *maxrs, int *maxl);

#ifdef __cplusplus
}
#endif

#endif
