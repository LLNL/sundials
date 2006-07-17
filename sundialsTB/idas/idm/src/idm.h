/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-17 16:49:51 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/idas/LICENSE.
 * -----------------------------------------------------------------
 * Header file for the IDAS Matlab interface.
 * -----------------------------------------------------------------
 */

#ifndef _IDM_H
#define _IDM_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <ida/ida_band.h>
#include <ida/ida_spgmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_sptfqmr.h>
#include <ida/ida_bbdpre.h>

  /*
   * ---------------------------------------------------------------------------------
   * Constants
   * ---------------------------------------------------------------------------------
   */

  /* Linear solver types */

  enum {LS_DENSE, LS_BAND, LS_SPGMR, LS_SPBCG, LS_SPTFQMR};

  /* Preconditioner modules */

  enum {PM_NONE, PM_BBDPRE};

  /*
   * ---------------------------------------------------------------------------------
   * Types for global data structures
   * ---------------------------------------------------------------------------------
   */

  typedef struct idm_IDASdataStruct {

    void *ida_mem;     /* IDAS solver memory */
    void *bp_data;     /* Preconditioner memory (BBDPre) */
    N_Vector yy;       /* solution vector yy */
    N_Vector yp;       /* solution vector yp */
    N_Vector yQ;       /* quadratures vector */
    N_Vector *yyS;     /* sensitivity vectors yyS */
    N_Vector *ypS;     /* sensitivity vectors ypS */
    int N;             /* problem dimension */
    int Nq;            /* number of quadratures */
    int Ng;            /* number of root functions */
    int Ns;            /* number of sensitivities */
    int Nd;            /* number of data points */
    int Nc;            /* number of check points */
    int ls;            /* linear solver type */
    int pm;            /* preconditioner module */
    int icopt;         /* IC calculation */
    int ism;           /* sensitivity method */
    
    void *idaadj_mem;  /* IDAS adjoint memory */
    int interp;
    N_Vector yyB;
    N_Vector ypB;
    N_Vector yQB;
    int NB;
    int NqB;
    int lsB;
    int pmB;
    int icoptB;

  } *idm_IDASdata;

  typedef struct idm_MATLABdataStruct {

    mxArray *mx_RESfct;
    mxArray *mx_QUADfct;
    mxArray *mx_JACfct;
    mxArray *mx_PSETfct;
    mxArray *mx_PSOLfct;
    mxArray *mx_GLOCfct;
    mxArray *mx_GCOMfct;
    
    mxArray *mx_Gfct;
    
    mxArray *mx_SRESfct;
    
    mxArray *mx_RESfctB;
    mxArray *mx_QUADfctB;
    mxArray *mx_JACfctB;
    mxArray *mx_PSETfctB;
    mxArray *mx_PSOLfctB;
    mxArray *mx_GLOCfctB;
    mxArray *mx_GCOMfctB;

    mxArray *mx_data;

    /* Monitor */
  
    mxArray *mx_MONfct;
    mxArray *mx_MONdata;
    
    mxArray *mx_MONfctB;
    mxArray *mx_MONdataB;

  } *idm_MATLABdata;

  /*
   * ---------------------------------------------------------------------------------
   * Declarations for global variables (defined in idm.c)
   * ---------------------------------------------------------------------------------
   */

  extern idm_IDASdata idm_Cdata;    /* IDAS data */
  extern booleantype idm_quad;      /* Forward quadratures? */
  extern booleantype idm_quadB;     /* Backward quadratures? */
  extern booleantype idm_asa;       /* Adjoint sensitivity? */
  extern booleantype idm_fsa;       /* Forward sensitivity? */
  extern booleantype idm_mon;       /* Forward monitoring? */ 
  extern booleantype idm_monB;      /* Backward monitoring? */ 

  extern idm_MATLABdata idm_Mdata;  /* MATLAB data */

  /*
   * ---------------------------------------------------------------------------------
   * Wrapper functions
   * ---------------------------------------------------------------------------------
   */
  void mtlb_IdaErrHandler(int error_code, 
                          const char *module, const char *function, 
                          char *msg, void *eh_data); 
  
  int mtlb_IdaRes(realtype tt, N_Vector yy, N_Vector yp,
                  N_Vector rr, void *res_data);

  int mtlb_IdaGfct(realtype t, N_Vector y, N_Vector yp,
                   realtype *gout, void *g_data);

  int mtlb_IdaQuadFct(realtype tres, N_Vector yy, N_Vector yp, N_Vector ypQ,
                      void *rdataQ);

  int mtlb_IdaSensRes(int Nsens, realtype tres, 
                      N_Vector yy, N_Vector yp, N_Vector resval,
                      N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                      void *rdataS,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  int mtlb_IdaDenseJac(long int Neq, realtype tt, 
                       N_Vector yy, N_Vector yp, N_Vector rr,
                       realtype c_j, void *jac_data, 
                       DenseMat Jac, 
                       N_Vector tmp1, N_Vector tmp2, 
                       N_Vector tmp3);

  int mtlb_IdaBandJac(long int Neq, long int mupper, 
                      long int mlower, realtype tt, 
                      N_Vector yy, N_Vector yp, N_Vector rr, 
                      realtype c_j, void *jac_data, BandMat Jac, 
                      N_Vector tmp1, N_Vector tmp2, 
                      N_Vector tmp3);

  int mtlb_IdaSpilsJac(realtype tt,
                       N_Vector yy, N_Vector yp, N_Vector rr,
                       N_Vector v, N_Vector Jv,
                       realtype c_j, void *jac_data,
                       N_Vector tmp1, N_Vector tmp2);
  int mtlb_IdaSpilsPset(realtype tt,
                        N_Vector yy, N_Vector yp, N_Vector rr,
                        realtype c_j, void *prec_data,
                        N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);
  int mtlb_IdaSpilsPsol(realtype tt,
                        N_Vector yy, N_Vector yp, N_Vector rr,
                        N_Vector rvec, N_Vector zvec,
                        realtype c_j, realtype delta, void *prec_data,
                        N_Vector tmp);

  int mtlb_IdaBBDgloc(long int Nlocal, realtype tt,
                      N_Vector yy, N_Vector yp, N_Vector gval,
                      void *res_data);
  int mtlb_IdaBBDgcom(long int Nlocal, realtype tt,
                      N_Vector yy, N_Vector yp,
                      void *res_data);

  int mtlb_IdaResB(realtype tt, 
                   N_Vector yy, N_Vector yp,
                   N_Vector yyB, N_Vector ypB, N_Vector rrB,
                   void *rdataB);

  int mtlb_IdaQuadfctB(realtype tt, 
                       N_Vector yy, N_Vector yp, 
                       N_Vector yyB, N_Vector ypB,
                       N_Vector ypQB, void *rdataQB);

  int mtlb_IdaDenseJacB(long int NeqB, realtype tt, 
                        N_Vector yy, N_Vector yp,
                        N_Vector yyB, N_Vector ypB, N_Vector rrB,
                        realtype c_jB, void *jac_dataB, 
                        DenseMat JacB, 
                        N_Vector tmp1B, N_Vector tmp2B, 
                        N_Vector tmp3B);

  int mtlb_IdaBandJacB(long int NeqB, 
                       long int mupperB, long int mlowerB, 
                       realtype tt, 
                       N_Vector yy, N_Vector yp,
                       N_Vector yyB, N_Vector ypB, N_Vector rrB,
                       realtype c_jB, void *jac_dataB,
                       BandMat JacB, 
                       N_Vector tmp1B, N_Vector tmp2B, 
                       N_Vector tmp3B);

  int mtlb_IdaSpilsJacB(realtype t,
                        N_Vector yy, N_Vector yp,
                        N_Vector yyB, N_Vector ypB, N_Vector rrB,
                        N_Vector vB, N_Vector JvB, 
                        realtype c_jB, void *jac_dataB, 
                        N_Vector tmp1B, N_Vector tmp2B);
  int mtlb_IdaSpilsPsetB(realtype tt, 
                         N_Vector yy, N_Vector yp,
                         N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                         realtype c_jB, void *prec_dataB,
                         N_Vector tmp1B, N_Vector tmp2B, 
                         N_Vector tmp3B);
  int mtlb_IdaSpilsPsolB(realtype tt, 
                         N_Vector yy, N_Vector yp,
                         N_Vector yyB, N_Vector ypB, N_Vector rrB, 
                         N_Vector rvecB, N_Vector zvecB,
                         realtype c_jB, realtype deltaB,
                         void *prec_dataB, N_Vector tmpB);
  
  int mtlb_IdaBBDglocB(long int NlocalB, realtype tt,
                       N_Vector yy, N_Vector yp, 
                       N_Vector yyB, N_Vector ypB, N_Vector gvalB,
                       void *res_dataB);
  int mtlb_IdaBBDgcomB(long int NlocalB, realtype tt,
                       N_Vector yy, N_Vector yp,
                       N_Vector yyB, N_Vector ypB,
                       void *res_dataB);

  void mtlb_IdaMonitor(int call, double t, 
                       N_Vector yy, N_Vector yp, 
                       N_Vector yQ, 
                       N_Vector *yyS, N_Vector *ypS);
  void mtlb_IdaMonitorB(int call, double tB, 
                        N_Vector yyB, N_Vector ypB, 
                        N_Vector yQB);

  /*
   * ---------------------------------------------------------------------------------
   * Option handling functions
   * ---------------------------------------------------------------------------------
   */

  int get_IntgrOptions(const mxArray *options, booleantype fwd,
                       int *maxord,
                       long int *mxsteps,
                       int *itol, realtype *reltol, double *Sabstol, double **Vabstol,
                       double *hin, double *hmax,
                       double *tstop, booleantype *tstopSet,
                       double **id, double **cnstr);

  int get_LinSolvOptions(const mxArray *options, booleantype fwd,
                         int *mupper, int *mlower,
                         int *mudq, int *mldq,
                         int *gstype, int *maxl);

  int get_QuadOptions(const mxArray *options, booleantype fwd,
                      double **yQ0, booleantype *errconQ,
                      int *itolQ, double *reltolQ, double *SabstolQ, double **VabstolQ);

  int get_FSAOptions(const mxArray *options, 
                     char **pfield_name, int **plist, double **pbar,
                     booleantype *userSRES, double *rho,
                     booleantype *errconS, int *itolS, double *reltolS, 
                     double **SabstolS, double **VabstolS);

#ifdef __cplusplus
}
#endif

#endif
