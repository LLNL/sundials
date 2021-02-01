/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Michael Wittman, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This file contains implementations of routines for an
 * interface to hypre's BoomerAMG preconditioner, for use with CVODE, a CVSPILS linear
 * solver, and the ParHyp implementation of NVECTOR.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_hypamgpre_impl.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE.h"
#include "_hypre_parcsr_mv.h"
#include <nvector/nvector_parhyp.h>
#include "cvode_spils_impl.h"

#include <cvode/cvode_sptfqmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_spgmr.h>

#include <sundials/sundials_math.h>
#include <math.h>

#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Prototypes of functions CVBoomerAMGSetup and CVBoomerAMGSolve */


static int CVBoomerAMGSetup(realtype t, N_Vector y, N_Vector fy, 
                          booleantype jok, booleantype *jcurPtr, 
                          realtype gamma, void *hypamg_data, 
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int CVBoomerAMGSolve(realtype t, N_Vector y, N_Vector fy, 
                          N_Vector r, N_Vector z, 
                          realtype gamma, realtype delta,
                          int lr, void *hypamg_data, N_Vector tmp);

/* Prototype for CVBoomerAMGFree */
static void CVBoomerAMGFree(CVodeMem cv_mem);
                      
static int CVParCsrCreateGammaIdentity(HYPRE_IJMatrix* id_mat, realtype gamma, sunindextype ilower, sunindextype iupper);

/*
 * -----------------------------------------------------------------
 * User-Callable Functions: initialization, reinit and free
 * -----------------------------------------------------------------
 */

int CVBoomerAMGInit(void *cvode_mem, int ilower, int iupper, int jlower, int jupper, int N)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBoomerAMGData pdata;
  sunindextype muk, mlk, storage_mu;
  int flag;
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVHYPRE_BOOMERAMG", "CVBoomerAMGInit", MSGBBD_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;
  /* Test if a linear solver has been attached */
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVHYPRE_BOOMERAMG", "CVBoomerAMGInit", MSGBBD_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  /* Test if the NVECTOR package is compatible with the preconditioner */
  if(NV_HYPRE_PARVEC_PH((cv_mem->cv_tempv)) == NULL) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVHYPRE_BOOMERAMG", "CVBoomerAMGInit", MSGBBD_BAD_NVECTOR);
    return(CVSPILS_ILL_INPUT);
  }
  /* Allocate data memory */
  pdata = NULL;
  pdata = (CVBoomerAMGData) malloc(sizeof *pdata);  
  if (pdata == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVHYPRE_BOOMERAMG", "CVBoomerAMGInit", MSGBBD_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }
  /*  */
  pdata->cvode_mem = cvode_mem;

  /* Set up solver structures */
 /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
  /* set up proper communicator from cvode_mem */
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, jlower, jupper, &(pdata->A));
   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(pdata->A, HYPRE_PARCSR);
   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(pdata->A);
/*   pdata->jacfn=jacobian_update;
   pdata->jacobian_data=jacobian_data;*/
   pdata->ilower=ilower;
   pdata->iupper=iupper;
   pdata->jlower=jlower;
   pdata->jupper=jupper;
   pdata->N     =N;
   
  /* Create solver */
//  HYPRE_BoomerAMGCreate(&(pdata->solver));

      /* Now set up the AMG preconditioner and specify any parameters */
  HYPRE_BoomerAMGCreate(&(pdata->solver));
  HYPRE_BoomerAMGSetPrintLevel((pdata->solver), 0); /* print amg solution info */
  HYPRE_BoomerAMGSetCoarsenType((pdata->solver), 6);
  HYPRE_BoomerAMGSetRelaxType((pdata->solver), 6); /* Sym G.S./Jacobi hybrid */
  HYPRE_BoomerAMGSetNumSweeps((pdata->solver), 1);
  HYPRE_BoomerAMGSetTol((pdata->solver), 0.0); /* conv. tolerance zero */
  HYPRE_BoomerAMGSetMaxIter((pdata->solver), 1); /* do only one iteration! */
//  HYPRE_BoomerAMGSetMaxRowSum((pdata->solver), .9);
   
   /*Calls for using as a solver*/
  /* Set some parameters (See Reference Manual for more parameters) */
//  HYPRE_BoomerAMGSetPrintLevel((pdata->solver), );  /* print solve info + parameters */
//  HYPRE_BoomerAMGSetCoarsenType(pdata->solver, 6); /* Falgout coarsening */
//  HYPRE_BoomerAMGSetRelaxType(pdata->solver, 3);   /* G-S/Jacobi hybrid relaxation */
//  HYPRE_BoomerAMGSetNumSweeps(pdata->solver, 1);   /* Sweeeps on each level */
//  HYPRE_BoomerAMGSetMaxLevels(pdata->solver, 20);  /* maximum number of levels */
//  HYPRE_BoomerAMGSetTol(pdata->solver, 1e-7);      /* conv. tolerance */   
  /* Overwrite the P_data field in the SPILS memory */
  cvspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  cvspils_mem->s_pfree = CVBoomerAMGFree;

  /* Attach preconditioner solve and setup functions */
  flag = CVSpilsSetPreconditioner(cvode_mem, CVBoomerAMGSetup, CVBoomerAMGSolve);
  return(flag);
}

/*
 * -----------------------------------------------------------------
 * Function : CVBoomerAMGSetup                                      
 * -----------------------------------------------------------------
 *
 * CVBoomerAMGSetup calculates a new J,if necessary, then calculates
 * P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of CVBoomerAMGSetup used here are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *         namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == SUNFALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == SUNTRUE  means that Jacobian data from the
 *                  previous CVBoomerAMGon call can be reused
 *                  (with the current value of gamma).
 *         A CVBoomerAMG call with jok == SUNTRUE should only occur
 *         after a call with jok == SUNFALSE.
 *         Currently this argument is superflous, since we cannot store a previous J
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         set by CVBoomerAMGon as follows:
 *           *jcurPtr = SUNTRUE if Jacobian data was recomputed.
 *           *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
 *                      but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * hypamg_data is a pointer to the preconditioner data set by
 *          CVBoomerAMGInit
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *           for NVectors which are be used by CVBoomerAMGSetup
 *           as temporary storage or work space.
 *
 * Return value:
 * The value returned by this CVBoomerAMGSetup function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried).
 * -----------------------------------------------------------------
 */

int CVBoomerAMGSetup(realtype t, N_Vector y, N_Vector fy, 
                     booleantype jok, booleantype *jcurPtr, 
                     realtype gamma, void *hypamg_data, 
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int counter=0;
  int nnz, i;
  double values[5];
  int cols[5];
//printf("Entered setup\n");
  long int ier;
  CVBoomerAMGData pdata;
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  int retval;

  pdata = (CVBoomerAMGData) hypamg_data;

  cv_mem = (CVodeMem) pdata->cvode_mem;

  *jcurPtr = SUNTRUE;
  
  /*Consider recomputing P in some Solve iterations if BoomerAMGSetup is much more
  expensive than computing P*/
  pdata->jacfn(pdata->N, pdata->ilower, pdata->iupper, pdata->jlower, pdata->jupper, gamma, t, y, fy, &(pdata->A), (void*) (cv_mem->cv_user_data), tmp1, tmp2, tmp3);
  HYPRE_IJMatrixGetObject(pdata->A, (void**) &(pdata->parcsr_A));
  
  pdata->par_b=(HYPRE_ParVector) (NV_HYPRE_PARVEC_PH(y));
  pdata->par_x=(HYPRE_ParVector) (NV_HYPRE_PARVEC_PH(y));
  HYPRE_BoomerAMGSetup(pdata->solver, pdata->parcsr_A, pdata->par_b, pdata->par_x);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVBoomerAMGSolve
 * -----------------------------------------------------------------
 * CVBoomerAMGSolve solves a linear system P z = r, with the
 * HYPRE_ParCSRMatrix generated by CVBoomerAMGSetup.
 *
 * The parameters of CVBoomerAMGSolve used here are as follows:
 *
 * r is the right-hand side vector of the linear system.
 *
 * hypamg_data is a pointer to the preconditioner data set by
 *   CVBoomerAMGInit.
 *
 * z is the output vector computed by CVBoomerAMGSolve.
 *
 * The value returned by the CVBoomerAMGSolve function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */

int CVBoomerAMGSolve(realtype t, N_Vector y, N_Vector fy, 
                          N_Vector r, N_Vector z, 
                          realtype gamma, realtype delta,
                          int lr, void *hypamg_data, N_Vector tmp)
{
  CVBoomerAMGData pdata;
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  
  pdata = (CVBoomerAMGData) hypamg_data;
  cv_mem = (CVodeMem) pdata->cvode_mem;

/*  jacfn(N, ilower, iupper, jlower, jupper, gamma, t, y, fy, &A, (void*) jacobian_data, tmp, tmp, tmp);*/
  pdata->par_b=(HYPRE_ParVector) (NV_HYPRE_PARVEC_PH(r));
  pdata->par_x=(HYPRE_ParVector) (NV_HYPRE_PARVEC_PH(z));
  HYPRE_BoomerAMGSolve(pdata->solver, pdata->parcsr_A, pdata->par_b, pdata->par_x);

  return(0);
}


static void CVBoomerAMGFree(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  CVBoomerAMGData pdata;
   int num_iterations;
   double final_res_norm;
   int nprocs, myid;
   int global_size;
  
  if (cv_mem->cv_lmem == NULL) return;
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  
  if (cvspils_mem->s_P_data == NULL) return;
  pdata = (CVBoomerAMGData) cvspils_mem->s_P_data;

  /* Need to return solver function so that users can get diagnostics out before CVodeFree()
  Run info - needed logging turned on */
/*  MPI_Comm_rank(hypre_ParCSRMatrixComm(parcsr_A),&myid);
  HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
  HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
  if (myid == 0)
  {
     printf("\n");
     printf("Iterations = %d\n", num_iterations);
     printf("Final Relative Residual Norm = %e\n", final_res_norm);
     printf("\n");
  }*/

  /* Destroy solver */
  HYPRE_BoomerAMGDestroy(pdata->solver);

  free(pdata);
  pdata = NULL;
}

/*
int JacTimes(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBoomerAMGData pdata;

  if (user_data->ode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) user_data->ode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  pdata = (CVBoomerAMGData) cvspils_mem->s_P_data;
 hypre_ParCSRMatrixMatvec ( 1.0, pdata->parcsr_A, NV_HYPRE_PARVEC_PH(v), 1.0, NV_HYPRE_PARVEC_PH(Jv));
}*/
/*
 * -----------------------------------------------------------------
 * CVParCsrSetParCsrJacFn
 * -----------------------------------------------------------------
 */

/*CVSpilsSetParCsrPCFn*/
int CVSpilsSetParCsrJacFn(void *cvode_mem, CVParCsrJacFn jparcsr)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBoomerAMGData pdata;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  pdata = (CVBoomerAMGData) cvspils_mem->s_P_data;

  if (jparcsr != NULL) {
    pdata->jacfn  = jparcsr;
/*    jtimesDQ = SUNFALSE;*/
  }/* else {
    jtimesDQ = SUNTRUE;
  }*/

  return(CVSPILS_SUCCESS);
}

int CVParCsrSetSpilsJacTimesVecFn(void *cvode_mem, CVSpilsJacTimesVecFn jparcsr)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBoomerAMGData pdata;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_CVMEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVSPILS", "CVSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  pdata = (CVBoomerAMGData) cvspils_mem->s_P_data;

  if (jparcsr != NULL) {
    cvspils_mem->s_jtimes  = jparcsr;
    cvspils_mem->s_jtimesDQ = SUNFALSE;
  }/* else {
    jtimesDQ = SUNTRUE;
  }*/

  return(CVSPILS_SUCCESS);
}

int CVParCsrCreateGammaIdentity(HYPRE_IJMatrix* id_mat, realtype gamma, sunindextype ilow, sunindextype iup)
{
  int nnz, i;
  double values[5];
  int cols[5];
  for (i = ilow; i <= iup; i++)
      {
         nnz = 0;
            
         cols[nnz] = i;
         values[nnz] = gamma;
         nnz++;
         HYPRE_IJMatrixSetValues(*id_mat, 1, &nnz, &i, cols, values);
      }
   

   /* Assemble after setting the coefficients */
     HYPRE_IJMatrixAssemble(*id_mat);
}
