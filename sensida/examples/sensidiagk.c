/***********************************************************************
 * File:       sensidiagk.c   
 * Written by: Steven L. Lee and Alan C. Hindmarsh
 * Version of: 17 September 2001
 *----------------------------------------------------------------------
 *
 * Sensitivity analysis version of example problem for IDA.
 * Diagonal test problem:
 *
 * y' + alpha*D*y = 0, alpha = 0.5, Nloc = # of ODEs per processor,
 *                     npes = # of processors, Ny = Nloc*npes, and
 *                     D = diag(1, 2, ..., Ny).
 *
 * y(0) = 1 for all components.
 *
 * The system is solved with SensIDA using the Krylov linear solver IDASPGMR. 
 * The preconditioner uses the diagonal elements of the Jacobian only.
 * Routines for preconditioning, required by IDASPGMR, are supplied here.
 * Output is taken at t = 0, 0.1, 0.2, ..., 1.0.
 *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"
#include "llnlmath.h"
#include "nvector.h"
#include "ida.h"
#include "idaspgmr.h"
#include "iterativ.h"
#include "mpi.h"
#include "sensida.h"

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define NOUT        10              /* Number of output times    */
#define Nloc         2              /* number of ODEs per processor */
#define Np           1              /* number of parameters      */
#define Ns           1              /* number of sensitivities   */

typedef struct {  
  real      *p;
  integer   Ny, thispe;
  /* vector of diagonal preconditioner elements */
  N_Vector  pp;    
  MPI_Comm  comm;
} *UserData;


/* User-supplied residual function */
int resdiag(integer Ny, real tres, N_Vector uu, N_Vector up,
            N_Vector res, void *rdata);

/* User-supplied preconditioner routines */
int Precond(integer local_N, real tt, N_Vector uu,
	    N_Vector up, N_Vector rr, real cj, ResFn res, void *rdata,
	    void *pdata, N_Vector ewt, real delta, N_Vector rvec,
	    N_Vector zvec, long int *nrePtr, N_Vector tempv);

int PrecondSetup(integer local_N, real tt, N_Vector yy,
		 N_Vector yp, N_Vector rr, real cj,
		 ResFn res, void *rdata, void *pdata,
		 N_Vector ewt, N_Vector constraints, real hh, 
		 real uround, long int *nrePtr,
		 N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

main(int argc, char *argv[])

{
  integer retval, i, j, iout, itol, itask, local_N, npes, thispe;
  integer Ny;
  long int iopt[OPT_SIZE];
  boole optIn;
  real ropt[OPT_SIZE], rtol, atol;
  real t0, t1, tout, tret, umax;
  void *mem;
  UserData data;
  N_Vector uu, up, res, id;
  N_Vector usoln, error;
  N_Vector *usolnsub, *errorsub;
  IDAMem idamem;
  MPI_Comm comm;
  machEnvType machEnv;
  N_Vector *uusub, *upsub, *idsub;
  integer Ntotal;
  real *pbar;
  real rhomax;
  real maxnorm;
  int *plist;

  /* Get processor number and total number of pe's. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &thispe);

  Ny = npes*Nloc;
  Ntotal = (1+Ns)*Ny;
  
  /* Set local length local_N. */
  local_N = Nloc;

  /* Set machEnv block. */
  machEnv = PVecInitMPI(comm, local_N, Ny, &argc, &argv);
  if (machEnv == NULL) return(1);

  /* Allocate and initialize the data structure and N-vectors. */
  data = (UserData) malloc(sizeof *data);

  /* Store nominal parameter values in p */
  data->p = (real *) malloc(Np * sizeof(real));
  data->p[0] = 0.5;

  /* Scaling factor for each sensitivity equation */
  pbar = (real *) malloc(Np * sizeof(real));
  pbar[0] = 1.0;

  /* Store ordering of parameters in plist */
  plist = (int *) malloc(Np * sizeof(int));
  plist[0] = 1;

  rhomax = 0.0;

  uu = N_VNew(Ntotal, machEnv); 
  up = N_VNew(Ntotal, machEnv);
  res = N_VNew(Ntotal, machEnv);
  id = N_VNew(Ntotal, machEnv);
  usoln = N_VNew(Ntotal, machEnv);
  error = N_VNew(Ntotal, machEnv);

  /* An N-vector to hold preconditioner. */
  data->pp = N_VNew(Ny, machEnv); 
  
  /* Create pointers to subvectors */
  uusub = N_VSUB(uu);
  upsub = N_VSUB(up);
  idsub = N_VSUB(id);
  usolnsub = N_VSUB(usoln);
  errorsub = N_VSUB(error);

  data->Ny = Ny;
  data->thispe = thispe;
  data->comm = comm;

  /**********************************************/
  /* Initialize uu, up; call IDACalcIC for up.  */
  /**********************************************/
  N_VConst(ONE, uu);
  N_VConst(ZERO, up);
  N_VConst(ONE, idsub[0]);

  /* Initialize the sensitivity variables */
  SensInitZero(uu, Ns);
  SensInitZero(up, Ns);

  /* Identify all sensitivity variables as differential variables */
  SensSetId(id, Ns);

  t0 = 0.0; t1 = 0.10;

  /* Scalar relative and absolute tolerance. */
  itol = SS;
  rtol = 0.0;
  atol = 1.e-5;

  optIn = TRUE;
  for (i = 0; i < OPT_SIZE; i++) {iopt[i] = 0; ropt[i] = ZERO; }

  /* Call SensIDAMalloc to initialize solution.  (NULL argument is errfp.) */
  itask = NORMAL;

  mem = SensIDAMalloc(Ny, Ns, Ntotal, resdiag, data, t0, uu, up, 
		      itol, &rtol, &atol, id, NULL, NULL, optIn, 
		      iopt, ropt, machEnv, data->p, pbar, plist, rhomax);
  if (mem == NULL) {
    if (thispe == 0) printf ("SensIDAMalloc failed.");
    return(1); }
  idamem = (IDAMem) mem;

  /* Call SensIDASpgmr to specify the linear solver. */
  retval = SensIDASpgmr(idamem, PrecondSetup, Precond, MODIFIED_GS, 
                    0, 0, 0.0, 0.0, data);

  if (retval != SUCCESS) {
    if (thispe == 0) printf("SensIDASpgmr failed, returning %d.\n",retval);
    return(1);
  }

  /* Call IDACalcIC (with default options) to correct the initial values. */
  retval = IDACalcIC(idamem, CALC_YA_YDP_INIT, t1, ZERO, 0,0,0,0, ZERO);
  if (retval != SUCCESS) {
    if (thispe == 0) printf("IDACalcIC failed. retval = %d\n", retval);
    return(1);
  }

  if (thispe == 0) { 
    printf("\n");
    printf("sensdiagk: Parallel example for SensIDA \n");
    printf("           Test problem: y' + alpha*D*y = 0. \n");
    printf("           Total system size: Ny = %d\n", Ny);
    printf("           D = diag(1,2,...,Ny). \n");
    printf("           alpha = %f\n", data->p[0]);
    printf("Number of sensitivities: Ns = %d\n", Ns);
    printf("Parameter value:      alpha = %9.2e\n", data->p[0]);
    printf("Scale factor:        pbar_1 = %9.2e\n", pbar[0]);
    printf("Finite difference:   rhomax = %9.2e\n", rhomax);
    printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
    printf("Constraints set to force all components of solution u >= 0. \n");
    printf("Linear solver: IDASPGMR  ");
    printf("Preconditioner: diagonal elements only.\n"); 

    /* Print output table heading and initial line of table. */
    printf("\n");
    printf("Output Summary:   err(u) = max-norm of error in solution u \n");
    printf("                err(s_i) = max-norm of error in sensitivity wrt p_i\n\n");
    printf("  time     k   nst  nni  nli   nre    h       npe   nps\n" );
    printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");

    printf(" %5.2f     %d  %3d  %3d  %3d  %4d %9.2e  %3d   %3d\n",
           t0, iopt[KUSED], iopt[NST], iopt[NNI], iopt[SPGMR_NLI], 
           iopt[NRE], ropt[HUSED], iopt[SPGMR_NPE], iopt[SPGMR_NPS]);
  }

    for (i = 0; i < Nloc; i++){
     usoln->data[i] = exp(-1.0*((thispe*Nloc)+i+1)*data->p[0]*t0);
     if (Ns > 0)
       usolnsub[1]->data[i]=-1.0*pbar[0]*((thispe*Nloc)+i+1)*t0*usoln->data[i]; 
    }
    N_VLinearSum(ONE, usoln, -ONE, uu, error);
    for (i = 0; i <= Ns; i++) {
      maxnorm = N_VMaxNorm(errorsub[i]);
      j = (plist == NULL) ? i : plist[i-1];
      if (thispe == 0) {
	if (i == 0) 
	  printf("err(u)   = %10.5e\n", maxnorm);
	else
	  printf("err(s_%d) = %10.5e\n", j, maxnorm);
      }
    }
    if (thispe == 0) printf("\n");

  /* Loop over tout, call IDASolve, print output. */
  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout += t1) {
 
    retval = IDASolve(idamem, tout, t0, &tret, uu, up, itask);

    if (retval < 0) {
      if (thispe == 0) printf("IDASolve returned %d.\n",retval);
      return(1);
    }

    if (thispe == 0) {
      printf(" %5.2f     %d  %3d  %3d  %3d  %4d %9.2e  %3d   %3d\n",
	     tret, iopt[KUSED], iopt[NST], iopt[NNI], iopt[SPGMR_NLI], 
	     iopt[NRE], ropt[HUSED], iopt[SPGMR_NPE], iopt[SPGMR_NPS]);
    }

    for (i = 0; i < Nloc; i++){
     usoln->data[i] = exp(-1.0*((thispe*Nloc)+i+1)*data->p[0]*tret);
     if (Ns > 0)
       usolnsub[1]->data[i]=-1.0*pbar[0]*((thispe*Nloc)+i+1)*tret*usoln->data[i]; 
    }
    N_VLinearSum(ONE, usoln, -ONE, uu, error);
    for (i = 0; i <= Ns; i++) {
      maxnorm = N_VMaxNorm(errorsub[i]);
      j = (plist == NULL) ? i : plist[i-1];
      if (thispe == 0)
	if (i == 0) 
	  printf("err(u)   = %10.5e\n", maxnorm);
	else
	  printf("err(s_%d) = %10.5e\n", j, maxnorm);
    }
    if (thispe == 0) printf("\n");

  }  /* End of tout loop. */

  /* Print remaining counters and free memory. */
  if (thispe == 0) printf("\n netf = %d,   ncfn = %d,   ncfl = %d \n", 
                   iopt[NETF], iopt[NCFN], iopt[SPGMR_NCFL]);

  SensIDAFree(idamem);
  N_VFree(uu);
  N_VFree(up);
  N_VFree(res);
  N_VFree(usoln);
  N_VFree(id);
  N_VFree(error);
  N_VFree(data->pp);

  if (plist != NULL) free(plist);
  free(pbar);
  free(data->p);
  free(data);
  PVecFreeMPI(machEnv);
  MPI_Finalize();
  return(0);

} /* End of sensdiagk main program. */

/***************** Functions called by the IDA solver ******************/

/*************************************************************************
 * resdiag: ODE system residual function                       
 ************************************************************************/

int resdiag(integer Ny, real tres, N_Vector uu, N_Vector up,
              N_Vector res, void *rdata)
{
  int retval;
  UserData data;
  integer i, thispe;
  real p1;
  real *uudata, *updata, *resdata;

  data = (UserData) rdata;
  thispe = data->thispe;
  p1 = data->p[0];

  /* compute the diagonal ODE residual function. */
  uudata = N_VDATA(uu);
  updata = N_VDATA(up);
  resdata = N_VDATA(res);
  for (i = 0; i < Nloc; i++) {
    resdata[i] = updata[i] + p1*((thispe*Nloc)+i+1)*uudata[i];
  }
  return(0);

} /* End of residual function resdiag. */

/*******************************************************************
 * PrecondSetup: setup for diagonal preconditioner for sensdiagk.  *
 *                                                                 *
 * The optional user-supplied functions PrecondSetup and           *
 * Precond together must define the left preconditioner            * 
 * matrix P approximating the system Jacobian matrix               *
 *                   J = dF/du + cj*dF/du'                         *
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   *
 * systems P z = r.   This is done in this case by using only the  *
 * 2nd term (cj*dF/du') of the matrix J above as a preconditioner. *
 * Since P is the diagonal matrix cj*I, its inverse is stored in a *
 * vector pp in PrecondSetup and subsequently used in Precond.     *
 *                                                                 *
 * In this instance, only cj and data (user data structure, with   * 
 * pp etc.) are used from the PrecondSetup argument list.          *
 ******************************************************************/
  
int PrecondSetup(integer local_N, real tt, N_Vector yy,
                  N_Vector yp, N_Vector rr, real cj,
                  ResFn res, void *rdata, void *pdata,
                  N_Vector ewt, N_Vector constraints, real hh, 
                  real uround, long int *nrePtr,
                  N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  integer i;
  real *ppv;
  UserData data;

  data = (UserData) pdata;
  ppv = N_VDATA(data->pp);

  N_VConst(ONE/cj, data->pp);

  return(SUCCESS);

} /* End of PrecondSetup. */


/******************************************************************
 * Precond: solve preconditioner linear system.                   *
 * This routine multiplies the input vector rvec by the vector pp *
 * containing the inverse diagonal Jacobian elements (previously  *
 * computed in PrecondSetup), returning the result in zvec.       */
  
int Precond(integer local_N, real tt, N_Vector uu,
                 N_Vector up, N_Vector rr, real cj, ResFn res, void *rdata,
                 void *pdata, N_Vector ewt, real delta, N_Vector rvec,
                 N_Vector zvec, long int *nrePtr, N_Vector tempv)
 {
  UserData data;

  data = (UserData) pdata;

  N_VProd(data->pp, rvec, zvec);

  return(SUCCESS);

} /* End of Precond. */
