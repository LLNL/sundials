/***********************************************************************
 * File:       sensiheatbbd.c   
 * Written by: Steven L. Lee and Alan C. Hindmarsh
 * Version of: 17 September 2001
 *----------------------------------------------------------------------
 *
 * Example problem, with sensitivity analysis.
 * SensIDA: 2D heat equation, parallel, GMRES, IDABBDPRE.
 *
 * This example solves a discretized 2D heat equation problem, along with the
 * sensitivity of the DAE solution with respect to several parameters.
 * This version uses the Krylov solver IDASpgmr and BBD preconditioning.
 *
 * The DAE system solved is a spatial discretization of the PDE 
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square.  The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MX x MY grid.
 * The values of u at the interior points satisfy ODEs, and equations
 * u = 0 at the boundaries are appended, to form a DAE system of size
 * N = MX * MY.  Here MX = MY = 10.
 *
 * The system is actually implemented on submeshes, processor by processor,
 * with an MXSUB by MYSUB mesh on each of NPEX * NPEY processors.
 *
 * The system is solved with SensIDA using the Krylov linear solver 
 * SENSIDASPGMR in conjunction with the preconditioner module IDABBDPRE. 
 * The preconditioner uses a tridiagonal approximation (half-bandwidths = 1).
 * The constraints u >= 0 are posed for all components.
 * Local error testing on the boundary values is suppressed. 
 * Output is taken at t = 0, .01, .02, .04, ..., 10.24.
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
#include "idabbdpre.h"
#include "mpi.h"
#include "sensida.h"

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)


#define NOUT         11             /* Number of output times */

#define NPEX         2              /* No. PEs in x direction of PE array */
#define NPEY         2              /* No. PEs in y direction of PE array */
                                    /* Total no. PEs = NPEX*NPEY */
#define MXSUB        5              /* No. x points per subgrid */
#define MYSUB        5              /* No. y points per subgrid */

#define MX           (NPEX*MXSUB)   /* MX = number of x mesh points */
#define MY           (NPEY*MYSUB)   /* MY = number of y mesh points */
                                    /* Spatial mesh is MX by MY */
#define NY           (MX*MY)        /* number of equations      */
#define NP           2              /* number of parameters     */
#define NS           2              /* number of sensitivities  */

typedef struct {  
  real      *p;
  integer   ny, thispe, mx, my, ixsub, jysub, npex, npey, mxsub, mysub;
  real      dx, dy, coeffx, coeffy, coeffxy;
  real      uext[(MXSUB+2)*(MYSUB+2)];
  N_Vector  pp;    /* vector of diagonal preconditioner elements */
  MPI_Comm  comm;
} *UserData;


/* Prototypes of private helper functions */

static int InitUserData(integer Ny, integer thispe, integer npes,
			MPI_Comm comm, UserData data);

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
			     N_Vector res, UserData data);


/* User-supplied residual function and supporting routines */


int heatres(integer Ny, real tres, N_Vector uu, N_Vector up,
            N_Vector res, void *rdata);

static int rescomm(N_Vector uu, N_Vector up, void *rdata);

static int reslocal(real tres, N_Vector uu, N_Vector up, 
		    N_Vector res,  void *rdata);

static int BSend(MPI_Comm comm, integer thispe, integer ixsub, integer jysub,
                 integer dsizex, integer dsizey, real uarray[]);

static int BRecvPost(MPI_Comm comm, MPI_Request request[], integer thispe,
		     integer ixsub, integer jysub,
		     integer dsizex, integer dsizey,
		     real uext[], real buffer[]);

static int BRecvWait(MPI_Request request[], integer ixsub, integer jysub,
		     integer dsizex, real uext[], real buffer[]);


main(int argc, char *argv[])

{
  integer retval, i, j, iout, itol, itask, local_N, npes, thispe;
  long int iopt[OPT_SIZE];
  integer mudq, mldq, mukeep, mlkeep;
  boole optIn;
  real ropt[OPT_SIZE], rtol, atol;
  real t0, t1, tout, tret, umax;
  void *mem;
  UserData data;
  N_Vector uu, up, constraints, id, res;
  IDAMem idamem;
  MPI_Comm comm;
  machEnvType machEnv;
  IBBDData P_data;
  N_Vector *uusub, *upsub;
  N_Vector *constraintssub;
  integer Ntotal;
  real *pbar;
  real rhomax;
  int *plist;

  Ntotal = (1+NS)*NY;

  /* Get processor number and total number of pe's. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &thispe);

  if (npes != NPEX*NPEY) {
    if (thispe == 0)
      printf("\n npes=%d is not equal to NPEX*NPEY=%d\n", npes,NPEX*NPEY);
    return(1);
  }

  /* Set local length local_N. */
  local_N = MXSUB*MYSUB;

  /* Set machEnv block. */
  machEnv = PVecInitMPI(comm, local_N, NY, &argc, &argv);
  if (machEnv == NULL) return(1);

  /* Allocate and initialize the data structure and N-vectors. */
  data = (UserData) malloc(sizeof *data);

  /* Store nominal parameter values in p */
  data->p = (real *) malloc(NP * sizeof(real));
  data->p[0] = 1.0;
  data->p[1] = 1.0;

  /* Scaling factor for each sensitivity equation */
  pbar = (real *) malloc(NP * sizeof(real));
  pbar[0] = 1.0;
  pbar[1] = 1.0;

  /* Store ordering of parameters in plist */
  plist = (int *) malloc(NP * sizeof(int));
  plist[0] = 1;
  plist[1] = 2;

  rhomax = 0.0;

  uu = N_VNew(Ntotal, machEnv); 
  up = N_VNew(Ntotal, machEnv);
  res = N_VNew(Ntotal, machEnv);
  constraints = N_VNew(Ntotal, machEnv);
  id = N_VNew(Ntotal, machEnv);

  /* Create pointers to subvectors */
  uusub = N_VSUB(uu);
  upsub = N_VSUB(up);
  constraintssub = N_VSUB(constraints);

  data->pp = N_VNew(NY, machEnv); /* An N-vector to hold preconditioner. */

  InitUserData(NY, thispe, npes, comm, data);

  /* Initialize the uu, up, id, and res profiles. */
  SetInitialProfile(uu, up, id, res, data);

  /* Set constraints to all 1's for nonnegative solution values in y. */
  N_VConst(ONE, constraintssub[0]);

  /* Initialize the sensitivity variables */
  SensInitZero(uu, NS);
  SensInitZero(up, NS);
  SensInitZero(constraints, NS);

  /* Identify all sensitivity variables as differential variables */
  SensSetId(id, NS);

  t0 = 0.0; t1 = 0.01;

  /* Scalar relative and absolute tolerance. */
  itol = SS;
  rtol = 0.0;
  atol = 1.e-3;

  /* Set option to suppress error testing of algebraic components. */
  optIn = TRUE;
  for (i = 0; i < OPT_SIZE; i++) {iopt[i] = 0; ropt[i] = ZERO; }
  iopt[SUPPRESSALG] = 1;

  /* Call SensIDAMalloc to initialize solution.  (NULL argument is errfp.) */
  itask = NORMAL;
  mem = SensIDAMalloc(NY, NS, Ntotal, heatres, data, t0, uu, up, 
		      itol, &rtol, &atol, id, constraints, NULL, optIn, 
		      iopt, ropt, machEnv, data->p, pbar, plist, rhomax);
  if (mem == NULL) {
    if (thispe == 0) printf ("SensIDAMalloc failed.");
    return(1); 
  }
  idamem = (IDAMem)mem;

  mudq = MXSUB;
  mldq = MXSUB;
  mukeep = 1;
  mlkeep = 1;

  /* Call IBBDAlloc to initialize BBD preconditioner. */
  P_data = IBBDAlloc(local_N, mudq, mldq, mukeep, mlkeep, ZERO, reslocal, 
                     rescomm, mem, data);

  /* Call SensIDASpgmr to specify the linear solver. */
  retval = SensIDASpgmr(idamem, IBBDPrecon, IBBDPSol, MODIFIED_GS, 
			0, 0, 0.0, 0.0, P_data);

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

  /* Compute the max norm of uu. */
  umax = N_VMaxNorm(uusub[0]);

  /* Print output heading (on processor 0 only). */

  if (thispe == 0) { 
    printf("sensheatbbd: Heat equation, parallel example problem for SensIDA \n");
    printf("             Discretized heat equation on 2D unit square. \n");
    printf("             Zero boundary conditions,");
    printf(" polynomial initial conditions.\n");
    printf("             Mesh dimensions: %d x %d", MX, MY);
    printf("             Total system size: %d\n\n", NY);
    printf("Subgrid dimensions: %d x %d", MXSUB, MYSUB);
    printf("             Processor array: %d x %d\n", NPEX, NPEY);
    printf("Number of sensitivities: Ns = %d\n", NS);
    printf("Parameter values:       p_1 = %9.2e,    p_2 = %9.2e\n", 
	   data->p[0], data->p[1]);
    printf("Scale factors:       pbar_1 = %9.2e, pbar_2 = %9.2e\n", 
	   pbar[0], pbar[1]);
    printf("Finite difference:   rhomax = %g\n", rhomax);
    printf("Tolerance parameters:  rtol = %g,  atol = %g\n", rtol, atol);
    printf("Constraints set to force all components of solution u >= 0. \n");
    printf("iopt[SUPPRESSALG] = 1 to suppress local error testing on");
    printf(" all boundary components. \n");
    printf("Linear solver: IDASPGMR  ");
    printf("Preconditioner: IDABBDPRE - Banded-block-diagonal.\n"); 
    printf("   Difference quotient half-bandwidths = %d",mudq);
    printf("   Retained matrix half-bandwidths = %d \n",mukeep);

    /* Print output table heading and initial line of table. */
    printf("\n");
    printf("Output Summary:   max(u) = max-norm of solution u \n");
    printf("                max(s_i) = max-norm of sensitivity vector s_i\n\n");
    printf("  time       max(u)     k   nst  nni  nli   nre    h       npe nps\n" );
    printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");

    printf(" %5.2f   %13.5e  %d  %3d  %3d  %3d  %4d %9.2e  %3d %3d\n",
           t0, umax, iopt[KUSED], iopt[NST], iopt[NNI], iopt[SPGMR_NLI], 
           iopt[NRE], ropt[HUSED], iopt[SPGMR_NPE], iopt[SPGMR_NPS]);
  }

  for (i = 1; i <= NS; i++){
    umax = N_VMaxNorm(uusub[i]);
    j = (plist == NULL) ? i : plist[i-1];
    if (thispe == 0) {
      printf("max(s_%d) = %11.5e\n", j, umax/pbar[j-1]);
      if (i == NS) printf("\n");
    }
  }

  /* Loop over tout, call IDASolve, print output. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    retval = IDASolve(idamem, tout, t0, &tret, uu, up, itask);

    umax = N_VMaxNorm(uusub[0]);
    if (thispe == 0) printf(" %5.2f   %13.5e  %d  %3d  %3d  %3d  %4d %9.2e  %3d %3d\n",
                     tret, umax, iopt[KUSED], iopt[NST], iopt[NNI], iopt[SPGMR_NLI], 
                     iopt[NRE], ropt[HUSED], iopt[SPGMR_NPE], iopt[SPGMR_NPS]);

    for (i = 1; i <= NS; i++){
      umax = N_VMaxNorm(uusub[i]);
      j = (plist == NULL) ? i: plist[i-1];
      if (thispe == 0) {
	printf("max(s_%d) = %11.5e\n", j, umax/pbar[j-1]);
	if (i == NS) printf("\n");
      }
    }

    if (retval < 0) {
      if (thispe == 0) printf("IDASolve returned %d.\n",retval);
      return(1);
    }

  }  /* End of tout loop. */

  /* Print remaining counters and free memory. */
  if (thispe == 0) printf("\n netf = %d,   ncfn = %d,   ncfl = %d \n", 
                   iopt[NETF], iopt[NCFN], iopt[SPGMR_NCFL]);

  IBBDFree(P_data);
  SensIDAFree(idamem);
  N_VFree(uu);
  N_VFree(up);
  N_VFree(constraints);
  N_VFree(id);
  N_VFree(res);
  N_VFree(data->pp);

  if (plist != NULL) free(plist);
  free(pbar);
  free(data->p);
  free(data);
  PVecFreeMPI(machEnv);
  MPI_Finalize();
  return(0);

} /* End of sensheatbbd main program. */


/********************************************************/
/* InitUserData initializes the user's data block data. */

static int InitUserData(integer Ny, integer thispe, integer npes,
		        MPI_Comm comm, UserData data)
{

  data->ny = Ny;
  data->thispe = thispe;
  data->dx = ONE/(MX-ONE);       /* Assumes a [0,1] interval in x. */
  data->dy = ONE/(MY-ONE);       /* Assumes a [0,1] interval in y. */
  data->coeffx  = ONE/(data->dx * data->dx);
  data->coeffy  = ONE/(data->dy * data->dy);
  data->coeffxy = TWO/(data->dx * data->dx) + TWO/(data->dy * data->dy) ;
  data->jysub   = thispe/NPEX;
  data->ixsub   = thispe - data->jysub * NPEX;
  data->npex    = NPEX;
  data->npey    = NPEY;
  data->mx      = MX;
  data->my      = MY;
  data->mxsub = MXSUB;
  data->mysub = MYSUB;
  data->comm    = comm;
  return(0);

} /* End of InitUserData. */


/**************************************************************/
/* SetInitialProfile sets the initial values for the problem. */

static int SetInitialProfile(N_Vector uu, N_Vector up,  N_Vector id, 
			     N_Vector res, UserData data)
{
  integer i, iloc, j, jloc, offset, loc, ixsub, jysub;
  integer ixbegin, ixend, jybegin, jyend;
  real xfact, yfact, *udata, *iddata, dx, dy;

  /* Initialize uu. */ 

  udata = N_VDATA(uu);
  iddata = N_VDATA(id);

  /* Set mesh spacings and subgrid indices for this PE. */
  dx = data->dx;
  dy = data->dy;
  ixsub = data->ixsub;
  jysub = data->jysub;

  /* Set beginning and ending locations in the global array corresponding 
     to the portion of that array assigned to this processor. */
  ixbegin = MXSUB*ixsub;
  ixend   = MXSUB*(ixsub+1) - 1;
  jybegin = MYSUB*jysub;
  jyend   = MYSUB*(jysub+1) - 1;

 /* Loop over the local array, computing the initial profile value.
    The global indices are (i,j) and the local indices are (iloc,jloc).
    Also set the id vector to zero for boundary points, one otherwise. */

  N_VConst(ONE, id);
  for (j = jybegin, jloc = 0; j <= jyend; j++, jloc++) {
    yfact = data->dy*j;
    offset= jloc*MXSUB;
    for (i = ixbegin, iloc = 0; i <= ixend; i++, iloc++) {
      xfact = data->dx * i;
      loc = offset + iloc;
      udata[loc] = 16. * xfact * (ONE - xfact) * yfact * (ONE - yfact);
      if (i == 0 || i == MX-1 || j == 0 || j == MY-1) iddata[loc] = ZERO;
      }
    }

  /* Initialize up. */

  N_VConst(ZERO, up);    /* Initially set up = 0. */

  /* heatres sets res to negative of ODE RHS values at interior points. */
  heatres(data->ny, ZERO, uu, up, res, data);

  /* Copy -res into up to get correct initial up values. */
  N_VScale(-ONE, res, up);

  return(SUCCESS);

} /* End of SetInitialProfiles. */


/***************** Functions called by the IDA solver ******************/


/*************************************************************************
 * heatres: heat equation system residual function                       
 * This uses 5-point central differencing on the interior points, and    
 * includes algebraic equations for the boundary values.                 
 * So for each interior point, the residual component has the form       
 *    res_i = u'_i - (central difference)_i                              
 * while for each boundary point, it is res_i = u_i. 
 *                    
 * This parallel implementation uses several supporting routines. 
 * First a call is made to rescomm to do communication of subgrid boundary
 * data into array uext.  Then reslocal is called to compute the residual
 * on individual processors and their corresponding domains.  The routines
 * BSend, BRecvPost, and BREcvWait handle interprocessor communication
 * of uu required to calculate the residual. */

int heatres(integer Ny, real tres, N_Vector uu, N_Vector up,
              N_Vector res, void *rdata)
{
  int retval;
  UserData data;

  data = (UserData) rdata;

  /* Call rescomm to do inter-processor communication. */
  retval = rescomm(uu, up, data);

  /* Call reslocal to calculate res. */
  retval = reslocal(tres, uu, up, res, data);
  
  return(0);

} /* End of residual function heatres. */


/*  Supporting functions for heatres. */


/*********************************************************************/
/* rescomm routine.  This routine performs all inter-processor
   communication of data in u needed to calculate G.                 */

static int rescomm(N_Vector uu, N_Vector up, void *rdata)
{
  UserData data;
  real *uarray, *uext, buffer[2*MYSUB];
  MPI_Comm comm;
  integer thispe, ixsub, jysub, mxsub, mysub;
  MPI_Request request[4];

  data = (UserData) rdata;
  uarray = N_VDATA(uu);

  /* Get comm, thispe, subgrid indices, data sizes, extended array uext. */
  comm = data->comm;  thispe = data->thispe;
  ixsub = data->ixsub;   jysub = data->jysub;
  mxsub = data->mxsub;   mysub = data->mysub;
  uext = data->uext;

  /* Start receiving boundary data from neighboring PEs. */
  BRecvPost(comm, request, thispe, ixsub, jysub, mxsub, mysub, uext, buffer);

  /* Send data from boundary of local grid to neighboring PEs. */
  BSend(comm, thispe, ixsub, jysub, mxsub, mysub, uarray);

  /* Finish receiving boundary data from neighboring PEs. */
  BRecvWait(request, ixsub, jysub, mxsub, uext, buffer);

  return(0);

} /* End of rescomm. */


/**************************************************************************/
/* reslocal routine.  Compute res = F(t, uu, up).  This routine assumes
   that all inter-processor communication of data needed to calculate F
    has already been done, and that this data is in the work array uext.  */

static int reslocal(real tres, N_Vector uu, N_Vector up, N_Vector res,
                    void *rdata)
{
  real *uext, *uuv, *upv, *resv;
  real termx, termy, termctr;
  integer i, lx, j, ly, offsetu, offsetue, locu, locue;
  integer ixsub, jysub, mxsub, mxsub2, mysub, npex, npey;
  integer ixbegin, ixend, jybegin, jyend;
  UserData data;
  real p1, p2;

  /* Get subgrid indices, array sizes, extended work array uext. */

  data = (UserData) rdata;
  p1 = data->p[0];
  p2 = data->p[1];
  uext = data->uext;
  uuv = N_VDATA(uu);
  upv = N_VDATA(up);
  resv = N_VDATA(res);
  ixsub = data->ixsub; jysub = data->jysub;
  mxsub = data->mxsub; mxsub2 = data->mxsub + 2;
  mysub = data->mysub; npex = data->npex; npey = data->npey;

  /* Initialize all elements of res to uu. This sets the boundary
     elements simply without indexing hassles. */

  N_VScale(ONE, uu, res);

  /* Copy local segment of u vector into the working extended array uext.
     This completes uext prior to the computation of the res vector.     */

  offsetu = 0;
  offsetue = mxsub2 + 1;
  for (ly = 0; ly < mysub; ly++) {
    for (lx = 0; lx < mxsub; lx++) uext[offsetue+lx] = uuv[offsetu+lx];
    offsetu = offsetu + mxsub;
    offsetue = offsetue + mxsub2;
  }

  /* Set loop limits for the interior of the local subgrid. */

  ixbegin = 0;
  ixend   = mxsub-1;
  jybegin = 0;
  jyend   = mysub-1;
  if (ixsub == 0) ixbegin++; if (ixsub == npex-1) ixend--;
  if (jysub == 0) jybegin++; if (jysub == npey-1) jyend--;
  
  /* Loop over all grid points in local subgrid. */

  for (ly = jybegin; ly <=jyend; ly++) {
    for (lx = ixbegin; lx <= ixend; lx++) {
      locu  = lx + ly*mxsub;
      locue = (lx+1) + (ly+1)*mxsub2;
      termx = p1 * data->coeffx *(uext[locue-1]      + uext[locue+1]);
      termy = p2 * data->coeffy *(uext[locue-mxsub2] + uext[locue+mxsub2]);
      termctr = (p1*(data->coeffx) + p2*(data->coeffy)) * TWO * uext[locue];
      resv[locu] = upv[locu] - (termx + termy - termctr);
   }
  }
  return(0);

} /* End of reslocal. */


/*************************************************************************/
/* Routine to send boundary data to neighboring PEs.                     */

static int BSend(MPI_Comm comm, integer thispe, integer ixsub, integer jysub,
                 integer dsizex, integer dsizey, real uarray[])
{
  integer ly, offsetu;
  real bufleft[MYSUB], bufright[MYSUB];

  /* If jysub > 0, send data from bottom x-line of u. */

  if (jysub != 0)
    MPI_Send(&uarray[0], dsizex, PVEC_REAL_MPI_TYPE, thispe-NPEX, 0, comm);

  /* If jysub < NPEY-1, send data from top x-line of u. */

  if (jysub != NPEY-1) {
    offsetu = (MYSUB-1)*dsizex;
    MPI_Send(&uarray[offsetu], dsizex, PVEC_REAL_MPI_TYPE, 
           thispe+NPEX, 0, comm);
  }

  /* If ixsub > 0, send data from left y-line of u (via bufleft). */

  if (ixsub != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*dsizex;
      bufleft[ly] = uarray[offsetu];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, thispe-1, 0, comm);   
  }

  /* If ixsub < NPEX-1, send data from right y-line of u (via bufright). */

  if (ixsub != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*MXSUB + (MXSUB-1);
      bufright[ly] = uarray[offsetu];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, thispe+1, 0, comm);   
  }
} /* End of BSend. */


/**************************************************************/
/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*MYSUB real entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in 
   both calls also. */

static int BRecvPost(MPI_Comm comm, MPI_Request request[], integer thispe,
		     integer ixsub, integer jysub,
		     integer dsizex, integer dsizey,
		     real uext[], real buffer[])
{
  integer offsetue;
  /* Have bufleft and bufright use the same buffer. */
  real *bufleft = buffer, *bufright = buffer+MYSUB;

  /* If jysub > 0, receive data for bottom x-line of uext. */
  if (jysub != 0)
    MPI_Irecv(&uext[1], dsizex, PVEC_REAL_MPI_TYPE,
    					 thispe-NPEX, 0, comm, &request[0]);

  /* If jysub < NPEY-1, receive data for top x-line of uext. */
  if (jysub != NPEY-1) {
    offsetue = (1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&uext[offsetue], dsizex, PVEC_REAL_MPI_TYPE,
                                         thispe+NPEX, 0, comm, &request[1]);
  }

  /* If ixsub > 0, receive data for left y-line of uext (via bufleft). */
  if (ixsub != 0) {
    MPI_Irecv(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         thispe-1, 0, comm, &request[2]);
  }

 /* If ixsub < NPEX-1, receive data for right y-line of uext (via bufright). */
  if (ixsub != NPEX-1) {
    MPI_Irecv(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         thispe+1, 0, comm, &request[3]);
  }

} /* End of BRecvPost. */


/****************************************************************/
/* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*MYSUB real entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have four entries, and should be passed in both 
      calls also. */

static int BRecvWait(MPI_Request request[], integer ixsub, integer jysub,
		      integer dsizex, real uext[], real buffer[])
{
  integer ly, dsizex2, offsetue;
  real *bufleft = buffer, *bufright = buffer+MYSUB;
  MPI_Status status;

  dsizex2 = dsizex + 2;

  /* If jysub > 0, receive data for bottom x-line of uext. */
  if (jysub != 0)
    MPI_Wait(&request[0],&status);

  /* If jysub < NPEY-1, receive data for top x-line of uext. */
  if (jysub != NPEY-1)
    MPI_Wait(&request[1],&status);

  /* If ixsub > 0, receive data for left y-line of uext (via bufleft). */
  if (ixsub != 0) {
    MPI_Wait(&request[2],&status);

    /* Copy the buffer to uext. */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetue = (ly+1)*dsizex2;
      uext[offsetue] = bufleft[ly];
    }
  }

  /* If ixsub < NPEX-1, receive data for right y-line of uext (via bufright). */
  if (ixsub != NPEX-1) {
    MPI_Wait(&request[3],&status);

    /* Copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetue = (ly+2)*dsizex2 - 1;
      uext[offsetue] = bufright[ly];
    }
  }
} /* End of BRecvWait. */
