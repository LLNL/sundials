/***********************************************************************
 * File       : iheatpk.c   
 * Written by : Allan G. Taylor Alan C. Hindmarsh, and Radu Serban
 * Version of : 23 July 2003
 *----------------------------------------------------------------------
 *
 * Example problem for IDAS: 2D heat equation, parallel, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver IDASpgmr.
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
 * The system is solved with IDAS using the Krylov linear solver IDASPGMR. 
 * The preconditioner uses the diagonal elements of the Jacobian only.
 * Routines for preconditioning, required by IDASPGMR, are supplied here.
 * The constraints u >= 0 are posed for all components.
 * Local error testing on the boundary values is suppressed. 
 * Output is taken at t = 0, .01, .02, .04, ..., 10.24.
 *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"
#include "sundialsmath.h"
#include "nvector_parallel.h" /* definitions of type N_Vector, macro NV_DATA_P   */
#include "idas.h"
#include "idasspgmr.h"
#include "iterativ.h"
#include "mpi.h"


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

typedef struct {  
  integertype thispe, mx, my, ixsub, jysub, npex, npey, mxsub, mysub;
  realtype    dx, dy, coeffx, coeffy, coeffxy;
  realtype    uext[(MXSUB+2)*(MYSUB+2)];
  N_Vector    pp;    /* vector of diagonal preconditioner elements */
  MPI_Comm    comm;
} *UserData;


/* Prototypes of private helper functions */

static int InitUserData(int thispe, MPI_Comm comm, UserData data);

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, UserData data);


/* User-supplied residual function and supporting routines */


int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector res, void *rdata);

static int rescomm(N_Vector uu, N_Vector up, void *rdata);

static int reslocal(realtype tres, N_Vector uu, N_Vector up, 
                    N_Vector res,  void *rdata);

static int BSend(MPI_Comm comm, integertype thispe, integertype ixsub, integertype jysub,
                 integertype dsizex, integertype dsizey, realtype uarray[]);

static int BRecvPost(MPI_Comm comm, MPI_Request request[], integertype thispe,
                     integertype ixsub, integertype jysub,
                     integertype dsizex, integertype dsizey,
                     realtype uext[], realtype buffer[]);

static int BRecvWait(MPI_Request request[], integertype ixsub, integertype jysub,
                     integertype dsizex, realtype uext[], realtype buffer[]);

/* User-supplied preconditioner routines */

int PSolveHeateq(realtype tt, N_Vector uu,
                 N_Vector up, N_Vector rr, 
                 N_Vector rvec, N_Vector zvec,
                 realtype cj, realtype delta,
                 void *pdata, N_Vector tempv);

int PrecondHeateq(realtype tt, N_Vector yy,
                  N_Vector yp, N_Vector rr, 
                  realtype cj, void *pdata,
                  N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);


int main(int argc, char *argv[])

{
  MPI_Comm comm;
  NV_Spec nvSpec;
  void *mem;
  int iout, itol, itask, npes, thispe;
  integertype Neq, local_N;
  realtype rtol, atol;
  realtype t0, t1, tout, tret, umax, hused;
  UserData data;
  N_Vector uu, up, constraints, id, res;
  int ier, kused, nst, nni, nre, netf, ncfn, ncfl, nli, npe, nps, nreS;

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
  
  /* Set local length local_N and global length Neq. */
  local_N = MXSUB*MYSUB;
  Neq     = MX * MY;
  
  /* Set nvSpec block. */
  nvSpec = NV_SpecInit_Parallel(comm, local_N, Neq, &argc, &argv);
  if (nvSpec == NULL) return(1);
  
  /* Allocate and initialize the data structure and N-vectors. */
  data = (UserData) malloc(sizeof *data);
  
  uu = N_VNew(nvSpec); 
  up = N_VNew(nvSpec);
  res = N_VNew(nvSpec);
  constraints = N_VNew(nvSpec);
  id = N_VNew(nvSpec);
  data->pp = N_VNew(nvSpec); /* An N-vector to hold preconditioner. */
  
  InitUserData(thispe, comm, data);
  
  /* Initialize the uu, up, id, and res profiles. */
  SetInitialProfile(uu, up, id, res, data);
  
  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);
  
  t0 = 0.0; t1 = 0.01;
  
  /* Scalar relative and absolute tolerance. */
  itol = SS;
  rtol = 0.0;
  atol = 1.e-3;


  /* Call IDACreate and IDAMalloc to initialize solution. */
  mem = IDACreate();
  if (mem == NULL) {
    if (thispe == 0) printf ("IDACreate failed.");
    return(1); 
  }
  ier = IDASetRdata(mem, data);
  ier = IDASetSuppressAlg(mem, TRUE);
  ier = IDASetID(mem, id);
  ier = IDASetConstraints(mem, constraints);
  ier = IDAMalloc(mem, heatres, t0, uu, up, itol, &rtol, &atol, nvSpec);
  
  /* Call IDASpgmr to specify the linear solver. */
  /* Call IDASpgmr to specify the linear solver. */
  ier = IDASpgmr(mem, 0.0);
  if (ier != SUCCESS) {
    if (thispe == 0) printf("IDASpgmr failed, returning %d.\n",ier);
    return(1);
  }
  ier = IDASpgmrSetPrecSetupFn(mem, PrecondHeateq);
  ier = IDASpgmrSetPrecSolveFn(mem, PSolveHeateq);
  ier = IDASpgmrSetPrecData(mem, data);
  
  /* Compute the max norm of uu. */
  umax = N_VMaxNorm(uu);
  
  /* Print output heading (on processor 0 only). */
  
  if (thispe == 0) { 
    printf("iheatpk: Heat equation, parallel example problem for IDAS \n");
    printf("         Discretized heat equation on 2D unit square. \n");
    printf("         Zero boundary conditions,");
    printf(" polynomial initial conditions.\n");
    printf("         Mesh dimensions: %d x %d", MX, MY);
    printf("        Total system size: %ld\n\n", Neq);
    printf("Subgrid dimensions: %d x %d", MXSUB, MYSUB);
    printf("        Processor array: %d x %d\n", NPEX, NPEY);
    printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
    printf("Constraints set to force all solution components >= 0. \n");
    printf("SUPPRESSALG = TRUE to suppress local error testing on ");
    printf("all boundary components. \n");
    printf("Linear solver: IDASPGMR  ");
    printf("Preconditioner: diagonal elements only.\n"); 

    /* Print output table heading and initial line of table. */
    printf("\n   Output Summary (umax = max-norm of solution) \n\n");
    printf("  time     umax       k  nst  nni  nli   nre   nreS    h      npe nps\n");
    printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");

    printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d %4d  %9.2e  %3d %3d\n",
           t0, umax, 0, 0, 0, 0, 0, 0, 0.0, 0, 0);
  }

  /* Loop over tout, call IDASolve, print output. */
  itask = NORMAL;
  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {
    
    ier = IDASolve(mem, tout, &tret, uu, up, itask);
    
    umax = N_VMaxNorm(uu);
    if (thispe == 0) {
      IDAGetLastOrder(mem, &kused);
      IDAGetNumSteps(mem, &nst);
      IDAGetNumNonlinSolvIters(mem, &nni);
      IDAGetNumResEvals(mem, &nre);
      IDAGetLastStep(mem, &hused);
      IDASpgmrGetNumLinIters(mem, &nli);
      IDASpgmrGetNumPrecEvals(mem, &npe);
      IDASpgmrGetNumPrecSolves(mem, &nps);
      IDASpgmrGetNumResEvals(mem, &nreS);
      printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d %4d  %9.2e  %3d %3d\n",
             tret, umax, kused, nst, nni, nli, nre, nreS, hused, npe, nps);
    }
    if (ier < 0) {
      if (thispe == 0) printf("IDASolve returned %d.\n",ier);
      return(1);
    }
    
  }  /* End of tout loop. */
  
  /* Print remaining counters and free memory. */
  IDAGetNumErrTestFails(mem, &netf);
  IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  IDASpgmrGetNumConvFails(mem, &ncfl);
  printf("\n netf = %d,   ncfn = %d,   ncfl = %d \n", netf, ncfn, ncfl);


  IDAFree(mem);
  N_VFree(data->pp);
  free(data);
  N_VFree(id);
  N_VFree(constraints);
  N_VFree(res);
  N_VFree(up);
  N_VFree(uu);
  NV_SpecFree_Parallel(nvSpec);
  MPI_Finalize();

  return(0);

} /* End of iheatpk main program. */


/********************************************************/
/* InitUserData initializes the user's data block data. */

static int InitUserData(int thispe, MPI_Comm comm, UserData data)
{
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
  integertype i, iloc, j, jloc, offset, loc, ixsub, jysub;
  integertype ixbegin, ixend, jybegin, jyend;
  realtype xfact, yfact, *udata, *iddata, dx, dy;
  
  /* Initialize uu. */ 

  udata = NV_DATA_P(uu);
  iddata = NV_DATA_P(id);
  
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
  
  N_VConst(ONE,id);
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
  heatres(ZERO, uu, up, res, data);
  
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

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector res, void *rdata)
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
  realtype *uarray, *uext, buffer[2*MYSUB];
  MPI_Comm comm;
  integertype thispe, ixsub, jysub, mxsub, mysub;
  MPI_Request request[4];
  
  data = (UserData) rdata;
  uarray = NV_DATA_P(uu);
  
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

static int reslocal(realtype tres, N_Vector uu, N_Vector up, N_Vector res,
                    void *rdata)
{
  realtype *uext, *uuv, *upv, *resv;
  realtype termx, termy, termctr;
  integertype lx, ly, offsetu, offsetue, locu, locue;
  integertype ixsub, jysub, mxsub, mxsub2, mysub, npex, npey;
  integertype ixbegin, ixend, jybegin, jyend;
  UserData data;
  
  /* Get subgrid indices, array sizes, extended work array uext. */
  
  data = (UserData) rdata;
  uext = data->uext;
  uuv = NV_DATA_P(uu);
  upv = NV_DATA_P(up);
  resv = NV_DATA_P(res);
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
      termx = data->coeffx *(uext[locue-1]      + uext[locue+1]);
      termy = data->coeffy *(uext[locue-mxsub2] + uext[locue+mxsub2]);
      termctr = data->coeffxy*uext[locue];
      resv[locu] = upv[locu] - (termx + termy - termctr);
   }
  }
  return(0);

} /* End of reslocal. */


/*************************************************************************/
/* Routine to send boundary data to neighboring PEs.                     */

static int BSend(MPI_Comm comm, integertype thispe, integertype ixsub, integertype jysub,
                 integertype dsizex, integertype dsizey, realtype uarray[])
{
  integertype ly, offsetu;
  realtype bufleft[MYSUB], bufright[MYSUB];

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

  return(0);

} /* End of BSend. */


/**************************************************************/
/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in 
   both calls also. */

static int BRecvPost(MPI_Comm comm, MPI_Request request[], integertype thispe,
                     integertype ixsub, integertype jysub,
                     integertype dsizex, integertype dsizey,
                     realtype uext[], realtype buffer[])
{
  integertype offsetue;
  /* Have bufleft and bufright use the same buffer. */
  realtype *bufleft = buffer, *bufright = buffer+MYSUB;
  
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

  return(0);
  
} /* End of BRecvPost. */


/****************************************************************/
/* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have four entries, and should be passed in both 
   calls also. */

static int BRecvWait(MPI_Request request[], integertype ixsub, integertype jysub,
                     integertype dsizex, realtype uext[], realtype buffer[])
{
  integertype ly, dsizex2, offsetue;
  realtype *bufleft = buffer, *bufright = buffer+MYSUB;
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

  return(0);

} /* End of BRecvWait. */


/*******************************************************************
 * PrecondHeateq: setup for diagonal preconditioner for heatsk.    *
 *                                                                 *
 * The optional user-supplied functions PrecondHeateq and          *
 * PSolveHeateq together must define the left preconditoner        *
 * matrix P approximating the system Jacobian matrix               *
 *                   J = dF/du + cj*dF/du'                         *
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   *
 * systems P z = r.   This is done in this case by keeping only    *
 * the diagonal elements of the J matrix above, storing them as    *
 * inverses in a vector pp, when computed in PrecondHeateq, for    *
 * subsequent use in PSolveHeateq.                                 *
 *                                                                 *
 * In this instance, only cj and data (user data structure, with   * 
 * pp etc.) are used from the PrecondHeateq argument list.         *
 ******************************************************************/
int PrecondHeateq(realtype tt, N_Vector yy,
                  N_Vector yp, N_Vector rr, 
                  realtype cj, void *pdata,
                  N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *ppv, pelinv;
  integertype lx, ly, ixbegin, ixend, jybegin, jyend, locu, mxsub, mysub;
  integertype ixsub, jysub, npex, npey;
  UserData data;

  data = (UserData) pdata;

  ppv = NV_DATA_P(data->pp);
  ixsub = data->ixsub;
  jysub = data->jysub;
  mxsub = data->mxsub;
  mysub = data->mysub;
  npex  = data->npex;
  npey  = data->npey;
  
  /* Initially set all pp elements to one. */
  N_VConst(ONE, data->pp);
  
  /* Prepare to loop over subgrid. */
  ixbegin = 0;
  ixend   = mxsub-1;
  jybegin = 0;
  jyend   = mysub-1;
  if (ixsub == 0) ixbegin++; if (ixsub == npex-1) ixend--;
  if (jysub == 0) jybegin++; if (jysub == npey-1) jyend--;
  pelinv = ONE/(cj + data->coeffxy); 
  
  /* Load the inverse of the preconditioner diagonal elements
     in loop over all the local subgrid.                           */
  
  for (ly = jybegin; ly <=jyend; ly++) {
    for (lx = ixbegin; lx <= ixend; lx++) {
      locu  = lx + ly*mxsub;
      ppv[locu] = pelinv;
    }
  }

  return(SUCCESS);

} /* End of PrecondHeateq. */


/******************************************************************
 * PSolveHeateq: solve preconditioner linear system.              *
 * This routine multiplies the input vector rvec by the vector pp *
 * containing the inverse diagonal Jacobian elements (previously  *
 * computed in PrecondHeateq), returning the result in zvec.      */

int PSolveHeateq(realtype tt, N_Vector uu,
                 N_Vector up, N_Vector rr, 
                 N_Vector rvec, N_Vector zvec,
                 realtype cj, realtype delta,
                 void *pdata, N_Vector tempv)
{
  UserData data;

  data = (UserData) pdata;
  
  N_VProd(data->pp,rvec, zvec);

  return(SUCCESS);

} /* End of PSolveHeateq. */
