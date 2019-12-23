/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *         Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, parallel, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver SUNLinSol_SPGMR.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MX x MY
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = MX * MY. Here MX = MY = 10.
 *
 * The system is actually implemented on submeshes, processor by
 * processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
 * processors.
 *
 * The system is solved with IDA using the Krylov linear solver
 * SUNLinSol_SPGMR. The preconditioner uses the diagonal elements of the
 * Jacobian only. Routines for preconditioning, required by
 * SUNLinSol_SPGMR, are supplied here. The constraints u >= 0 are posed
 * for all components. Local error testing on the boundary values
 * is suppressed. Output is taken at t = 0, .01, .02, .04,
 * ..., 10.24.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <idas/idas.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>

#include <mpi.h>

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
  int thispe, npex, npey, ixsub, jysub;
  sunindextype mx, my, mxsub, mysub;
  realtype    dx, dy, coeffx, coeffy, coeffxy;
  realtype    uext[(MXSUB+2)*(MYSUB+2)];
  N_Vector    pp;    /* vector of diagonal preconditioner elements */
  MPI_Comm    comm;
} *UserData;

/* User-supplied residual function and supporting routines */

int resHeat(realtype tt, N_Vector uu, N_Vector up,
            N_Vector rr, void *user_data);

static int rescomm(N_Vector uu, N_Vector up, void *user_data);

static int reslocal(realtype tt, N_Vector uu, N_Vector up,
                    N_Vector res,  void *user_data);

static int BSend(MPI_Comm comm, int thispe, int ixsub, int jysub,
                 sunindextype dsizex, sunindextype dsizey, realtype uarray[]);

static int BRecvPost(MPI_Comm comm, MPI_Request request[], int thispe,
                     int ixsub, int jysub, sunindextype dsizex,
                     sunindextype dsizey, realtype uext[], realtype buffer[]);

static int BRecvWait(MPI_Request request[], int ixsub, int jysub,
                     sunindextype dsizex, realtype uext[], realtype buffer[]);

/* User-supplied preconditioner routines */

int PsolveHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
               N_Vector rvec, N_Vector zvec, realtype c_j,
               realtype delta, void *user_data);

int PsetupHeat(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
               realtype c_j, void *user_data);

/* Private function to check function return values */

static int InitUserData(int thispe, MPI_Comm comm, UserData data);

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, UserData data);

static void PrintHeader(sunindextype Neq, realtype rtol, realtype atol);

static void PrintOutput(int id, void *ida_mem, realtype t, N_Vector uu);

static void PrintFinalStats(void *ida_mem);

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  MPI_Comm comm;
  void *ida_mem;
  SUNLinearSolver LS;
  UserData data;
  int iout, thispe, retval, npes;
  sunindextype Neq, local_N;
  realtype rtol, atol, t0, t1, tout, tret;
  N_Vector uu, up, constraints, id, res;

  ida_mem = NULL;
  LS = NULL;
  data = NULL;
  uu = up = constraints = id = res = NULL;

  /* Get processor number and total number of pe's. */

  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &thispe);

  if (npes != NPEX*NPEY) {
    if (thispe == 0)
      fprintf(stderr,
              "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n",
              npes,NPEX*NPEY);
    MPI_Finalize();
    return(1);
  }

  /* Set local length local_N and global length Neq. */

  local_N = MXSUB*MYSUB;
  Neq     = MX * MY;

  /* Allocate and initialize the data structure and N-vectors. */

  data = (UserData) malloc(sizeof *data);
  data->pp = NULL;
  if(check_retval((void *)data, "malloc", 2, thispe))
    MPI_Abort(comm, 1);

  uu = N_VNew_Parallel(comm, local_N, Neq);
  if(check_retval((void *)uu, "N_VNew_Parallel", 0, thispe))
    MPI_Abort(comm, 1);

  up = N_VNew_Parallel(comm, local_N, Neq);
  if(check_retval((void *)up, "N_VNew_Parallel", 0, thispe))
    MPI_Abort(comm, 1);

  res = N_VNew_Parallel(comm, local_N, Neq);
  if(check_retval((void *)res, "N_VNew_Parallel", 0, thispe))
    MPI_Abort(comm, 1);

  constraints = N_VNew_Parallel(comm, local_N, Neq);
  if(check_retval((void *)constraints, "N_VNew_Parallel", 0, thispe))
    MPI_Abort(comm, 1);

  id = N_VNew_Parallel(comm, local_N, Neq);
  if(check_retval((void *)id, "N_VNew_Parallel", 0, thispe))
    MPI_Abort(comm, 1);

  /* An N-vector to hold preconditioner. */
  data->pp = N_VNew_Parallel(comm, local_N, Neq);
  if(check_retval((void *)data->pp, "N_VNew_Parallel", 0, thispe))
    MPI_Abort(comm, 1);

  InitUserData(thispe, comm, data);

  /* Initialize the uu, up, id, and res profiles. */

  SetInitialProfile(uu, up, id, res, data);

  /* Set constraints to all 1's for nonnegative solution values. */

  N_VConst(ONE, constraints);

  t0 = ZERO; t1 = RCONST(0.01);

  /* Scalar relative and absolute tolerance. */

  rtol = ZERO;
  atol = RCONST(1.0e-3);

  /* Call IDACreate and IDAMalloc to initialize solution. */

  ida_mem = IDACreate();
  if(check_retval((void *)ida_mem, "IDACreate", 0, thispe)) MPI_Abort(comm, 1);

  retval = IDASetUserData(ida_mem, data);
  if(check_retval(&retval, "IDASetUserData", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetSuppressAlg(ida_mem, SUNTRUE);
  if(check_retval(&retval, "IDASetSuppressAlg", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetId(ida_mem, id);
  if(check_retval(&retval, "IDASetId", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetConstraints(ida_mem, constraints);
  if(check_retval(&retval, "IDASetConstraints", 1, thispe)) MPI_Abort(comm, 1);
  N_VDestroy(constraints);

  retval = IDAInit(ida_mem, resHeat, t0, uu, up);
  if(check_retval(&retval, "IDAInit", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASStolerances(ida_mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1, thispe)) MPI_Abort(comm, 1);

  /* Call SUNLinSol_SPGMR and IDASetLinearSolver to specify the linear solver. */

  LS = SUNLinSol_SPGMR(uu, PREC_LEFT, 0);  /* use default maxl */
  if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0, thispe)) MPI_Abort(comm, 1);

  retval = IDASetLinearSolver(ida_mem, LS, NULL);
  if(check_retval(&retval, "IDASetLinearSolver", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetPreconditioner(ida_mem, PsetupHeat, PsolveHeat);
  if(check_retval(&retval, "IDASetPreconditioner", 1, thispe)) MPI_Abort(comm, 1);

  /* Print output heading (on processor 0 only) and intial solution  */

  if (thispe == 0) PrintHeader(Neq, rtol, atol);
  PrintOutput(thispe, ida_mem, t0, uu);

  /* Loop over tout, call IDASolve, print output. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    retval = IDASolve(ida_mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1, thispe)) MPI_Abort(comm, 1);

    PrintOutput(thispe, ida_mem, tret, uu);

  }

  /* Print remaining counters. */

  if (thispe == 0) PrintFinalStats(ida_mem);

  /* Free memory */

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);

  N_VDestroy(id);
  N_VDestroy(res);
  N_VDestroy(up);
  N_VDestroy(uu);

  N_VDestroy(data->pp);
  free(data);

  MPI_Finalize();

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resHeat: heat equation system residual function
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
 * of uu required to calculate the residual.
 */

int resHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
            void *user_data)
{
  int retval;

  /* Call rescomm to do inter-processor communication. */
  retval = rescomm(uu, up, user_data);

  /* Call reslocal to calculate res. */
  retval = reslocal(tt, uu, up, rr, user_data);

  return(retval);

}

/*
 * PsetupHeat: setup for diagonal preconditioner for heatsk.
 *
 * The optional user-supplied functions PsetupHeat and
 * PsolveHeat together must define the left preconditoner
 * matrix P approximating the system Jacobian matrix
 *                   J = dF/du + cj*dF/du'
 * (where the DAE system is F(t,u,u') = 0), and solve the linear
 * systems P z = r.   This is done in this case by keeping only
 * the diagonal elements of the J matrix above, storing them as
 * inverses in a vector pp, when computed in PsetupHeat, for
 * subsequent use in PsolveHeat.
 *
 * In this instance, only cj and data (user data structure, with
 * pp etc.) are used from the PsetupHeat argument list.
 *
 */

int PsetupHeat(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
               realtype c_j, void *user_data)
{
  realtype *ppv, pelinv;
  sunindextype lx, ly, ixbegin, ixend, jybegin, jyend, locu, mxsub, mysub;
  int ixsub, jysub, npex, npey;
  UserData data;

  data = (UserData) user_data;

  ppv = N_VGetArrayPointer(data->pp);
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
  if (ixsub == 0) ixbegin++;
  if (ixsub == npex-1) ixend--;
  if (jysub == 0) jybegin++;
  if (jysub == npey-1) jyend--;
  pelinv = ONE/(c_j + data->coeffxy);

  /* Load the inverse of the preconditioner diagonal elements
     in loop over all the local subgrid. */

  for (ly = jybegin; ly <=jyend; ly++) {
    for (lx = ixbegin; lx <= ixend; lx++) {
      locu  = lx + ly*mxsub;
      ppv[locu] = pelinv;
    }
  }

  return(0);

}

/*
 * PsolveHeat: solve preconditioner linear system.
 * This routine multiplies the input vector rvec by the vector pp
 * containing the inverse diagonal Jacobian elements (previously
 * computed in PsetupHeat), returning the result in zvec.
 */

int PsolveHeat(realtype tt, N_Vector uu, N_Vector up,
               N_Vector rr, N_Vector rvec, N_Vector zvec,
               realtype c_j, realtype delta, void *user_data)
{
  UserData data;

  data = (UserData) user_data;

  N_VProd(data->pp, rvec, zvec);

  return(0);

}

/*
 *--------------------------------------------------------------------
 * SUPPORTING FUNCTIONS
 *--------------------------------------------------------------------
 */


/*
 * rescomm routine.  This routine performs all inter-processor
 * communication of data in u needed to calculate G.
 */

static int rescomm(N_Vector uu, N_Vector up, void *user_data)
{
  UserData data;
  realtype *uarray, *uext, buffer[2*MYSUB];
  MPI_Comm comm;
  int thispe, ixsub, jysub;
  sunindextype mxsub, mysub;
  MPI_Request request[4];

  data = (UserData) user_data;
  uarray = N_VGetArrayPointer(uu);

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

}

/*
 * reslocal routine.  Compute res = F(t, uu, up).  This routine assumes
 * that all inter-processor communication of data needed to calculate F
 * has already been done, and that this data is in the work array uext.
 */

static int reslocal(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
                    void *user_data)
{
  realtype *uext, *uuv, *upv, *resv;
  realtype termx, termy, termctr;
  sunindextype lx, ly, offsetu, offsetue, locu, locue;
  int ixsub, jysub, npex, npey;
  sunindextype mxsub, mxsub2, mysub;
  sunindextype ixbegin, ixend, jybegin, jyend;
  UserData data;

  /* Get subgrid indices, array sizes, extended work array uext. */

  data = (UserData) user_data;
  uext = data->uext;
  uuv = N_VGetArrayPointer(uu);
  upv = N_VGetArrayPointer(up);
  resv = N_VGetArrayPointer(rr);
  ixsub = data->ixsub; jysub = data->jysub;
  mxsub = data->mxsub; mxsub2 = data->mxsub + 2;
  mysub = data->mysub; npex = data->npex; npey = data->npey;

  /* Initialize all elements of rr to uu. This sets the boundary
     elements simply without indexing hassles. */

  N_VScale(ONE, uu, rr);

  /* Copy local segment of u vector into the working extended array uext.
     This completes uext prior to the computation of the rr vector.     */

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
  if (ixsub == 0) ixbegin++;
  if (ixsub == npex-1) ixend--;
  if (jysub == 0) jybegin++;
  if (jysub == npey-1) jyend--;

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

}

/*
 * Routine to send boundary data to neighboring PEs.
 */

static int BSend(MPI_Comm comm, int thispe, int ixsub, int jysub,
                 sunindextype dsizex, sunindextype dsizey, realtype uarray[])
{
  sunindextype ly, offsetu;
  realtype bufleft[MYSUB], bufright[MYSUB];

  /* If jysub > 0, send data from bottom x-line of u. */

  if (jysub != 0)
    MPI_Send(&uarray[0], (int) dsizex, MPI_SUNREALTYPE, thispe-NPEX, 0, comm);

  /* If jysub < NPEY-1, send data from top x-line of u. */

  if (jysub != NPEY-1) {
    offsetu = (MYSUB-1)*dsizex;
    MPI_Send(&uarray[offsetu], (int) dsizex, MPI_SUNREALTYPE,
             thispe+NPEX, 0, comm);
  }

  /* If ixsub > 0, send data from left y-line of u (via bufleft). */

  if (ixsub != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*dsizex;
      bufleft[ly] = uarray[offsetu];
    }
    MPI_Send(&bufleft[0], (int) dsizey, MPI_SUNREALTYPE, thispe-1, 0, comm);
  }

  /* If ixsub < NPEX-1, send data from right y-line of u (via bufright). */

  if (ixsub != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*MXSUB + (MXSUB-1);
      bufright[ly] = uarray[offsetu];
    }
    MPI_Send(&bufright[0], (int) dsizey, MPI_SUNREALTYPE, thispe+1, 0, comm);
  }

  return(0);

}

/*
 * Routine to start receiving boundary data from neighboring PEs.
 * Notes:
 *   1) buffer should be able to hold 2*MYSUB realtype entries, should be
 *      passed to both the BRecvPost and BRecvWait functions, and should not
 *      be manipulated between the two calls.
 *   2) request should have 4 entries, and should be passed in
 *      both calls also.
 */

static int BRecvPost(MPI_Comm comm, MPI_Request request[], int thispe,
                     int ixsub, int jysub, sunindextype dsizex,
                     sunindextype dsizey, realtype uext[], realtype buffer[])
{
  sunindextype offsetue;
  /* Have bufleft and bufright use the same buffer. */
  realtype *bufleft = buffer, *bufright = buffer+MYSUB;

  /* If jysub > 0, receive data for bottom x-line of uext. */
  if (jysub != 0)
    MPI_Irecv(&uext[1], (int) dsizex, MPI_SUNREALTYPE,
              thispe-NPEX, 0, comm, &request[0]);

  /* If jysub < NPEY-1, receive data for top x-line of uext. */
  if (jysub != NPEY-1) {
    offsetue = (1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&uext[offsetue], (int) dsizex, MPI_SUNREALTYPE,
              thispe+NPEX, 0, comm, &request[1]);
  }

  /* If ixsub > 0, receive data for left y-line of uext (via bufleft). */
  if (ixsub != 0) {
    MPI_Irecv(&bufleft[0], (int) dsizey, MPI_SUNREALTYPE,
              thispe-1, 0, comm, &request[2]);
  }

  /* If ixsub < NPEX-1, receive data for right y-line of uext (via bufright). */
  if (ixsub != NPEX-1) {
    MPI_Irecv(&bufright[0], (int) dsizey, MPI_SUNREALTYPE,
              thispe+1, 0, comm, &request[3]);
  }

  return(0);

}

/*
 * Routine to finish receiving boundary data from neighboring PEs.
 * Notes:
 *   1) buffer should be able to hold 2*MYSUB realtype entries, should be
 *      passed to both the BRecvPost and BRecvWait functions, and should not
 *      be manipulated between the two calls.
 *   2) request should have four entries, and should be passed in both
 *      calls also.
 */

static int BRecvWait(MPI_Request request[], int ixsub, int jysub,
                     sunindextype dsizex, realtype uext[], realtype buffer[])
{
  sunindextype ly, dsizex2, offsetue;
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

}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * InitUserData initializes the user's data block data.
 */

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

}

/*
 * SetInitialProfile sets the initial values for the problem.
 */

static int SetInitialProfile(N_Vector uu, N_Vector up,  N_Vector id,
                             N_Vector res, UserData data)
{
  int ixsub, jysub;
  sunindextype i, iloc, j, jloc, offset, loc;
  sunindextype ixbegin, ixend, jybegin, jyend;
  realtype xfact, yfact, *udata, *iddata;

  /* Initialize uu. */

  udata = N_VGetArrayPointer(uu);
  iddata = N_VGetArrayPointer(id);

  /* Set mesh spacings and subgrid indices for this PE. */
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
      udata[loc] = RCONST(16.0) * xfact * (ONE - xfact) * yfact * (ONE - yfact);
      if (i == 0 || i == MX-1 || j == 0 || j == MY-1) iddata[loc] = ZERO;
    }
  }

  /* Initialize up. */

  N_VConst(ZERO, up);    /* Initially set up = 0. */

  /* resHeat sets res to negative of ODE RHS values at interior points. */
  resHeat(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct initial up values. */
  N_VScale(-ONE, res, up);

  return(0);
}

/*
 * Print first lines of output and table heading
 */

static void PrintHeader(sunindextype Neq, realtype rtol, realtype atol)
{
  printf("\nidasHeat2D_p: Heat equation, parallel example problem for IDA\n");
  printf("                Discretized heat equation on 2D unit square.\n");
  printf("                Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("                Mesh dimensions: %d x %d", MX, MY);
  printf("        Total system size: %ld\n\n", (long int) Neq);
  printf("Subgrid dimensions: %d x %d", MXSUB, MYSUB);
  printf("        Processor array: %d x %d\n", NPEX, NPEY);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");
  printf("SUPPRESSALG = SUNTRUE to suppress local error testing on ");
  printf("all boundary components. \n");
  printf("Linear solver: SUNLinSol_SPGMR  ");
  printf("Preconditioner: diagonal elements only.\n");

  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nli   nre   nreLS    h      npe nps\n");
  printf("----------------------------------------------------------------------\n");
}

/*
 * PrintOutput: print max norm of solution and current solver statistics
 */

static void PrintOutput(int id, void *ida_mem, realtype t, N_Vector uu)
{
  realtype hused, umax;
  long int nst, nni, nje, nre, nreLS, nli, npe, nps;
  int kused, retval;

  umax = N_VMaxNorm(uu);

  if (id == 0) {

    retval = IDAGetLastOrder(ida_mem, &kused);
    check_retval(&retval, "IDAGetLastOrder", 1, id);
    retval = IDAGetNumSteps(ida_mem, &nst);
    check_retval(&retval, "IDAGetNumSteps", 1, id);
    retval = IDAGetNumNonlinSolvIters(ida_mem, &nni);
    check_retval(&retval, "IDAGetNumNonlinSolvIters", 1, id);
    retval = IDAGetNumResEvals(ida_mem, &nre);
    check_retval(&retval, "IDAGetNumResEvals", 1, id);
    retval = IDAGetLastStep(ida_mem, &hused);
    check_retval(&retval, "IDAGetLastStep", 1, id);
    retval = IDAGetNumJtimesEvals(ida_mem, &nje);
    check_retval(&retval, "IDAGetNumJtimesEvals", 1, id);
    retval = IDAGetNumLinIters(ida_mem, &nli);
    check_retval(&retval, "IDAGetNumLinIters", 1, id);
    retval = IDAGetNumLinResEvals(ida_mem, &nreLS);
    check_retval(&retval, "IDAGetNumLinResEvals", 1, id);
    retval = IDAGetNumPrecEvals(ida_mem, &npe);
    check_retval(&retval, "IDAGetPrecEvals", 1, id);
    retval = IDAGetNumPrecSolves(ida_mem, &nps);
    check_retval(&retval, "IDAGetNumPrecSolves", 1, id);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2Le  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#else
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#endif

  }
}

/*
 * Print some final integrator statistics
 */

static void PrintFinalStats(void *ida_mem)
{
  long int netf, ncfn, ncfl;

  IDAGetNumErrTestFails(ida_mem, &netf);
  IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  IDAGetNumLinConvFails(ida_mem, &ncfl);

  printf("\nError test failures            = %ld\n", netf);
  printf("Nonlinear convergence failures = %ld\n", ncfn);
  printf("Linear convergence failures    = %ld\n", ncfl);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id)
{
  int *retval;

  if (opt == 0 && returnvalue == NULL) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if retval < 0 */
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR(%d): %s() failed with retval = %d\n\n",
              id, funcname, *retval);
      return(1);
    }
  } else if (opt == 2 && returnvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  }

  return(0);
}
