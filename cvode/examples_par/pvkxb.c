/************************************************************************
 * File       : pvkxb.c                                                 *
 * Programmers: S. D. Cohen, A. C. Hindmarsh, M. R. Wittman, and        *
 *              Radu Serban  @LLNL                                      *
 * Version of : 11 July 2003                                            *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * An ODE system is generated from the following 2-species diurnal      *
 * kinetics advection-diffusion PDE system in 2 space dimensions:       *
 *                                                                      *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)    *
 *                 + Ri(c1,c2,t)      for i = 1,2,   where              *
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,       *
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,                    *
 *   Kv(y) = Kv0*exp(y/5) ,                                             *
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)        *
 * vary diurnally.   The problem is posed on the square                 *
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),                      *
 * with homogeneous Neumann boundary conditions, and for time t in      *
 *   0 <= t <= 86400 sec (1 day).                                       *
 * The PDE system is treated by central differences on a uniform        *
 * mesh, with simple polynomial initial profiles.                       *
 *                                                                      *
 * The problem is solved by CVODE on NPE processors, treated as a       *
 * rectangular process grid of size NPEX by NPEY, with NPE = NPEX*NPEY. *
 * Each processor contains a subgrid of size MXSUB by MYSUB of the      *
 * (x,y) mesh.  Thus the actual mesh sizes are MX = MXSUB*NPEX and      *
 * MY = MYSUB*NPEY, and the ODE system size is neq = 2*MX*MY.           *
 *                                                                      *
 * The solution with CVODE is done with the BDF/GMRES method (i.e.      *
 * using the CVSPGMR linear solver) and a block-diagonal matrix with    *
 * banded blocks as a preconditioner, using the CVBBDPRE module.        *
 * Each block is generated using difference quotients, with             *
 * half-bandwidths mudq = mldq = 2*MXSUB, but the retained banded       *
 * blocks have half-bandwidths mukeep = mlkeep = 2.                     *
 * A copy of the approximate Jacobian is saved and conditionally reused *
 * within the Precond routine.                                          *
 *                                                                      *
 * The problem is solved twice -- with left and right preconditioning.  *
 *                                                                      *
 * Performance data and sampled solution values are printed at selected *
 * output times, and all performance counters are printed on completion.*
 *                                                                      *
 * This version uses MPI for user routines, and the CVODE solver.       *   
 * Execute with number of processors = NPEX*NPEY (see constants below). *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h" /* definitions of realtype, integertype            */
#include "cvode.h"         /* main CVODE header file                          */
#include "iterative.h"     /* contains the enum for types of preconditioning  */
#include "cvspgmr.h"       /* use CVSPGMR linear solver                       */
#include "cvbbdpre.h"      /* band preconditioner function prototypes         */
#include "nvector_parallel.h" /* definitions of type N_Vector, macro NV_DATA_P*/
#include "sundialsmath.h"  /* contains SQR macro                              */
#include "mpi.h"           /* MPI data types and prototypes                   */


/* Problem Constants */

#define NVARS        2             /* number of species         */
#define KH           4.0e-6        /* horizontal diffusivity Kh */
#define VEL          0.001         /* advection velocity V      */
#define KV0          1.0e-8        /* coefficient in Kv(y)      */
#define Q1           1.63e-16      /* coefficients q1, q2, c3   */ 
#define Q2           4.66e-16
#define C3           3.7e16
#define A3           22.62         /* coefficient in expression for q3(t) */
#define A4           7.601         /* coefficient in expression for q4(t) */
#define C1_SCALE     1.0e6         /* coefficients in initial profiles    */
#define C2_SCALE     1.0e12

#define T0           0.0           /* initial time */
#define NOUT         12            /* number of output times */
#define TWOHR        7200.0        /* number of seconds in two hours  */
#define HALFDAY      4.32e4        /* number of seconds in a half day */
#define PI       3.1415926535898   /* pi */ 

#define XMIN          0.0          /* grid boundaries in x  */
#define XMAX         20.0           
#define YMIN         30.0          /* grid boundaries in y  */
#define YMAX         50.0

#define NPEX         2              /* no. PEs in x direction of PE array */
#define NPEY         2              /* no. PEs in y direction of PE array */
                                    /* Total no. PEs = NPEX*NPEY */
#define MXSUB        5              /* no. x points per subgrid */
#define MYSUB        5              /* no. y points per subgrid */

#define MX           (NPEX*MXSUB)   /* MX = number of x mesh points */
#define MY           (NPEY*MYSUB)   /* MY = number of y mesh points */
                                    /* Spatial mesh is MX by MY */

/* CVodeMalloc Constants */

#define RTOL    1.0e-5            /* scalar relative tolerance */
#define FLOOR   100.0             /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */


/* Type : UserData 
   contains problem constants, extended dependent variable array,
   grid constants, processor indices, MPI communicator */

typedef struct {
  realtype q4, om, dx, dy, hdco, haco, vdco;
  realtype uext[NVARS*(MXSUB+2)*(MYSUB+2)];
  integertype my_pe, isubx, isuby, nvmxsub, nvmxsub2;
  integertype Nlocal;
  MPI_Comm comm;
} *UserData;


/* Prototypes of private helper functions */

static void InitUserData(int my_pe, integertype local_N, MPI_Comm comm, UserData data);
static void SetInitialProfiles(N_Vector u, UserData data);
static void PrintOutput(void *cvode_mem, integertype my_pe, MPI_Comm comm,
                        N_Vector u, realtype t);
static void PrintFinalStats(void *cvode_mem);
static void BSend(MPI_Comm comm, integertype my_pe, integertype isubx,
                  integertype isuby, integertype dsizex, integertype dsizey,
                  realtype uarray[]);
static void BRecvPost(MPI_Comm comm, MPI_Request request[], integertype my_pe,
		      integertype isubx, integertype isuby,
		      integertype dsizex, integertype dsizey,
		      realtype uext[], realtype buffer[]);
static void BRecvWait(MPI_Request request[], integertype isubx,
		      integertype isuby, integertype dsizex, realtype uext[],
                      realtype buffer[]);

static void fucomm(realtype t, N_Vector u, void *f_data);


/* Prototype of function called by the CVODE solver */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data);


/* Prototype of functions called by the CVBBDPRE module */

static void flocal(integertype Nlocal, realtype t, N_Vector u,
                   N_Vector udot, void *f_data);

static void ucomm(integertype Nlocal, realtype t, N_Vector u, void *f_data);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  NV_Spec nvSpec;
  UserData data;
  void *cvode_mem;
  void *pdata;
  realtype abstol, reltol, t, tout;
  N_Vector u;
  int iout, flag, my_pe, npes, jpre;
  long int leniwBBDP, lenrwBBDP;
  int ngevalsBBDP;
  integertype neq, local_N, mudq, mldq, mukeep, mlkeep;
  MPI_Comm comm;

  /* Set problem size neq */
  neq = NVARS*MX*MY;

  /* Get processor number and total number of pe's */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  if (npes != NPEX*NPEY) {
    if (my_pe == 0)
      printf("\n npes=%d is not equal to NPEX*NPEY=%d\n", npes,NPEX*NPEY);
    return(1);
  }

  /* Set local length */
  local_N = NVARS*MXSUB*MYSUB;

  nvSpec = NV_SpecInit_Parallel(comm, local_N, neq, &argc, &argv);

  /* Allocate and load user data block */
  data = (UserData) malloc(sizeof *data);
  InitUserData(my_pe, local_N, comm, data);

  /* Allocate and initialize u, and set tolerances */ 
  u = N_VNew(nvSpec);
  SetInitialProfiles(u, data);
  abstol = ATOL; reltol = RTOL;

  /* 
     Call CVodeCreate to create CVODE memory:
     
     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration

     A pointer to CVODE problem memory is returned and stored in cvode_mem.
  */

  cvode_mem = CVodeCreate(BDF, NEWTON);
  if (cvode_mem == NULL) { printf("CVodeCreate failed.\n"); return(1); }

  /* Set the pointer to user-defined data */
  flag = CVodeSetFdata(cvode_mem, data);
  if (flag != SUCCESS) { printf("CVodeSetFdata failed.\n"); return(1); }

  /* 
     Call CVodeMalloc to initialize CVODE memory: 

     cvode_mem is the pointer to CVODE memory returned by CVodeCreate
     f       is the user's right hand side function in y'=f(t,y)
     T0      is the initial time
     u       is the initial dependent variable vector
     SS      specifies scalar relative and absolute tolerances
     &reltol and &abstol are pointers to the scalar tolerances
     nvSpec  is the vector specification object 
  */

  flag = CVodeMalloc(cvode_mem, f, T0, u, SS, &reltol, &abstol, nvSpec);
  if (flag != SUCCESS) { printf("CVodeMalloc failed.\n"); return(1); }

  /* Allocate preconditioner block */
  mudq = mldq = NVARS*MXSUB;
  mukeep = mlkeep = NVARS;
  pdata = CVBBDPrecAlloc(cvode_mem, local_N, mudq, mldq, 
                         mukeep, mlkeep, 0.0, flocal, ucomm);
  if (pdata == NULL) { printf("CVBBDPrecAlloc failed."); return(1); }

  /* Call CVBBDSpgmr to specify the CVODE linear solver CVSPGMR 
     with left preconditioning, the maximum Krylov dimension maxl,
     and using the CVBBDPRE preconditioner */
  flag = CVBBDSpgmr(cvode_mem, LEFT, 0, pdata);
  if (flag != SUCCESS) { printf("CVBBDSpgmr failed."); return(1); }

  /* Print heading */
  if (my_pe == 0) {
    printf("2-species diurnal advection-diffusion problem\n");
    printf("  %d by %d mesh on %d processors\n",MX,MY,npes);
    printf("  Using CVBBDPRE preconditioner module\n");
    printf("    Difference-quotient half-bandwidths are");
    printf(" mudq = %ld,  mldq = %ld\n", mudq, mldq);
    printf("    Retained band block half-bandwidths are");
    printf(" mukeep = %ld,  mlkeep = %ld", mukeep, mlkeep); 
  }

  /* Loop over jpre (= LEFT, RIGHT), and solve the problem */
  for (jpre = LEFT; jpre <= RIGHT; jpre++) {

  /* On second run, re-initialize u, CVODE, CVBBDPRE, and CVSPGMR */

  if (jpre == RIGHT) {

    SetInitialProfiles(u, data);

    flag = CVodeReInit(cvode_mem, f, T0, u, SS, &reltol, &abstol);
    if (flag != SUCCESS) { printf("CVodeReInit failed."); return(1); }

    flag = CVBBDPrecReInit(pdata, mudq, mldq, 0.0, flocal, ucomm);

    flag = CVSpgmrResetPrecType(cvode_mem, RIGHT);

    if (my_pe == 0) {
      printf("\n\n-------------------------------------------------------");
      printf("------------\n");
    }
  }

  if (my_pe == 0) {
    printf("\n\nPreconditioner type is:  jpre = %s\n\n",
           (jpre == LEFT) ? "LEFT" : "RIGHT");
  }

  /* In loop over output points, call CVode, print results, test for error */

  for (iout = 1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    PrintOutput(cvode_mem, my_pe, comm, u, t);
    if (flag != SUCCESS) {
      if (my_pe == 0) printf("CVode failed, flag =%d.\n", flag);
      break;
    }
  }

  /* Print final statistics */

  if (my_pe == 0) {
    PrintFinalStats(cvode_mem);
    CVBBDPrecGetIntWorkSpace(pdata, &leniwBBDP);
    CVBBDPrecGetRealWorkSpace(pdata, &lenrwBBDP);
    CVBBDPrecGetNumGfnEvals(pdata, &ngevalsBBDP);
    printf("In CVBBDPRE: real/integer local work space sizes = %ld, %ld\n",
           lenrwBBDP, leniwBBDP);  
    printf("             no. flocal evals. = %d\n",ngevalsBBDP);
  }

  } /* End of jpre loop */

  /* Free memory */

  N_VFree(u);
  CVBBDPrecFree(pdata);
  free(data);
  CVodeFree(cvode_mem);
  NV_SpecFree_Parallel(nvSpec);
  MPI_Finalize();

  return(0);
}


/*********************** Private Helper Functions ************************/

/* Load constants in data */

static void InitUserData(int my_pe, integertype local_N, MPI_Comm comm, UserData data)
{
  integertype isubx, isuby;

  /* Set problem constants */
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/((realtype)(MX-1));
  data->dy = (YMAX-YMIN)/((realtype)(MY-1));
  data->hdco = KH/SQR(data->dx);
  data->haco = VEL/(2.0*data->dx);
  data->vdco = (1.0/SQR(data->dy))*KV0;

  /* Set machine-related constants */
  data->comm = comm;
  data->my_pe = my_pe;
  data->Nlocal = local_N;
  /* isubx and isuby are the PE grid indices corresponding to my_pe */
  isuby = my_pe/NPEX;
  isubx = my_pe - isuby*NPEX;
  data->isubx = isubx;
  data->isuby = isuby;
  /* Set the sizes of a boundary x-line in u and uext */
  data->nvmxsub = NVARS*MXSUB;
  data->nvmxsub2 = NVARS*(MXSUB+2);

}

/* Set initial conditions in u */

static void SetInitialProfiles(N_Vector u, UserData data)
{
  integertype isubx, isuby, lx, ly, jx, jy, offset;
  realtype dx, dy, x, y, cx, cy, xmid, ymid;
  realtype *uarray;

  /* Set pointer to data array in vector u */

  uarray = NV_DATA_P(u);

  /* Get mesh spacings, and subgrid indices for this PE */

  dx = data->dx;         dy = data->dy;
  isubx = data->isubx;   isuby = data->isuby;

  /* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. */

  offset = 0;
  xmid = .5*(XMIN + XMAX);
  ymid = .5*(YMIN + YMAX);
  for (ly = 0; ly < MYSUB; ly++) {
    jy = ly + isuby*MYSUB;
    y = YMIN + jy*dy;
    cy = SQR(0.1*(y - ymid));
    cy = 1.0 - cy + 0.5*SQR(cy);
    for (lx = 0; lx < MXSUB; lx++) {
      jx = lx + isubx*MXSUB;
      x = XMIN + jx*dx;
      cx = SQR(0.1*(x - xmid));
      cx = 1.0 - cx + 0.5*SQR(cx);
      uarray[offset  ] = C1_SCALE*cx*cy; 
      uarray[offset+1] = C2_SCALE*cx*cy;
      offset = offset + 2;
    }
  }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, integertype my_pe, MPI_Comm comm, 
                        N_Vector u, realtype t)
{
  int nst, qu;
  realtype hu, *uarray, tempu[2];
  integertype npelast, i0, i1;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;
  uarray = NV_DATA_P(u);

  /* Send c1,c2 at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&uarray[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      tempu[0] = uarray[i0];
      tempu[1] = uarray[i1];
    }
  }

  /* On PE 0, receive c1,c2 at top right, then print performance data
     and sampled solution values */ 
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&tempu[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    CVodeGetNumSteps(cvode_mem, &nst);
    CVodeGetLastOrder(cvode_mem, &qu);
    CVodeGetLastStep(cvode_mem, &hu);
    printf("t = %.2e   no. steps = %d   order = %d   stepsize = %.2e\n",
           t, nst, qu, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", uarray[0], uarray[1]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
  }
}

/* Print final statistics contained in iopt */

static void PrintFinalStats(void *cvode_mem)
{
  long int lenrw, leniw ;
  long int lenrwSPGMR, leniwSPGMR;
  int nst, nfe, nsetups, nni, ncfn, netf;
  int nli, npe, nps, ncfl, nfeSPGMR;
  

  CVodeGetIntWorkSpace(cvode_mem, &leniw);
  CVodeGetRealWorkSpace(cvode_mem, &lenrw);
  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetNumRhsEvals(cvode_mem, &nfe);
  CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  CVodeGetNumErrTestFails(cvode_mem, &netf);
  CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  CVSpgmrGetIntWorkSpace(cvode_mem, &leniwSPGMR);
  CVSpgmrGetRealWorkSpace(cvode_mem, &lenrwSPGMR);
  CVSpgmrGetNumLinIters(cvode_mem, &nli);
  CVSpgmrGetNumPrecEvals(cvode_mem, &npe);
  CVSpgmrGetNumPrecSolves(cvode_mem, &nps);
  CVSpgmrGetNumConvFails(cvode_mem, &ncfl);
  CVSpgmrGetNumRhsEvals(cvode_mem, &nfeSPGMR);

  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %5ld     leniw = %5ld\n", lenrw, leniw);
  printf("llrw    = %5ld     lliw  = %5ld\n", lenrwSPGMR, leniwSPGMR);
  printf("nst     = %5d\n"                  , nst);
  printf("nfe     = %5d     nfel  = %5d\n"  , nfe, nfeSPGMR);
  printf("nni     = %5d     nli   = %5d\n"  , nni, nli);
  printf("nsetups = %5d     netf  = %5d\n"  , nsetups, netf);
  printf("npe     = %5d     nps   = %5d\n"  , npe, nps);
  printf("ncfn    = %5d     ncfl  = %5d\n\n", ncfn, ncfl); 
}
 
/* Routine to send boundary data to neighboring PEs */

static void BSend(MPI_Comm comm, integertype my_pe, integertype isubx,
                  integertype isuby, integertype dsizex, integertype dsizey,
                  realtype uarray[])
{
  int i, ly;
  integertype offsetu, offsetbuf;
  realtype bufleft[NVARS*MYSUB], bufright[NVARS*MYSUB];

  /* If isuby > 0, send data from bottom x-line of u */

  if (isuby != 0)
    MPI_Send(&uarray[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If isuby < NPEY-1, send data from top x-line of u */

  if (isuby != NPEY-1) {
    offsetu = (MYSUB-1)*dsizex;
    MPI_Send(&uarray[offsetu], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If isubx > 0, send data from left y-line of u (via bufleft) */

  if (isubx != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = ly*dsizex;
      for (i = 0; i < NVARS; i++)
        bufleft[offsetbuf+i] = uarray[offsetu+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);   
  }

  /* If isubx < NPEX-1, send data from right y-line of u (via bufright) */

  if (isubx != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = offsetbuf*MXSUB + (MXSUB-1)*NVARS;
      for (i = 0; i < NVARS; i++)
        bufright[offsetbuf+i] = uarray[offsetu+i];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe+1, 0, comm);   
  }

}
 
/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvPost(MPI_Comm comm, MPI_Request request[], integertype my_pe,
		      integertype isubx, integertype isuby,
		      integertype dsizex, integertype dsizey,
		      realtype uext[], realtype buffer[])
{
  integertype offsetue;
  /* Have bufleft and bufright use the same buffer */
  realtype *bufleft = buffer, *bufright = buffer+NVARS*MYSUB;

  /* If isuby > 0, receive data for bottom x-line of uext */
  if (isuby != 0)
    MPI_Irecv(&uext[NVARS], dsizex, PVEC_REAL_MPI_TYPE,
    					 my_pe-NPEX, 0, comm, &request[0]);

  /* If isuby < NPEY-1, receive data for top x-line of uext */
  if (isuby != NPEY-1) {
    offsetue = NVARS*(1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&uext[offsetue], dsizex, PVEC_REAL_MPI_TYPE,
                                         my_pe+NPEX, 0, comm, &request[1]);
  }

  /* If isubx > 0, receive data for left y-line of uext (via bufleft) */
  if (isubx != 0) {
    MPI_Irecv(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         my_pe-1, 0, comm, &request[2]);
  }

  /* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Irecv(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         my_pe+1, 0, comm, &request[3]);
  }

}

/* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvWait(MPI_Request request[], integertype isubx,
		      integertype isuby, integertype dsizex, realtype uext[],
                      realtype buffer[])
{
  int i, ly;
  integertype dsizex2, offsetue, offsetbuf;
  realtype *bufleft = buffer, *bufright = buffer+NVARS*MYSUB;
  MPI_Status status;

  dsizex2 = dsizex + 2*NVARS;

  /* If isuby > 0, receive data for bottom x-line of uext */
  if (isuby != 0)
    MPI_Wait(&request[0],&status);

  /* If isuby < NPEY-1, receive data for top x-line of uext */
  if (isuby != NPEY-1)
    MPI_Wait(&request[1],&status);

  /* If isubx > 0, receive data for left y-line of uext (via bufleft) */
  if (isubx != 0) {
    MPI_Wait(&request[2],&status);

    /* Copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetue = (ly+1)*dsizex2;
      for (i = 0; i < NVARS; i++)
        uext[offsetue+i] = bufleft[offsetbuf+i];
    }
  }

  /* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Wait(&request[3],&status);

    /* Copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetue = (ly+2)*dsizex2 - NVARS;
      for (i = 0; i < NVARS; i++)
	uext[offsetue+i] = bufright[offsetbuf+i];
    }
  }

}

/* fucomm routine.  This routine performs all inter-processor
   communication of data in u needed to calculate f.         */

static void fucomm(realtype t, N_Vector u, void *f_data)
{

  UserData data;
  realtype *uarray, *uext, buffer[2*NVARS*MYSUB];
  MPI_Comm comm;
  integertype my_pe, isubx, isuby, nvmxsub, nvmysub;
  MPI_Request request[4];

  data = (UserData) f_data;
  uarray = NV_DATA_P(u);


  /* Get comm, my_pe, subgrid indices, data sizes, extended array uext */

  comm = data->comm;  my_pe = data->my_pe;
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;
  nvmysub = NVARS*MYSUB;
  uext = data->uext;

  /* Start receiving boundary data from neighboring PEs */

  BRecvPost(comm, request, my_pe, isubx, isuby, nvmxsub, nvmysub, uext, buffer);

  /* Send data from boundary of local grid to neighboring PEs */

  BSend(comm, my_pe, isubx, isuby, nvmxsub, nvmysub, uarray);

  /* Finish receiving boundary data from neighboring PEs */

  BRecvWait(request, isubx, isuby, nvmxsub, uext, buffer);

}


/***************** Function called by the CVODE solver *******************/

/* f routine.  Evaluate f(t,y).  First call fucomm to do communication of 
   subgrid boundary data into uext.  Then calculate f by a call to flocal. */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data)
{
  
  UserData data;

  data = (UserData) f_data;

  /* Call fucomm to do inter-processor communication */

  fucomm (t, u, f_data);

  /* Call flocal to calculate all right-hand sides */

  flocal (data->Nlocal, t, u, udot, f_data);

}


/***************** Functions called by the CVBBDPRE module ****************/

/* flocal routine.  Compute f(t,y).  This routine assumes that all
   inter-processor communication of data needed to calculate f has already
   been done, and this data is in the work array uext.                    */

static void flocal(integertype Nlocal, realtype t, N_Vector u,
                   N_Vector udot, void *f_data)
{
  realtype *uext;
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  realtype q4coef, dely, verdco, hordco, horaco;
  int i, lx, ly, jx, jy;
  integertype isubx, isuby, nvmxsub, nvmxsub2, offsetu, offsetue;
  UserData data;
  realtype *uarray, *duarray;

  uarray = NV_DATA_P(u);
  duarray = NV_DATA_P(udot);

  /* Get subgrid indices, array sizes, extended work array uext */

  data = (UserData) f_data;
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub; nvmxsub2 = data->nvmxsub2;
  uext = data->uext;

  /* Copy local segment of u vector into the working extended array uext */

  offsetu = 0;
  offsetue = nvmxsub2 + NVARS;
  for (ly = 0; ly < MYSUB; ly++) {
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = uarray[offsetu+i];
    offsetu = offsetu + nvmxsub;
    offsetue = offsetue + nvmxsub2;
  }

  /* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of u to uext */

  /* If isuby = 0, copy x-line 2 of u to uext */
  if (isuby == 0) {
    for (i = 0; i < nvmxsub; i++) uext[NVARS+i] = uarray[nvmxsub+i];
  }

  /* If isuby = NPEY-1, copy x-line MYSUB-1 of u to uext */
  if (isuby == NPEY-1) {
    offsetu = (MYSUB-2)*nvmxsub;
    offsetue = (MYSUB+1)*nvmxsub2 + NVARS;
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = uarray[offsetu+i];
  }

  /* If isubx = 0, copy y-line 2 of u to uext */
  if (isubx == 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*nvmxsub + NVARS;
      offsetue = (ly+1)*nvmxsub2;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = uarray[offsetu+i];
    }
  }

  /* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext */
  if (isubx == NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = (ly+1)*nvmxsub - 2*NVARS;
      offsetue = (ly+2)*nvmxsub2 - NVARS;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = uarray[offsetu+i];
    }
  }

  /* Make local copies of problem variables, for efficiency */

  dely = data->dy;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Set diurnal rate coefficients as functions of t, and save q4 in 
  data block for use by preconditioner evaluation routine            */

  s = sin((data->om)*t);
  if (s > 0.0) {
    q3 = exp(-A3/s);
    q4coef = exp(-A4/s);
  } else {
    q3 = 0.0;
    q4coef = 0.0;
  }
  data->q4 = q4coef;


  /* Loop over all grid points in local subgrid */

  for (ly = 0; ly < MYSUB; ly++) {

    jy = ly + isuby*MYSUB;

    /* Set vertical diffusion coefficients at jy +- 1/2 */

    ydn = YMIN + (jy - .5)*dely;
    yup = ydn + dely;
    cydn = verdco*exp(0.2*ydn);
    cyup = verdco*exp(0.2*yup);
    for (lx = 0; lx < MXSUB; lx++) {

      jx = lx + isubx*MXSUB;

      /* Extract c1 and c2, and set kinetic rate terms */

      offsetue = (lx+1)*NVARS + (ly+1)*nvmxsub2;
      c1 = uext[offsetue];
      c2 = uext[offsetue+1];
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms */

      c1dn = uext[offsetue-nvmxsub2];
      c2dn = uext[offsetue-nvmxsub2+1];
      c1up = uext[offsetue+nvmxsub2];
      c2up = uext[offsetue+nvmxsub2+1];
      vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms */

      c1lt = uext[offsetue-2];
      c2lt = uext[offsetue-1];
      c1rt = uext[offsetue+2];
      c2rt = uext[offsetue+3];
      hord1 = hordco*(c1rt - 2.0*c1 + c1lt);
      hord2 = hordco*(c2rt - 2.0*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);

      /* Load all terms into duarray */

      offsetu = lx*NVARS + ly*nvmxsub;
      duarray[offsetu]   = vertd1 + hord1 + horad1 + rkin1; 
      duarray[offsetu+1] = vertd2 + hord2 + horad2 + rkin2;
    }
  }

}

/* ucomm routine.  This routine is empty, as communication needed
   to evaluate flocal has been done in prior call to f */

static void ucomm(integertype Nlocal, realtype t, N_Vector u, void *f_data)
{
}
