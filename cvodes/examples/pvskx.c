/************************************************************************
 *                                                                      *
 * File: pvkx.c                                                         *
 * Programmers: S. D. Cohen, A. C. Hindmarsh, Radu Serban, and          *
 *              M. R. Wittman @ LLNL                                    *
 * Version of 14 November 2001                                          *
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
 * The problem is solved by PVODE on NPE processors, treated as a       *
 * rectangular process grid of size NPEX by NPEY, with NPE = NPEX*NPEY. *
 * Each processor contains a subgrid of size MXSUB by MYSUB of the      *
 * (x,y) mesh.  Thus the actual mesh sizes are MX = MXSUB*NPEX and      *
 * MY = MYSUB*NPEY, and the ODE system size is neq = 2*MX*MY.           *
 *                                                                      *
 * The solution with PVODE is done with the BDF/GMRES method (i.e.      *
 * using the CVSPGMR linear solver) and the block-diagonal part of the  *
 * Newton matrix as a left preconditioner. A copy of the block-diagonal *
 * part of the Jacobian is saved and conditionally reused within the    *
 * Precond routine.                                                     *
 *                                                                      *
 * Performance data and sampled solution values are printed at selected *
 * output times, and all performance counters are printed on completion.*
 *                                                                      *
 * Optionally, PVODES can compute sensitivities with respect to the     *
 * problem parameters q1 and q2.                                        *
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and       *
 * STAGGERED1) can be used and sensitivities may be included in the     *
 * error test or not (error control set on FULL or PARTIAL,             *
 * respectively).                                                       *
 *                                                                      *
 * Execution:                                                           *
 *                                                                      *
 * NOTE: This version uses MPI for user routines, and the MPI_PVODE     *
 *       solver. In what follows, N is the number of processors,        *
 *       N = NPEX*NPEY (see constants below) and it is assumed that     *
 *       the MPI script mpirun is used to run a paralles application.   *
 * If no sensitivities are desired:                                     *
 *    % mpirun -np N pvskx -nosensi                                     *
 * If sensitivities are to be computed:                                 *
 *    % mpirun -np N pvskx -sensi sensi_meth err_con                    *
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of    *
 * {full, partial}.                                                     *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"   /* definitions of real, integer, boole, TRUE,FALSE */
#include "cvodes.h"     /* main CVODE header file                          */
#include "iterativ.h"   /* contains the enum for types of preconditioning  */
#include "cvsspgmr.h"   /* use CVSPGMR linear solver each internal step    */
#include "smalldense.h" /* use generic DENSE solver in preconditioning     */
#include "nvector.h"    /* definitions of type N_Vector, macro N_VDATA     */
#include "llnlmath.h"   /* contains SQR macro                              */
#include "mpi.h"


/* Problem Constants */

#define NVARS        2             /* number of species                    */
#define C1_SCALE     1.0e6         /* coefficients in initial profiles     */
#define C2_SCALE     1.0e12

#define T0           0.0           /* initial time                         */
#define NOUT         12            /* number of output times               */
#define TWOHR        7200.0        /* number of seconds in two hours       */
#define HALFDAY      4.32e4        /* number of seconds in a half day      */
#define PI       3.1415926535898   /* pi                                   */ 

#define XMIN          0.0          /* grid boundaries in x                 */
#define XMAX         20.0           
#define YMIN         30.0          /* grid boundaries in y                 */
#define YMAX         50.0

#define NPEX         3              /* no. PEs in x direction of PE array  */
#define NPEY         3              /* no. PEs in y direction of PE array  */
                                    /* Total no. PEs = NPEX*NPEY           */
#define MXSUB        5              /* no. x points per subgrid            */
#define MYSUB        5              /* no. y points per subgrid            */

#define MX           (NPEX*MXSUB)   /* MX = number of x mesh points        */
#define MY           (NPEY*MYSUB)   /* MY = number of y mesh points        */
                                    /* Spatial mesh is MX by MY            */
/* CVodeMalloc Constants */

#define RTOL    1.0e-5            /* scalar relative tolerance             */
#define FLOOR   100.0             /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance             */

/* Sensitivity constants */
#define NP    8                   /* number of problem parameters          */
#define NS    2                   /* number of sensitivities               */

#define ZERO  RCONST(0.0)


/* User-defined matrix accessor macro: IJth */

/* IJth is defined in order to write code which indexes into small dense
   matrices with a (row,column) pair, where 1 <= row,column <= NVARS.   

   IJth(a,i,j) references the (i,j)th entry of the small matrix real **a,
   where 1 <= i,j <= NVARS. The small matrix routines in dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */

#define IJth(a,i,j)        (a[j-1][i-1])


/* Type : UserData 
   contains problem constants, preconditioner blocks, pivot arrays, 
   grid constants, and processor indices */

typedef struct {
  real *p;
  real q4, om, dx, dy, hdco, haco, vdco;
  real uext[NVARS*(MXSUB+2)*(MYSUB+2)];
  integer my_pe, isubx, isuby, nvmxsub, nvmxsub2;
  MPI_Comm comm;
} *UserData;

typedef struct {
  void *f_data;
  real **P[MXSUB][MYSUB], **Jbd[MXSUB][MYSUB];
  integer *pivot[MXSUB][MYSUB];
} *PreconData;


/* Private Helper Functions */

static void WrongArgs(integer my_pe, char *argv[]);
static PreconData AllocPreconData(UserData data);
static void InitUserData(integer my_pe, MPI_Comm comm, UserData data);
static void FreePreconData(PreconData pdata);
static void SetInitialProfiles(N_Vector u, UserData data);
static void PrintOutput(integer my_pe, MPI_Comm comm, long int iopt[],
                        real ropt[], real t, N_Vector u);
static void PrintOutputS(integer my_pe, MPI_Comm comm, N_Vector *uS);
static void PrintFinalStats(boole sensi, int sensi_meth, int err_con, long int iopt[]);
static void BSend(MPI_Comm comm, integer my_pe, integer isubx, integer isuby,
                  integer dsizex, integer dsizey, real udata[]);
static void BRecvPost(MPI_Comm comm, MPI_Request request[], integer my_pe,
                      integer isubx, integer isuby,
                      integer dsizex, integer dsizey,
                      real uext[], real buffer[]);
static void BRecvWait(MPI_Request request[], integer isubx, integer isuby,
                      integer dsizex, real uext[], real buffer[]);
static void ucomm(integer N, real t, N_Vector u, UserData data);
static void fcalc(integer N, real t, real udata[], real dudata[], UserData data);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data);

static int Precond(integer N, real tn, N_Vector u, N_Vector fu, boole jok,
                   boole *jcurPtr, real gamma, N_Vector ewt, real h,
                   real uround, long int *nfePtr, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int PSolve(integer N, real tn, N_Vector u, N_Vector fu, N_Vector vtemp,
                  real gamma, N_Vector ewt, real delta, long int *nfePtr,
                  N_Vector r, int lr, void *P_data, N_Vector z);


/***************************** Main Program ******************************/

main(int argc, char *argv[])
{
  real abstol, reltol, t, tout, ropt[OPT_SIZE];
  long int iopt[OPT_SIZE];
  N_Vector u;
  UserData data;
  PreconData predata;
  void *cvode_mem;
  int iout, flag, my_pe, npes;
  integer neq, local_N;
  machEnvType machEnv;
  MPI_Comm comm;

  real *pbar, rhomax;
  integer is, *plist;
  N_Vector *uS;
  boole sensi;
  int sensi_meth, err_con, ifS;

  /* Set problem size neq */
  neq = NVARS*MX*MY;

  /* Get processor number and total number of pe's */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  /* Process arguments */
  if (argc < 2)
    WrongArgs(my_pe,argv);

  if (strcmp(argv[1],"-nosensi") == 0)
    sensi = FALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    sensi = TRUE;
  else
    WrongArgs(my_pe,argv);

  if (sensi) {

    if (argc != 4)
      WrongArgs(my_pe,argv);

    if (strcmp(argv[2],"sim") == 0)
      sensi_meth = SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      sensi_meth = STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      sensi_meth = STAGGERED1;
    else 
      WrongArgs(my_pe,argv);

    if (strcmp(argv[3],"full") == 0)
      err_con = FULL;
    else if (strcmp(argv[3],"partial") == 0)
      err_con = PARTIAL;
    else
      WrongArgs(my_pe,argv);

  }

  if (npes != NPEX*NPEY) {
    if (my_pe == 0)
      printf("\n npes=%d is not equal to NPEX*NPEY=%d\n", npes,NPEX*NPEY);
    return(1);
  }

  /* Set local length */
  local_N = NVARS*MXSUB*MYSUB;

  /* Allocate and load user data block; allocate preconditioner block */
  data = (UserData) malloc(sizeof *data);
  data->p = (real *) malloc(NP*sizeof(real));
  InitUserData(my_pe, comm, data);
  predata = AllocPreconData (data);

  /* Set machEnv block */
  machEnv = PVecInitMPI(comm, local_N, neq, &argc, &argv);
  if (machEnv == NULL) return(1);

  /* Allocate u, and set initial values and tolerances */ 
  u = N_VNew(neq, machEnv);
  SetInitialProfiles(u, data);
  abstol = ATOL; reltol = RTOL;

  for (is=0; is<OPT_SIZE; is++) {
    iopt[is] = 0;
    ropt[is] = 0.0;
  }
  iopt[MXSTEP] = 2000;

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(neq, f, T0, u, BDF, NEWTON, SS, &reltol,
                          &abstol, data, NULL, TRUE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) {
    printf("CVodeMalloc failed.");
    return(1);
  }

  /* CVSPGMR */
  flag = CVSpgmr(cvode_mem, LEFT, MODIFIED_GS, 0, 0.0, Precond, PSolve, predata, NULL, NULL);
  if (flag != SUCCESS) { printf("CVSpgmr failed.\n"); return(1); }

  /* SENSITIVTY */
  if(sensi) {
    pbar = (real *) malloc(NP*sizeof(real));
    for(is=0; is<NP; is++) pbar[is] = data->p[is];
    plist = (integer *) malloc(NS * sizeof(integer));
    for(is=0; is<NS; is++) plist[is] = is+1;

    uS = N_VNew_S(NS,neq,machEnv);
    for(is=0;is<NS;is++)
      N_VConst(ZERO,uS[is]);

    rhomax = ZERO;

    ifS = ALLSENS;
    if(sensi_meth==STAGGERED1) ifS = ONESENS;

    flag = CVodeSensMalloc(cvode_mem,NS,sensi_meth,data->p,pbar,plist,
                           ifS,NULL,err_con,rhomax,uS,NULL,NULL);
    if (flag != SUCCESS) {
      if (my_pe == 0) printf("CVodeSensMalloc failed, flag=%d\n",flag);
      return(1);
    }

  }

  if (my_pe == 0) {
    printf("\n2-species diurnal advection-diffusion problem\n\n");
    printf("========================================================================\n");
    printf("     T     Q       H      NST                    Bottom left  Top right \n");
    printf("========================================================================\n");
  }

  /* In loop over output points, call CVode, print results, test for error */
  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    if (flag != SUCCESS) {
      if (my_pe == 0) printf("CVode failed, flag=%d.\n", flag);
      break;
    }
    PrintOutput(my_pe, comm, iopt, ropt, t, u);
    if (sensi) {
      flag = CVodeSensExtract(cvode_mem, t, uS);
      if (flag != SUCCESS) { 
        printf("CVodeSensExtract failed, flag=%d.\n", flag); 
        break; 
      }
      PrintOutputS(my_pe, comm, uS);
    }

    if (my_pe == 0)
      printf("------------------------------------------------------------------------\n");

  }

  /* Print final statistics */  
  if (my_pe == 0) PrintFinalStats(sensi,sensi_meth,err_con,iopt);

  /* Free memory */
  N_VFree(u);
  if(sensi) 
    N_VFree_S(NS, uS);
  free(data->p);  
  free(data);
  FreePreconData(predata);
  CVodeFree(cvode_mem);
  PVecFreeMPI(machEnv);
  MPI_Finalize();

  return(0);
}


/*********************** Private Helper Functions ************************/

/* ======================================================================= */
/* Exit if arguments are incorrect */

static void WrongArgs(integer my_pe, char *argv[])
{
  if (my_pe == 0) {
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",argv[0]);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = full or partial\n");
  }  
  MPI_Finalize();
  exit(0);
}


/* ======================================================================= */
/* Allocate memory for data structure of type UserData */

static PreconData AllocPreconData(UserData fdata)
{
  int lx, ly;
  PreconData pdata;

  pdata = (PreconData) malloc(sizeof *pdata);

  pdata->f_data = fdata;

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      (pdata->P)[lx][ly] = denalloc(NVARS);
      (pdata->Jbd)[lx][ly] = denalloc(NVARS);
      (pdata->pivot)[lx][ly] = denallocpiv(NVARS);
    }
  }

  return(pdata);
}

/* ======================================================================= */
/* Load constants in data */

static void InitUserData(integer my_pe, MPI_Comm comm, UserData data)
{
  integer isubx, isuby;
  real Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Set problem parameters */
  Q1 = 1.63e-16; /* Q1  coefficients q1, q2, c3             */
  Q2 = 4.66e-16; /* Q2                                      */
  C3 = 3.7e16;   /* C3                                      */
  A3 = 22.62;    /* A3  coefficient in expression for q3(t) */
  A4 = 7.601;    /* A4  coefficient in expression for q4(t) */
  KH = 4.0e-6;   /* KH  horizontal diffusivity Kh           */ 
  VEL = 0.001;   /* VEL advection velocity V                */
  KV0 = 1.0e-8;  /* KV0 coefficient in Kv(z)                */ 

  /* Set problem constants */
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/((real)(MX-1));
  data->dy = (YMAX-YMIN)/((real)(MY-1));
  data->hdco = KH/SQR(data->dx);
  data->haco = VEL/(2.0*data->dx);
  data->vdco = (1.0/SQR(data->dy))*KV0;

  /* Problem parameters */
  data->p[0] = Q1;
  data->p[1] = Q2;
  data->p[2] = C3;
  data->p[3] = A3;
  data->p[4] = A4;
  data->p[5] = KH;
  data->p[6] = VEL;
  data->p[7] = KV0;

  /* Set machine-related constants */
  data->comm = comm;
  data->my_pe = my_pe;
  /* isubx and isuby are the PE grid indices corresponding to my_pe */
  isuby = my_pe/NPEX;
  isubx = my_pe - isuby*NPEX;
  data->isubx = isubx;
  data->isuby = isuby;
  /* Set the sizes of a boundary x-line in u and uext */
  data->nvmxsub = NVARS*MXSUB;
  data->nvmxsub2 = NVARS*(MXSUB+2);

}

/* ======================================================================= */
/* Free data memory */

static void FreePreconData(PreconData pdata)
{
  int lx, ly;

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      denfree((pdata->P)[lx][ly]);
      denfree((pdata->Jbd)[lx][ly]);
      denfreepiv((pdata->pivot)[lx][ly]);
    }
  }

  free(pdata);
}

/* ======================================================================= */
/* Set initial conditions in u */

static void SetInitialProfiles(N_Vector u, UserData data)
{
  integer isubx, isuby, lx, ly, jx, jy, offset;
  real dx, dy, x, y, cx, cy, xmid, ymid;
  real *udata;

  /* Set pointer to data array in vector u */

  udata = N_VDATA(u);

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
      udata[offset  ] = C1_SCALE*cx*cy; 
      udata[offset+1] = C2_SCALE*cx*cy;
      offset = offset + 2;
    }
  }
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(integer my_pe, MPI_Comm comm, long int iopt[], 
                        real ropt[], real t, N_Vector u)
{
  real *udata, tempu[2];
  integer npelast, i0, i1;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;
  udata = N_VDATA(u);

  /* Send c at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&udata[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      tempu[0] = udata[i0];
      tempu[1] = udata[i1];
    }
  }

  /* On PE 0, receive c at top right, then print performance data
     and sampled solution values */ 
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&tempu[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);

    printf("%8.3e %2d  %8.3e %5ld\n", t,iopt[QU],ropt[HU],iopt[NST]);
    printf("                                Solution       ");
    printf("%12.4e %12.4e \n", udata[0], tempu[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", udata[1], tempu[1]);

  }
}

/* ======================================================================= */
/* Print sampled sensitivity values */

static void PrintOutputS(integer my_pe, MPI_Comm comm, N_Vector *uS)
{
  real *sdata, temps[2];
  integer npelast, i0, i1;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;

  sdata = N_VDATA(uS[0]);
  /* Send s1 at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&sdata[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      temps[0] = sdata[i0];
      temps[1] = sdata[i1];
    }
  }
  /* On PE 0, receive s1 at top right, then print sampled sensitivity values */ 
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&temps[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    printf("                                ----------------------------------------\n");
    printf("                                Sensitivity 1  ");
    printf("%12.4e %12.4e \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", sdata[1], temps[1]);
  }

  sdata = N_VDATA(uS[1]);
  /* Send s2 at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&sdata[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      temps[0] = sdata[i0];
      temps[1] = sdata[i1];
    }
  }
  /* On PE 0, receive s2 at top right, then print sampled sensitivity values */ 
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&temps[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    printf("                                ----------------------------------------\n");
    printf("                                Sensitivity 2  ");
    printf("%12.4e %12.4e \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", sdata[1], temps[1]);
  }

}

/* ======================================================================= */
/* Print final statistics contained in iopt */

static void PrintFinalStats(boole sensi, int sensi_meth, int err_con, long int iopt[])
{

  printf("\n\n========================================================");
  printf("\nFinal Statistics");
  printf("\nSensitivity: ");

  if(sensi) {
    printf("YES ");
    if(sensi_meth == SIMULTANEOUS)   
      printf("( SIMULTANEOUS +");
    else 
      if(sensi_meth == STAGGERED) printf("( STAGGERED +");
      else                        printf("( STAGGERED1 +");                      
    if(err_con == FULL) printf(" FULL ERROR CONTROL )");
    else                printf(" PARTIAL ERROR CONTROL )");
  } else {
    printf("NO");
  }

  printf("\n\n");
  /*
  printf("lenrw   = %5ld    leniw = %5ld\n", iopt[LENRW], iopt[LENIW]);
  printf("llrw    = %5ld    lliw  = %5ld\n", iopt[SPGMR_LRW], iopt[SPGMR_LIW]);
  */
  printf("nst     = %5ld                \n\n", iopt[NST]);
  printf("nfe     = %5ld    nfSe  = %5ld  \n", iopt[NFE],  iopt[NFSE]);
  printf("nni     = %5ld    nniS  = %5ld  \n", iopt[NNI],  iopt[NNIS]);
  printf("ncfn    = %5ld    ncfnS = %5ld  \n", iopt[NCFN], iopt[NCFNS]);
  printf("netf    = %5ld    netfS = %5ld\n\n", iopt[NETF], iopt[NETFS]);
  printf("nsetups = %5ld                  \n", iopt[NSETUPS]);
  printf("nli     = %5ld    ncfl  = %5ld  \n", iopt[SPGMR_NLI], iopt[SPGMR_NCFL]);
  printf("npe     = %5ld    nps   = %5ld  \n", iopt[SPGMR_NPE], iopt[SPGMR_NPS]);

  printf("========================================================\n");

}

/* ======================================================================= */
/* Routine to send boundary data to neighboring PEs */

static void BSend(MPI_Comm comm, integer my_pe, integer isubx, integer isuby,
                  integer dsizex, integer dsizey, real udata[])
{
  int i, ly;
  integer offsetu, offsetbuf;
  real bufleft[NVARS*MYSUB], bufright[NVARS*MYSUB];

  /* If isuby > 0, send data from bottom x-line of u */

  if (isuby != 0)
    MPI_Send(&udata[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If isuby < NPEY-1, send data from top x-line of u */

  if (isuby != NPEY-1) {
    offsetu = (MYSUB-1)*dsizex;
    MPI_Send(&udata[offsetu], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If isubx > 0, send data from left y-line of u (via bufleft) */

  if (isubx != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = ly*dsizex;
      for (i = 0; i < NVARS; i++)
        bufleft[offsetbuf+i] = udata[offsetu+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);   
  }

  /* If isubx < NPEX-1, send data from right y-line of u (via bufright) */

  if (isubx != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = offsetbuf*MXSUB + (MXSUB-1)*NVARS;
      for (i = 0; i < NVARS; i++)
        bufright[offsetbuf+i] = udata[offsetu+i];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe+1, 0, comm);   
  }

}
 
/* ======================================================================= */
/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB real entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvPost(MPI_Comm comm, MPI_Request request[], integer my_pe,
                      integer isubx, integer isuby,
                      integer dsizex, integer dsizey,
                      real uext[], real buffer[])
{
  integer offsetue;
  /* Have bufleft and bufright use the same buffer */
  real *bufleft = buffer, *bufright = buffer+NVARS*MYSUB;

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

/* ======================================================================= */
/* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB real entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvWait(MPI_Request request[], integer isubx, integer isuby,
                      integer dsizex, real uext[], real buffer[])
{
  int i, ly;
  integer dsizex2, offsetue, offsetbuf;
  real *bufleft = buffer, *bufright = buffer+NVARS*MYSUB;
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

/* ======================================================================= */
/* ucomm routine.  This routine performs all communication 
   between processors of data needed to calculate f. */

static void ucomm(integer N, real t, N_Vector u, UserData data)
{

  real *udata, *uext, buffer[2*NVARS*MYSUB];
  MPI_Comm comm;
  integer my_pe, isubx, isuby, nvmxsub, nvmysub;
  MPI_Request request[4];

  udata = N_VDATA(u);


  /* Get comm, my_pe, subgrid indices, data sizes, extended array uext */

  comm = data->comm;  my_pe = data->my_pe;
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;
  nvmysub = NVARS*MYSUB;
  uext = data->uext;

  /* Start receiving boundary data from neighboring PEs */

  BRecvPost(comm, request, my_pe, isubx, isuby, nvmxsub, nvmysub, uext, buffer);

  /* Send data from boundary of local grid to neighboring PEs */

  BSend(comm, my_pe, isubx, isuby, nvmxsub, nvmysub, udata);

  /* Finish receiving boundary data from neighboring PEs */

  BRecvWait(request, isubx, isuby, nvmxsub, uext, buffer);



}

/* ======================================================================= */
/* fcalc routine. Compute f(t,y).  This routine assumes that communication 
   between processors of data needed to calculate f has already been done,
   and this data is in the work array uext. */

static void fcalc(integer N, real t, real udata[], real dudata[], UserData data)
{
  real *uext;
  real q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  real c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2;
  real qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  real q4coef, dely, verdco, hordco, horaco;
  int i, lx, ly, jx, jy;
  integer isubx, isuby, nvmxsub, nvmxsub2, offsetu, offsetue;
  real Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Get subgrid indices, data sizes, extended work array uext */

  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub; nvmxsub2 = data->nvmxsub2;
  uext = data->uext;

  /* Load problem coefficients and parameters */

  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];
  A3 = data->p[3];
  A4 = data->p[4];
  KH = data->p[5];
  VEL = data->p[6];
  KV0 = data->p[7];

  /* Copy local segment of u vector into the working extended array uext */

  offsetu = 0;
  offsetue = nvmxsub2 + NVARS;
  for (ly = 0; ly < MYSUB; ly++) {
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = udata[offsetu+i];
    offsetu = offsetu + nvmxsub;
    offsetue = offsetue + nvmxsub2;
  }

  /* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of u to uext */

  /* If isuby = 0, copy x-line 2 of u to uext */
  if (isuby == 0) {
    for (i = 0; i < nvmxsub; i++) uext[NVARS+i] = udata[nvmxsub+i];
  }

  /* If isuby = NPEY-1, copy x-line MYSUB-1 of u to uext */
  if (isuby == NPEY-1) {
    offsetu = (MYSUB-2)*nvmxsub;
    offsetue = (MYSUB+1)*nvmxsub2 + NVARS;
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = udata[offsetu+i];
  }

  /* If isubx = 0, copy y-line 2 of u to uext */
  if (isubx == 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*nvmxsub + NVARS;
      offsetue = (ly+1)*nvmxsub2;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = udata[offsetu+i];
    }
  }

  /* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext */
  if (isubx == NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = (ly+1)*nvmxsub - 2*NVARS;
      offsetue = (ly+2)*nvmxsub2 - NVARS;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = udata[offsetu+i];
    }
  }

  /* Make local copies of problem variables, for efficiency */

  dely = data->dy;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Set diurnal rate coefficients as functions of t, and save q4 in 
  data block for use by preconditioner evaluation routine */

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

      /* Load all terms into dudata */

      offsetu = lx*NVARS + ly*nvmxsub;
      dudata[offsetu]   = vertd1 + hord1 + horad1 + rkin1; 
      dudata[offsetu+1] = vertd2 + hord2 + horad2 + rkin2;
    }
  }

}


/***************** Functions Called by the CVODE Solver ******************/

/* ======================================================================= */
/* f routine.  Evaluate f(t,y).  First call ucomm to do communication of 
   subgrid boundary data into uext.  Then calculate f by a call to fcalc. */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data)
{
  real *udata, *dudata;
  UserData data;

  udata = N_VDATA(u);
  dudata = N_VDATA(udot);
  data = (UserData) f_data;


  /* Call ucomm to do inter-processor communicaiton */

  ucomm (N, t, u, data);

  /* Call fcalc to calculate all right-hand sides */

  fcalc (N, t, udata, dudata, data);

}


/* ======================================================================= */
/* Preconditioner setup routine. Generate and preprocess P. */

static int Precond(integer N, real tn, N_Vector u, N_Vector fu, boole jok,
                   boole *jcurPtr, real gamma, N_Vector ewt, real h,
                   real uround, long int *nfePtr, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  real c1, c2, cydn, cyup, diag, ydn, yup, q4coef, dely, verdco, hordco;
  real **(*P)[MYSUB], **(*Jbd)[MYSUB];
  integer nvmxsub, *(*pivot)[MYSUB], ier, offset;
  int lx, ly, jx, jy, isubx, isuby;
  real *udata, **a, **j;
  PreconData predata;
  UserData data;
  real Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Make local copies of pointers in P_data, pointer to u's data,
     and PE index pair */

  predata = (PreconData) P_data;
  data = (UserData) (predata->f_data);
  P = predata->P;
  Jbd = predata->Jbd;
  pivot = predata->pivot;
  udata = N_VDATA(u);
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;

  /* Load problem coefficients and parameters */
  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];
  A3 = data->p[3];
  A4 = data->p[4];
  KH = data->p[5];
  VEL = data->p[6];
  KV0 = data->p[7];

  if (jok) {

  /* jok = TRUE: Copy Jbd to P */

    for (ly = 0; ly < MYSUB; ly++)
      for (lx = 0; lx < MXSUB; lx++)
        dencopy(Jbd[lx][ly], P[lx][ly], NVARS);

  *jcurPtr = FALSE;

  }

  else {

  /* jok = FALSE: Generate Jbd from scratch and copy to P */

  /* Make local copies of problem variables, for efficiency */

  q4coef = data->q4;
  dely = data->dy;
  verdco = data->vdco;
  hordco  = data->hdco;

  /* Compute 2x2 diagonal Jacobian blocks (using q4 values 
     computed on the last f call).  Load into P. */

    for (ly = 0; ly < MYSUB; ly++) {
      jy = ly + isuby*MYSUB;
      ydn = YMIN + (jy - .5)*dely;
      yup = ydn + dely;
      cydn = verdco*exp(0.2*ydn);
      cyup = verdco*exp(0.2*yup);
      diag = -(cydn + cyup + 2.0*hordco);
      for (lx = 0; lx < MXSUB; lx++) {
        jx = lx + isubx*MXSUB;
        offset = lx*NVARS + ly*nvmxsub;
        c1 = udata[offset];
        c2 = udata[offset+1];
        j = Jbd[lx][ly];
        a = P[lx][ly];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + q4coef;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4coef) + diag;
        dencopy(j, a, NVARS);
      }
    }

  *jcurPtr = TRUE;

  }

  /* Scale by -gamma */

    for (ly = 0; ly < MYSUB; ly++)
      for (lx = 0; lx < MXSUB; lx++)
        denscale(-gamma, P[lx][ly], NVARS);

  /* Add identity matrix and do LU decompositions on blocks in place */

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      denaddI(P[lx][ly], NVARS);
      ier = gefa(P[lx][ly], NVARS, pivot[lx][ly]);
      if (ier != 0) return(1);
    }
  }

  return(0);
}


/* ======================================================================= */
/* Preconditioner solve routine */

static int PSolve(integer N, real tn, N_Vector u, N_Vector fu, N_Vector vtemp,
                  real gamma, N_Vector ewt, real delta, long int *nfePtr,
                  N_Vector r, int lr, void *P_data, N_Vector z)
{
  real **(*P)[MYSUB];
  integer nvmxsub, *(*pivot)[MYSUB];
  int lx, ly;
  real *zdata, *v;
  PreconData predata;
  UserData data;

  /* Extract the P and pivot arrays from P_data */

  predata = (PreconData) P_data;
  data = (UserData) (predata->f_data);
  P = predata->P;
  pivot = predata->pivot;

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z.
     First copy vector r to z. */

  N_VScale(1.0, r, z);

  nvmxsub = data->nvmxsub;
  zdata = N_VDATA(z);

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      v = &(zdata[lx*NVARS + ly*nvmxsub]);
      gesl(P[lx][ly], NVARS, pivot[lx][ly], v);
    }
  }

  return(0);
}
