/*
 * -----------------------------------------------------------------
 * $Revision: 1.7.2.1 $
 * $Date: 2005-03-18 23:26:21 $
 * ----------------------------------------------------------------- 
 * Programmer(s): S. D. Cohen, A. C. Hindmarsh, M. R. Wittman, and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Parallel CVODES sensitivity test problem.
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions: 
 * 
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + Kh*(d/dy)^2 c(i)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where 
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 * Kh, Kh, V, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally.   The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * mesh, except for the advection term, which is treated with a biased
 * 3-point difference formula.
 * The initial profiles are proportional to a simple polynomial in x
 * and a tanh function in z.
 * 
 * The problem is solved by CVODES on NPE processors, treated as a 
 * rectangular process grid of size NPEX by NPEY, with NPE = NPEX*NPEY.
 * Each processor contains a subgrid of size MXSUB by MYSUB of the 
 * (x,y) mesh.  Thus the actual mesh sizes are MX = MXSUB*NPEX and 
 * MY = MYSUB*NPEY, and the ODE system size is neq = 2*MX*MY.
 * 
 * The solution with CVODES is done with the BDF/GMRES method (i.e.
 * using the CVSPGMR linear solver) and the block-diagonal part of the 
 * Newton matrix as a left preconditioner. A copy of the block-diagonal
 * part of the Jacobian is saved and conditionally reused within the
 * Precond routine.
 * 
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters Kh and Kv.
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on TRUE or FALSE
 * respectively).
 * 
 * Performance data and sampled solution values are printed at selected
 * output times, and all performance counters are printed on completion.
 * 
 * This version uses MPI for user routines. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"  /* definitions of realtype and booleantype           */
#include "cvodes.h"         /* main CVODES header file                           */
#include "iterative.h"      /* contains the enum for types of preconditioning    */
#include "cvspgmr.h"        /* use CVSPGMR linear solver each internal step      */
#include "smalldense.h"     /* use generic DENSE solver in preconditioning       */
#include "nvector_parallel.h" /* definitions of type N_Vector, macro NV_DATA_P   */
#include "sundialsmath.h"   /* contains SQR macro                                */
#include "mpi.h"

#include <memory.h>

/* Problem Constants */

#define NVARS        2             /* number of species         */
#define C1_SCALE     RCONST(1.0e6)         /* coefficients in initial profiles    */
#define C2_SCALE     RCONST(1.0e12)

#define T0           RCONST(0.0)           /* initial time */
#define NOUT         12            /* number of output times */
#define TWOHR        RCONST(7200.0)        /* number of seconds in two hours  */
#define HALFDAY      RCONST(4.32e4)        /* number of seconds in a half day */
#define PI       RCONST(3.1415926535898)   /* pi */

#define XMIN         RCONST(0.0)          /* grid boundaries in x  */
#define XMAX         RCONST(20.0)
#define YMIN         RCONST(30.0)          /* grid boundaries in y  */
#define YMAX         RCONST(50.0)

/* no. PEs in x and y directions of the PE matrix */
/* Total no. PEs = NPEX*NPEY */
#define NPEX         2
#define NPEY         2

#define MX           160   /* MX = number of x mesh points */
#define MY            40   /* MY = number of y mesh points */
                           /* Spatial mesh is MX by MY */

/* no. x and y points per subgrid */
#define MXSUB        MX/NPEX
#define MYSUB        MY/NPEY

/* These are (the only) testing-related definitions */
#define P5_KH        RCONST(4.0e-7)
#define P6_VEL       RCONST(1.0e-4)
#define P7_KV        RCONST(4.0e-7)
#define DELTA        RCONST(1.0)

#define OUTPUT_PATH  "data"         /* Directory where grid output goes */
#define PRINT_GRIDINFO -1
#define PRINT_BINARY -1
#define MAXL	     10

/* Sensitivity constants */
#define NP           8              /* number of problem parameters         */
#define NS           2              /* number of sensitivities              */
#define ZERO         RCONST(0.0)

/* CVodeMalloc Constants */
#define RTOL    RCONST(1.0e-5)            /* scalar relative tolerance */
#define FLOOR   RCONST(100.0)             /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */

/* User-defined matrix accessor macro: IJth */
#define IJth(a,i,j) (a[j-1][i-1])

/* Type : UserData
   contains problem constants, preconditioner blocks, pivot arrays,
   grid constants, and processor indices */

typedef struct {
  realtype *p;
  realtype q4, om, dx, dy;
  realtype **P[MXSUB][MYSUB], **Jbd[MXSUB][MYSUB];
  long int *pivot[MXSUB][MYSUB];
  long int my_pe, isubx, isuby, nvmxsub, nvmxsub2;
  MPI_Comm comm;
} *UserData;

/* Private Helper Functions */

static void ProcessArgs(int argc, char *argv[], int my_pe,
                        booleantype *sensi, int *sensi_meth, booleantype *err_con);
static void WrongArgs(int my_pe, char *name);

static UserData AllocUserData(void);
static void InitUserData(long int my_pe, MPI_Comm comm, UserData data);
static void FreeUserData(UserData data);

static void SetInitialProfiles(N_Vector u, UserData data);

static void PrintOutput(void *cvode_mem, int my_pe, MPI_Comm comm,
                        N_Vector u, realtype t);
static void PrintOutputS(int my_pe, MPI_Comm comm, N_Vector *uS);
#if PRINT_GRIDINFO == 1
static void PrintGrids(long int *xpes, long int *ypes, long int my_pe, 
                       MPI_Comm comm, N_Vector u);
static void PrintRow(realtype *row, FILE* f1, FILE* f2);
#endif
static void PrintFinalStats(void *cvode_mem, booleantype sensi);

static void BSend(MPI_Comm comm, long int my_pe, long int isubx, long int isuby,
                  long int dsizex, long int dsizey, realtype udata[]);
static void BRecvPost(MPI_Comm comm, MPI_Request request[], long int my_pe,
                      long int isubx, long int isuby,
                      long int dsizex, long int dsizey,
                      realtype uext[], realtype buffer[]);
static void BRecvWait(MPI_Request request[], long int isubx, long int isuby,
                      long int dsizex, realtype uext[], realtype buffer[]);

/* Functions Called by the CVODE Solver */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data);

static int Precond(realtype tn, N_Vector u, N_Vector fu, 
                   booleantype jok, booleantype *jcurPtr, 
                   realtype gamma, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
static int PSolve(realtype tn, N_Vector u, N_Vector fu, 
                  N_Vector r, N_Vector z, 
                  realtype gamma, realtype delta,
                  int lr, void *P_data, N_Vector vtemp);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt, int id);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  UserData data;

  int  my_pe, npes;
  MPI_Comm comm;

  void *cvode_mem;

  realtype abstol, reltol, t, tout;
  N_Vector u;
  long int neq, local_N;
  
  realtype *pbar;
  int is, *plist;
  N_Vector *uS;
  booleantype sensi, err_con;
  int sensi_meth;

  int iout, flag;

  /* following vars for timing purposes only */
  FILE *times;
  double mpi_time, mpi_time_start;

#if PRINT_GRIDINFO == 1
  /* for output of timestep and grid data */
  FILE *timesteps;
#endif

  /* Initialize pointers to NULL */
  data = NULL;
  cvode_mem = NULL;
  u = NULL;
  uS = NULL;
  pbar = NULL;
  plist = NULL;

  /* Set problem size neq */
  neq = NVARS*MX*MY;

  /* Get processor number and total number of pe's */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  mpi_time_start = MPI_Wtime();

  if (npes != NPEX*NPEY) {
    if (my_pe == 0)
      fprintf(stderr, "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n", npes,NPEX*NPEY);
    MPI_Finalize();
    return(1);
  }

  /* Process arguments */
  ProcessArgs(argc, argv, my_pe, &sensi, &sensi_meth, &err_con);

  /* Set local length */
  local_N = NVARS*MXSUB*MYSUB;

  /* Allocate user data block and load it */
  data = AllocUserData();
  if(check_flag((void *)data, "AllocUserData", 2, my_pe)) MPI_Abort(comm, 1);
  InitUserData(my_pe, comm, data);

  /* Allocate u, and set initial values and tolerances */

  u = N_VNew_Parallel(comm, local_N, neq);
  if(check_flag((void *)u, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  SetInitialProfiles(u, data);
  abstol = ATOL; reltol = RTOL;

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0, my_pe)) MPI_Abort(comm, 1);

  flag = CVodeSetFdata(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetFdata", 1, my_pe)) MPI_Abort(comm, 1);

  flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1, my_pe)) MPI_Abort(comm, 1);

  flag = CVodeMalloc(cvode_mem, f, T0, u, CV_SS, &reltol, &abstol);
  if(check_flag(&flag, "CVodeMalloc", 1, my_pe)) MPI_Abort(comm, 1);

  flag = CVSpgmr(cvode_mem, PREC_LEFT, MAXL);
  if(check_flag(&flag, "CVSpgmr", 1, my_pe)) MPI_Abort(comm, 1);

  flag = CVSpgmrSetPrecSetupFn(cvode_mem, Precond);
  if(check_flag(&flag, "CVSpgmrSetPrecSetupFn", 1, my_pe)) MPI_Abort(comm, 1);

  flag = CVSpgmrSetPrecSolveFn(cvode_mem, PSolve);
  if(check_flag(&flag, "CVSpgmrSetPrecSolveFn", 1, my_pe)) MPI_Abort(comm, 1);

  flag = CVSpgmrSetPrecData(cvode_mem, data);
  if(check_flag(&flag, "CVSpgmrSetPrecData", 1, my_pe)) MPI_Abort(comm, 1);

  if (my_pe == 0) {
    printf("MXxMY = %dx%d (%dx%d, %dx%d), PRINT_GRIDINFO = %d, MAXL = %d\n",
	   NPEX*MXSUB,NPEY*MYSUB,NPEX,NPEY,MXSUB,MYSUB,PRINT_GRIDINFO,MAXL);
    printf("\n2-species diurnal advection-diffusion problem\n\n");
  }

  /* SENSITIVTY */
  if( sensi) {

    pbar = (realtype *) malloc(NP*sizeof(realtype));
    if (check_flag((void *)pbar, "malloc", 2, my_pe)) MPI_Abort(comm, 1);

    for (is=0; is<NP; is++) pbar[is] = data->p[is];

    plist = (int *) malloc(NS * sizeof(int));
    if (check_flag((void *)plist, "malloc", 2, my_pe)) MPI_Abort(comm, 1);

    plist[0] = 6;
    plist[1] = 8;

    uS = N_VNewVectorArray_Parallel(NS, comm, local_N, neq);
    if (check_flag((void *)uS, "N_VNew", 0, my_pe)) MPI_Abort(comm, 1);
    for (is = 0; is < NS; is++)
      N_VConst(ZERO,uS[is]);

    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    if (check_flag(&flag, "CVodeSetSensErrCon", 1, my_pe)) MPI_Abort(comm, 1);

    flag = CVodeSetSensRho(cvode_mem, ZERO);
    if (check_flag(&flag, "CVodeSetSensRho", 1, my_pe)) MPI_Abort(comm, 1);

    flag = CVodeSetSensPbar(cvode_mem, pbar);
    if (check_flag(&flag, "CVodeSetSensPbar", 1, my_pe)) MPI_Abort(comm, 1);

    flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, data->p, plist, uS);
    if (check_flag(&flag, "CVodeSensMalloc", 1, my_pe)) MPI_Abort(comm, 1);

    if(my_pe == 0) {
      printf("Sensitivity: YES ");
      if(sensi_meth == CV_SIMULTANEOUS)   
        printf("( SIMULTANEOUS +");
      else 
        if(sensi_meth == CV_STAGGERED) printf("( STAGGERED +");
        else                           printf("( STAGGERED1 +");   
      if(err_con) printf(" FULL ERROR CONTROL )");
      else        printf(" PARTIAL ERROR CONTROL )");
    }

  } else {

    if(my_pe == 0) printf("Sensitivity: NO ");

  }

  if (my_pe == 0) {
    printf("\n\n");
    printf("========================================================================\n");
    printf("     T     Q       H      NST                    Bottom left  Top right \n");
    printf("========================================================================\n");
  }

#if PRINT_GRIDINFO == 1
  if (my_pe == 0) {
    FILE *griddata;
    char filename[30];

    sprintf(filename,"%s/timesteps",OUTPUT_PATH);
    timesteps = fopen(filename,"w");

    sprintf(filename,"%s/gridinfo",OUTPUT_PATH);
    griddata = fopen(filename,"w");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(griddata,"%d %d %d %d %Lg %Lg %Lg %Lg\n",NPEX,NPEY,MXSUB,MYSUB,
	    XMIN,XMAX,YMIN,YMAX);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(griddata,"%d %d %d %d %lg %lg %lg %lg\n",NPEX,NPEY,MXSUB,MYSUB,
	    XMIN,XMAX,YMIN,YMAX);
#else
    fprintf(griddata,"%d %d %d %d %g %g %g %g\n",NPEX,NPEY,MXSUB,MYSUB,
	    XMIN,XMAX,YMIN,YMAX);
#endif
    fclose(griddata);
  }
#endif

  /* In loop over output points, call CVode, print results, test for error */
  
  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_flag(&flag, "CVode", 1, my_pe)) break;
    
    PrintOutput(cvode_mem, my_pe, comm, u, t);
    
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, t, uS);
      if (check_flag(&flag, "CVodeGetSens", 1, my_pe)) break;
      PrintOutputS(my_pe, comm, uS);
    }
    if (my_pe == 0)
      printf("------------------------------------------------------------------------\n");

#if PRINT_GRIDINFO == 1
    {
      long int xpes[] = {0, 3}, ypes[] = {4, 4};
      
      /* Print both grids to files -
         Use xpes and ypes to specify PE ranges in x and y to print, or
         use NULL in their place to print out from all processors.
      */
      PrintGrids(NULL,NULL, my_pe, comm, u);
      
      if (my_pe == 0)
#if defined(SUNDIALS_EXTENDED_PRECISION)
        fprintf(timesteps,"%8.3Le\n",tout);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        fprintf(timesteps,"%8.3le\n",tout);
#else
        fprintf(timesteps,"%8.3e\n",tout);
#endif
    }
#endif

  }
  
#if PRINT_GRIDINFO == 1
  if (my_pe == 0)
    fclose(timesteps);
#endif

  /* Print final statistics */  

  if (my_pe == 0) PrintFinalStats(cvode_mem, sensi);

  /* do timing */
  mpi_time = MPI_Wtime() - mpi_time_start;
  if (my_pe == 0) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("seconds used: %Lg\n",mpi_time);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("seconds used: %lg\n",mpi_time);
#else
    printf("seconds used: %g\n",mpi_time);
#endif

    if (PRINT_GRIDINFO == 0) {
      times = fopen("test_times","a");
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(times,"%d %d %d %d %Lg\n",NPEX,NPEY,MXSUB,MYSUB,
	      mpi_time);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(times,"%d %d %d %d %lg\n",NPEX,NPEY,MXSUB,MYSUB,
	      mpi_time);
#else
      fprintf(times,"%d %d %d %d %g\n",NPEX,NPEY,MXSUB,MYSUB,
	      mpi_time);
#endif
      fclose(times);
    }
  }

  /* Free memory */
  N_VDestroy_Parallel(u);
  FreeUserData(data);
  CVodeFree(cvode_mem);

  MPI_Finalize();

  return(0);
}

/*********************** Private Helper Functions ************************/

/* Exit if arguments are incorrect */

static void ProcessArgs(int argc, char *argv[], int my_pe,
                        booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = FALSE;
  *sensi_meth = -1;
  *err_con = FALSE;

  if (argc < 2) WrongArgs(my_pe, argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = FALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = TRUE;
  else
    WrongArgs(my_pe, argv[0]);
  
  if (*sensi) {

    if (argc != 4)
      WrongArgs(my_pe, argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else 
      WrongArgs(my_pe, argv[0]);

    if (strcmp(argv[3],"t") == 0)
      *err_con = TRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = FALSE;
    else
      WrongArgs(my_pe, argv[0]);
  }

}

static void WrongArgs(int my_pe, char *name)
{
  if (my_pe == 0) {
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");
  }  
  MPI_Finalize();
  exit(0);
}


/* Allocate memory for data structure of type UserData */

static UserData AllocUserData(void)
{
  int lx, ly;
  UserData data;

  data = (UserData) malloc(sizeof *data);

  data->p = (realtype *) malloc(NP*sizeof(realtype));

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      (data->P)[lx][ly] = denalloc(NVARS);
      (data->Jbd)[lx][ly] = denalloc(NVARS);
      (data->pivot)[lx][ly] = denallocpiv(NVARS);
    }
  }

  return(data);
}

/* Load constants in data */

static void InitUserData(long int my_pe, MPI_Comm comm, UserData data)
{
  long int isubx, isuby;

  /* Set problem parameters */
  data->p[0]  = RCONST(1.63e-16); /* Q1  coeffs. q1, q2, c3             */
  data->p[1]  = RCONST(4.66e-16); /* Q2                                 */
  data->p[2]  = RCONST(3.7e16);   /* C3                                 */
  data->p[3]  = RCONST(22.62);    /* A3  coeff. in expression for q3(t) */
  data->p[4]  = RCONST(7.601);    /* A4  coeff. in expression for q4(t) */
  data->p[5]  = P5_KH;    /* KH  horizontal diffusivity Kh      */ 
  data->p[6]  = P6_VEL;   /* VEL advection velocity V           */
  data->p[7]  = P7_KV;    /* KV  vertical diffusivity Kv        */ 

  /* Set problem constants */
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/((realtype)(MX-1));
  data->dy = (YMAX-YMIN)/((realtype)(MY-1));

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

/* Free data memory */

static void FreeUserData(UserData data)
{
  int lx, ly;

  free(data->p);

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      denfree((data->P)[lx][ly]);
      denfree((data->Jbd)[lx][ly]);
      denfreepiv((data->pivot)[lx][ly]);
    }
  }

  free(data);
}

/* Set initial conditions in u */

static void SetInitialProfiles(N_Vector u, UserData data)
{
  long int isubx, isuby, lx, ly, jx, jy, offset;
  realtype dx, dy, x, y, cx, cy, xmid, ymid;
  realtype *udata;

  /* Set pointer to data array in vector u */

  udata = NV_DATA_P(u);

  /* Get mesh spacings, and subgrid indices for this PE */

  dx = data->dx;         dy = data->dy;
  isubx = data->isubx;   isuby = data->isuby;

  /* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. */

  offset = 0;
  xmid = RCONST(0.5)*(XMIN + XMAX);
  ymid = RCONST(0.5)*(YMIN + YMAX);
  for (ly = 0; ly < MYSUB; ly++) {
    jy = ly + isuby*MYSUB;
    y = YMIN + jy*dy;
    cy = RCONST(0.1)*(y - ymid);
    cy = RCONST(0.75) + RCONST(0.25)*tanh(cy/DELTA);
    for (lx = 0; lx < MXSUB; lx++) {
      jx = lx + isubx*MXSUB;
      x = XMIN + jx*dx;
      cx = SQR(RCONST(0.1)*(x - xmid));
      cx = RCONST(1.0) - cx + RCONST(0.5)*SQR(cx);
      udata[offset  ] = C1_SCALE*cx*cy;
      udata[offset+1] = C2_SCALE*cx*cy;
      offset = offset + 2;
    }
  }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, int my_pe, MPI_Comm comm,
                        N_Vector u, realtype t)
{
  long int nst;
  int qu;
  realtype hu, *udata, tempu[2];
  long int npelast, i0, i1;
  MPI_Status status;
  int flag;

  npelast = NPEX*NPEY - 1;
  udata = NV_DATA_P(u);

  /* Send c1,c2 at top right mesh point to PE 0 */
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

  /* On PE 0, receive c1,c2 at top right, then print performance data
     and sampled solution values */
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&tempu[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    flag = CVodeGetLastOrder(cvode_mem, &qu);
    flag = CVodeGetLastStep(cvode_mem, &hu);
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%8.3Le %2d  %8.3Le %5ld\n", t,qu,hu,nst);
    printf("                                Solution       ");
    printf("%12.4Le %12.4Le \n", udata[0], tempu[0]); 
    printf("                                               ");
    printf("%12.4Le %12.4Le \n", udata[1], tempu[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%8.3le %2d  %8.3le %5ld\n", t,qu,hu,nst);
    printf("                                Solution       ");
    printf("%12.4le %12.4le \n", udata[0], tempu[0]); 
    printf("                                               ");
    printf("%12.4le %12.4le \n", udata[1], tempu[1]);
#else
    printf("%8.3e %2d  %8.3e %5ld\n", t,qu,hu,nst);
    printf("                                Solution       ");
    printf("%12.4e %12.4e \n", udata[0], tempu[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", udata[1], tempu[1]);
#endif
  }
}

/* ======================================================================= */
/* Print sampled sensitivity values */

static void PrintOutputS(int my_pe, MPI_Comm comm, N_Vector *uS)
{
  realtype *sdata, temps[2];
  long int npelast, i0, i1;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;

  sdata = NV_DATA_P(uS[0]);

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
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le %12.4Le \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4Le %12.4Le \n", sdata[1], temps[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4le %12.4le \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4le %12.4le \n", sdata[1], temps[1]);
#else
    printf("%12.4e %12.4e \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", sdata[1], temps[1]);
#endif
  }

  sdata = NV_DATA_P(uS[1]);

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
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le %12.4Le \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4Le %12.4Le \n", sdata[1], temps[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4le %12.4le \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4le %12.4le \n", sdata[1], temps[1]);
#else
    printf("%12.4e %12.4e \n", sdata[0], temps[0]); 
    printf("                                               ");
    printf("%12.4e %12.4e \n", sdata[1], temps[1]);
#endif
  }
}



#if PRINT_GRIDINFO == 1

/* If xpes == NULL, print all pe's in the x direction.  Similarly for y.
   If xpes != NULL, print pe's with xpes[0] <= xpe number <= xpes[1].
   Similarly for y. */
static void PrintGrids(long int *xpes, long int *ypes, long int my_pe,
		       MPI_Comm comm, N_Vector u)
{
  realtype *udata;
  MPI_Status status;
  long int isubx, isuby;
  FILE *f1, *f2;
  char file1[30], file2[30];
  int rownum;
  static int called_count = 0;	/* used to number the time steps */

  isuby = my_pe/NPEX;
  isubx = my_pe - isuby*NPEX;

  /* check if this pe is among those to be printed out */
  if ((xpes == NULL || isubx >= xpes[0] && isubx <= xpes[1]) &&
      (ypes == NULL || isuby >= ypes[0] && isuby <= ypes[1])) {
    
    /* make filenames */
    sprintf(file1,"%s/c1.%ld.%ld.%d",OUTPUT_PATH,isubx,isuby,called_count);
    sprintf(file2,"%s/c2.%ld.%ld.%d",OUTPUT_PATH,isubx,isuby,called_count);
    
    /* open files */
    f1 = fopen(file1,"w");
    if (f1 == NULL) {
      fprintf(stderr,"\nFILE_ERROR(%ld): unable to open file %s\n\n", my_pe, file1);
      exit(1);
    }
    f2 = fopen(file2,"w");
    if (f2 == NULL) {
      fprintf(stderr,"\nFILE_ERROR(%ld): unable to open file %s\n\n", my_pe, file2);
      exit(1);
    }
    
    called_count++;
    
    udata = NV_DATA_P(u);

    for (rownum = 0; rownum < MYSUB; rownum++) {
      PrintRow(udata+rownum*2*MXSUB,f1,f2);
      
#if PRINT_BINARY != 1
      fprintf(f1,"\n");
      fprintf(f2,"\n");
#endif
    }

    /* close files */
    if (fclose(f1) == EOF) {
      fprintf(stderr,"\nFILE_ERROR(%ld): error closing file %s\n\n", my_pe, file1);
      exit(1);
    }
    if (fclose(f2) == EOF) {
      fprintf(stderr,"\nFILE_ERROR(%ld): error closing file %s\n\n", my_pe, file2);
      exit(1);
    }
  }
}

/* print one row of a subgrid (i.e. y constant, x changing) */
static void PrintRow(realtype *row, FILE* f1, FILE* f2)
{
  int i;
  
  for (i = 0; i < MXSUB; i++) {
#if PRINT_BINARY == 1
    fwrite(row+2*i,sizeof(realtype),1,f1);
    fwrite(row+2*i+1,sizeof(realtype),1,f2);
#else
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(f1,"%8.3Le ",row[2*i]);
    fprintf(f2,"%8.3Le ",row[2*i+1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(f1,"%8.3le ",row[2*i]);
    fprintf(f2,"%8.3le ",row[2*i+1]);
#else
    fprintf(f1,"%8.3e ",row[2*i]);
    fprintf(f2,"%8.3e ",row[2*i+1]);
#endif
#endif
  }
}

#endif


/* Print final statistics */

static void PrintFinalStats(void *cvode_mem, booleantype sensi) 
{
  long int lenrw, leniw ;
  long int lenrwSPGMR, leniwSPGMR;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeSPGMR;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int flag;

  flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_flag(&flag, "CVodeGetWorkSpace", 1, 0);
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1, 0);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1, 0);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1, 0);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1, 0);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1, 0);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1, 0);

  flag = CVSpgmrGetWorkSpace(cvode_mem, &lenrwSPGMR, &leniwSPGMR);
  check_flag(&flag, "CVSpgmrGetWorkSpace", 1, 0);
  flag = CVSpgmrGetNumLinIters(cvode_mem, &nli);
  check_flag(&flag, "CVSpgmrGetNumLinIters", 1, 0);
  flag = CVSpgmrGetNumPrecEvals(cvode_mem, &npe);
  check_flag(&flag, "CVSpgmrGetNumPrecEvals", 1, 0);
  flag = CVSpgmrGetNumPrecSolves(cvode_mem, &nps);
  check_flag(&flag, "CVSpgmrGetNumPrecSolves", 1, 0);
  flag = CVSpgmrGetNumConvFails(cvode_mem, &ncfl);
  check_flag(&flag, "CVSpgmrGetNumConvFails", 1, 0);
  flag = CVSpgmrGetNumRhsEvals(cvode_mem, &nfeSPGMR);
  check_flag(&flag, "CVSpgmrGetNumRhsEvals", 1, 0);

  if (sensi) {
    flag = CVodeGetNumSensRhsEvals(cvode_mem, &nfSe);
    check_flag(&flag, "CVodeGetNumSensRhsEvals", 1, 0);
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1, 0);
    flag = CVodeGetNumSensLinSolvSetups(cvode_mem, &nsetupsS);
    check_flag(&flag, "CVodeGetNumSensLinSolvSetups", 1, 0);
    flag = CVodeGetNumSensErrTestFails(cvode_mem, &netfS);
    check_flag(&flag, "CVodeGetNumSensErrTestFails", 1, 0);
    flag = CVodeGetNumSensNonlinSolvIters(cvode_mem, &nniS);
    check_flag(&flag, "CVodeGetNumSensNonlinSolvIters", 1, 0);
    flag = CVodeGetNumSensNonlinSolvConvFails(cvode_mem, &ncfnS);
    check_flag(&flag, "CVodeGetNumSensNonlinSolvConvFails", 1, 0);
  }  

  printf("\nFinal Statistics\n\n");
  printf("(Px*Mx)*(Py*My)*NVARS = (%d*%d)*(%d*%d)*%d = %d\n\n",
         NPEX,MXSUB,NPEY,MYSUB,NVARS,NPEX*MXSUB*NPEY*MYSUB*NVARS);
  printf("lrw      = %5ld    liw      = %5ld\n",   lenrw,      leniw);
  printf("lrwSPGMR = %5ld    liwSPGMR = %5ld\n\n", lenrwSPGMR, leniwSPGMR);
  printf("nst      = %5ld\n"                  , nst);
  printf("nfe      = %5ld    nfel     = %5ld\n" , nfe, nfeSPGMR);
  printf("nni      = %5ld    nli      = %5ld\n" , nni, nli);
  printf("nsetups  = %5ld    netf     = %5ld\n" , nsetups, netf);
  printf("npe      = %5ld    nps      = %5ld\n" , npe, nps);
  printf("ncfn     = %5ld    ncfl     = %5ld\n" , ncfn, ncfl); 

  if(sensi) {
    printf("\n");
    printf("nfSe     = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs    = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS     = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }  

  printf("\n");

}

/* Routine to send boundary data to neighboring PEs */

static void BSend(MPI_Comm comm, long int my_pe, long int isubx, long int isuby,
                  long int dsizex, long int dsizey, realtype udata[])
{
  int i, ly;
  long int offsetu, offsetbuf;
  realtype bufleft[NVARS*MYSUB], bufright[NVARS*MYSUB];

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

/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   buffer should be able to hold 2*NVARS*MYSUB real entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   request should have 4 entries, and should be passed both calls also. */

static void BRecvPost(MPI_Comm comm, MPI_Request request[], long int my_pe,
                      long int isubx, long int isuby,
                      long int dsizex, long int dsizey,
                      realtype uext[], realtype buffer[])
{
  long int offsetue;
  /* have bufleft and bufright use the same buffer */
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

/* Routine to finish receiving boundary data from neighboring PEs
   Notes:
   buffer should be able to hold 2*NVARS*MYSUB real entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   request should have 4 entries, and should be passed both calls also. */

static void BRecvWait(MPI_Request request[], long int isubx, long int isuby,
                      long int dsizex, realtype uext[], realtype buffer[])
{
  int i, ly;
  long int dsizex2, offsetue, offsetbuf;
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
    
    /* copy the buffer to uext */
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
    
    /* copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetue = (ly+2)*dsizex2 - NVARS;
      for (i = 0; i < NVARS; i++)
        uext[offsetue+i] = bufright[offsetbuf+i];
    }
  }
  
}

/***************** Functions Called by the CVODE Solver ******************/

/* f routine. Compute f(t,y). */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data)
{
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2;
  realtype q4coef, delx, dely, verdco, hordco, horaco, bcfac;
  realtype *udata, *dudata;
  int i, lx, ly, jx, jy;
  UserData data;
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV;
  MPI_Comm comm;
  long int my_pe, isubx, isuby, nvmxsub, nvmxsub2, nvmysub,
           offsetu, offsetue;
  realtype uext[NVARS*(MXSUB+2)*(MYSUB+2)];
  MPI_Request request[4];
  realtype buffer[2*NVARS*MYSUB];
  
  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);
  
  data = (UserData) f_data;

  /* Load problem coefficients and parameters */
  Q1  = data->p[0];
  Q2  = data->p[1];
  C3  = data->p[2];
  A3  = data->p[3];
  A4  = data->p[4];
  KH  = data->p[5];
  VEL = data->p[6];
  KV  = data->p[7];
  
  /* Get comm, my_pe, subgrid indices, data sizes */
  
  comm = data->comm;  my_pe = data->my_pe;
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;
  nvmxsub2 = data->nvmxsub2;
  nvmysub = NVARS*MYSUB;

  /* Set diurnal rate coefficients as functions of t */

  s = sin((data->om)*t);
  if (s > ZERO) {
    q3 = exp(-A3/s);
    data->q4 = exp(-A4/s);
  } else {
    q3 = ZERO;
    data->q4 = ZERO;
  }
  
  /* Make local copies of problem variables, for efficiency */
  
  q4coef = data->q4;
  delx = data->dx;
  dely = data->dy;
  
  verdco =  KV/SQR(dely);
  hordco  = KH/SQR(delx);
  horaco  = VEL/(RCONST(2.0)*delx);

  /* Copy local segment of u vector into the working array uext */

  offsetu = 0;
  offsetue = nvmxsub2 + NVARS;
  for (ly = 0; ly < MYSUB; ly++) {
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = udata[offsetu+i];
    offsetu = offsetu + nvmxsub;
    offsetue = offsetue + nvmxsub2;
  }

  /* Start receiving boundary data from neighboring PEs */

  BRecvPost(comm, request, my_pe, isubx, isuby, nvmxsub, nvmysub, uext, buffer);

  /* Send data from boundary of local grid to neighboring PEs */

  BSend(comm, my_pe, isubx, isuby, nvmxsub, nvmysub, udata);
  
  /* Finish receiving boundary data from neighboring PEs */
  
  BRecvWait(request, isubx, isuby, nvmxsub, uext, buffer);

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

  /* Loop over all grid points in local subgrid */

  for (ly = 0; ly < MYSUB; ly++) {

    jy = ly + isuby*MYSUB;
    
    /* Set vertical diffusion coefficients at jy +- 1/2 */

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
      rkin1 = -qq1 - qq2 + RCONST(2.0)*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms */

      c1dn = uext[offsetue-nvmxsub2];
      c2dn = uext[offsetue-nvmxsub2+1];
      c1up = uext[offsetue+nvmxsub2];
      c2up = uext[offsetue+nvmxsub2+1];
      vertd1 = verdco*(c1up - RCONST(2.0)*c1 + c1dn);
      vertd2 = verdco*(c2up - RCONST(2.0)*c2 + c2dn);

      /* Set horizontal diffusion and advection terms */

      c1lt = uext[offsetue-2];
      c2lt = uext[offsetue-1];
      c1rt = uext[offsetue+2];
      c2rt = uext[offsetue+3];
      hord1 = hordco*(c1rt - RCONST(2.0)*c1 + c1lt);
      hord2 = hordco*(c2rt - RCONST(2.0)*c2 + c2lt);
      bcfac = (jx==0 || jx==MX-1) ? ZERO : RCONST(1.0);
      horad1 = horaco*(RCONST(1.5)*c1rt - c1 - RCONST(0.5)*c1lt)*bcfac;
      horad2 = horaco*(RCONST(1.5)*c2rt - c2 - RCONST(0.5)*c2lt)*bcfac;

      /* Load all terms into udot */

      offsetu = lx*NVARS + ly*nvmxsub;
      dudata[offsetu]   = vertd1 + hord1 + horad1 + rkin1;
      dudata[offsetu+1] = vertd2 + hord2 + horad2 + rkin2;
    }
  }

}

/* Preconditioner setup routine. Generate and preprocess P. */
static int Precond(realtype tn, N_Vector u, N_Vector fu, 
                   booleantype jok, booleantype *jcurPtr, 
                   realtype gamma, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype c1, c2, diagd, diag;
  realtype q4coef, delx, dely, verdco, hordco, horaco;
  realtype **(*P)[MYSUB], **(*Jbd)[MYSUB];
  long int nvmxsub, *(*pivot)[MYSUB], ier, offset;
  int lx, ly, jx, jy, isubx, isuby;
  realtype *udata, **a, **j;
  UserData data;
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV;

  /* Make local copies of pointers in P_data, pointer to u's data,
     and PE index pair */

  data = (UserData) P_data;

  /* Load problem coefficients and parameters */
  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];
  A3 = data->p[3];
  A4 = data->p[4];
  KH = data->p[5];
  VEL = data->p[6];
  KV  = data->p[7];

  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  udata = NV_DATA_P(u);
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;

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
    delx = data->dx;
    dely = data->dy;
    
    verdco  = KV/SQR(dely);
    hordco  = KH/SQR(delx);
    horaco  = VEL/(RCONST(2.0)*delx);
    
    /* Compute 2x2 diagonal Jacobian blocks (using q4 values
       computed on the last f call).  Load into P. */
    
    diagd = -(RCONST(2.0)*verdco + RCONST(2.0)*hordco);
    for (ly = 0; ly < MYSUB; ly++) {
      jy = ly + isuby*MYSUB;
      for (lx = 0; lx < MXSUB; lx++) {
        jx = lx + isubx*MXSUB;
        diag = (jx==0 || jx==MX-1) ? diagd : diagd - horaco;
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


/* Preconditioner solve routine */
static int PSolve(realtype tn, N_Vector u, N_Vector fu, 
                  N_Vector r, N_Vector z, 
                  realtype gamma, realtype delta,
                  int lr, void *P_data, N_Vector vtemp)
{
  realtype **(*P)[MYSUB];
  long int nvmxsub, *(*pivot)[MYSUB];
  int lx, ly;
  realtype *zdata, *v;
  UserData data;
  
  /* Extract the P and pivot arrays from P_data */
  
  data = (UserData) P_data;
  nvmxsub = data->nvmxsub;
  P = data->P;
  pivot = data->pivot;
  zdata = NV_DATA_P(z);

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z.
     First copy vector r to z. */

  N_VScale(1.0, r, z);
  
  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      v = &(zdata[lx*NVARS + ly*nvmxsub]);
      gesl(P[lx][ly], NVARS, pivot[lx][ly], v);
    }
  }
  
  return(0);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, char *funcname, int opt, int id)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n", id, funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with flag = %d\n\n", id, funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n", id, funcname);
    return(1); }

  return(0);
}
