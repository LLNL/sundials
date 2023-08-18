/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds, Jean Sexton @ SMU
 *                Slaven Peles @ LLNL
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *-----------------------------------------------------------------
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(y) = Kv0*exp(y/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * mesh, with simple polynomial initial profiles.
 *
 * The problem is solved by ARKODE on nprocs processors, treated
 * as a rectangular process grid of size nprocsx by nprocsy, with
 * nprocs = nprocsx*nprocsy. Each processor contains a subgrid of
 * size local_Mx by local_My of the (x,y) mesh. The global mesh
 * size is Mx by My, and the ODE system size is N = 2*Mx*My.
 *
 * The solution is done with the DIRK/GMRES method (i.e. using the
 * SUNLinSol_SPGMR linear solver) and the block-diagonal part of
 * the Newton matrix as a left preconditioner. A copy of the
 * block-diagonal part of the Jacobian is saved and conditionally
 * reused within the preconditioner routine.
 *
 * Performance data and sampled solution values are printed at
 * selected output times, and all performance counters are printed
 * on completion.
 *
 * This example uses Hypre vector and MPI parallelization. User is
 * expected to be familiar with the Hypre library.
 *-----------------------------------------------------------------
 * Data flow: SERIAL backend
 *                  ╭                          ╮
 *                  ┊   ╭─<───[+dt]───<─╮      ┊
 *                  ┊ udata╶─>─[f]─>─╴udotdata ┊
 *                  ┊   ┊       ╰──<─╮         ┊
 *         [Send] ^ ┊   ╰┄┄┄┄>┄┄┄┄┄┄╮│         ┊ v [Irecv]
 *              <╶┼────────────<──╴bufs╶──<──────┼╴<
 *                v ╰                          ╯ ^
 *
 * Data flow: CUDA/HIP backend
 *               ╭                                ╮
 *               ┊ [fill]  ╭───<─[+dt]─<─udotdata ┊
 *               ┊  ╭┄┄╴udata╶─>──[f]──>───╯      ┊
 *               ┊  ╰>bufs_dev╶>───╯      DEVICE  ┊
 *        [Send] ┊- - - -┊╰────<───╮- - - - - - - ┊ [Irecv]
 *             ^ ┊       ╰┄┄┄┄┄>┄┄╮│       HOST   ┊ v 
 *           <╶┼─────────────<──╴bufs╶──<───────────┼╴<
 *             v ╰                                ╯ ^
 *                                              
 * Note: bufs[_dev] is a combined buffer for sends, receives, and
 *       host-device transfers (with CUDA/HIP backends). It is
 *       conceptually subdivided into 4 send and 4 receive buffers,
 *       [sendB,sendT,sendL,sendR,recvB,recvT,recvL,recvR], and its
 *       length is 2*NVARS*(2*local_Mx + 2*local_My).
 *-----------------------------------------------------------------
 * Execution: mpiexec -n N ark_diurnal_kry_ph nprocsx nprocsy Mx My
 * 
 * Four (4) arguments required:
 * - nprocsx (# of procs in x)
 * - nprocsy (# of procs in y)
 * - Mx      (# of interior grid points in x)
 * - My      (# of interior grid points in y)
 * 
 * Arguments must obey:
 *  (1) nprocsx*nprocsy = N
 *  (2) Mx >= nprocsx (at least one interior grid point per proc)
 *  (3) My >= nprocsy (at least one interior grid point per proc)
 * 
 * See other problem constants below.
 *-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>        /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_parhyp.h>       /* declaration of N_Vector              */
#include <sunlinsol/sunlinsol_spgmr.h>    /* access to SPGMR SUNLinearSolver      */
#include <sundials/sundials_dense.h>      /* prototypes for small dense fcts.     */
#include <sundials/sundials_types.h>      /* definitions of realtype, booleantype */
#include <sundials/sundials_math.h>       /* definition of macros SUNSQR and EXP  */
#include <mpi.h>                          /* MPI constants and types              */

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>

/* Problem Constants */
#define NVARS        2                    /* number of species         */
#define KH           RCONST(4.0e-6)       /* horizontal diffusivity Kh */
#define VEL          RCONST(0.001)        /* advection velocity V      */
#define KV0          RCONST(1.0e-8)       /* coefficient in Kv(y)      */
#define Q1           RCONST(1.63e-16)     /* coefficients q1, q2, c3   */
#define Q2           RCONST(4.66e-16)
#define C3           RCONST(3.7e16)
#define A3           RCONST(22.62)        /* coefficient in expression for q3(t) */
#define A4           RCONST(7.601)        /* coefficient in expression for q4(t) */
#define C1_SCALE     RCONST(1.0e6)        /* coefficients in initial profiles    */
#define C2_SCALE     RCONST(1.0e12)

#define T0           RCONST(0.0)          /* initial time */
#define NOUT         12                   /* number of output times */
#define HALFHR       RCONST(1800.0)
#define TWOHR        RCONST(7200.0)       /* number of seconds in two hours  */
#define HALFDAY      RCONST(4.32e4)       /* number of seconds in a half day */
#define PI        RCONST(3.1415926535898) /* pi */

#define XMIN         RCONST(0.0)          /* grid boundaries in x  */
#define XMAX         RCONST(20.0)
#define YMIN         RCONST(30.0)         /* grid boundaries in y  */
#define YMAX         RCONST(50.0)

/* initialization constants */
#define RTOL    RCONST(1.0e-5)            /* scalar relative tolerance */
#define FLOOR   RCONST(100.0)             /* value of C1 or C2 at which tolerances */
                                          /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)              /* scalar absolute tolerance */

/* --- Backend-specific definitions --- */

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA)
#define NV_ADD_LANG_PREFIX_PH(token) cuda##token // token pasting; expands to ```cuda[token]```
#define NV_VERIFY_CALL_PH SUNDIALS_CUDA_VERIFY

#elif defined(SUNDIALS_HYPRE_BACKENDS_HIP)
#define NV_ADD_LANG_PREFIX_PH(token) hip##token // token pasting; expands to ```hip[token]```
#define NV_VERIFY_CALL_PH SUNDIALS_HIP_VERIFY
#endif

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA) || defined(SUNDIALS_HYPRE_BACKENDS_HIP)
#define SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP
#endif

/* --- Error-checking macros --- */

#define CHECK_LAST_ERROR() NV_VERIFY_CALL_PH(NV_ADD_LANG_PREFIX_PH(GetLastError)());

/* User-defined matrix accessor macro: IJth

   IJth is defined in order to write code which indexes into dense
   matrices with a (row,column) pair, where 1 <= row,column <= NVARS.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NVARS. The small matrix routines in sundials_dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */
#define IJth(a,i,j) (a[j-1][i-1])

/* User-defined MPI assert macro */
#define MPI_ASSERT(expr,msg,comm,myproc,code)                              \
  if(!(expr)) {                                                            \
    fprintf(stderr, "ERROR in %s (%s line %d): %s",                        \
        __func__, __FILE__, __LINE__, "\n──> " msg "\n──> Aborting...\n"); \
    MPI_Abort(comm,code);                                                  \
  }



/* Type : UserData
   contains problem constants, preconditioner blocks, pivot arrays,
   grid constants, and processor indices, as well as data needed
   for the preconditiner */
typedef struct {
  realtype q4, om, dx, dy, hdco, haco, vdco;  /* quantities (q4 depends on time of day)    */
  int M, Mx, My, local_M, local_Mx, local_My; /* # of interior grid points M=Mx*My         */
  int N, local_N;                             /* # of equations N=2*M                      */
  int mybase, mybasex, mybasey;               /* global offset of first local grid point   */
  int nprocs, nprocsx, nprocsy;               /* # of processes nprocs=nprocsx*nprocsy     */
  int myproc, myprocx, myprocy;               /* my process indices, flat & Cartesian      */
  int dsizex, dsizey;//, dsizex2, dsizey2;    /* x/y size of u & uext                      */
  int isbottom, istop, isleft, isright;       /* (bool) does process abut given boundary   */
  MPI_Comm comm;

  realtype *bufs;                             /* combined send/receieve buffers            */
  int      bufslen;                           /* total length of bufs                      */
  int      sendB, sendT, sendL, sendR;        /* index in bufs of host send buffers        */
  int      recvB, recvT, recvL, recvR;        /* index in bufs of host receive buffers     */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  /* for preconditioner (SERIAL backend only, for now) */
  realtype ****P, ****Jbd;                    /* 2D arrays of dense matrices (** to **)    */
  sunindextype ***pivot;                      /* 2D arrays of index arrays (** to * )      */
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  /* device arrays (CUDA/HIP backends, only) */
  void     *data_dev;                         /* give device a copy of UserData            */
  realtype *bufs_dev;                         /* give device a buffer for sends/receives   */
#endif
} *UserData;

/* Private Helper Functions */
static void InitUserData(UserData data, MPI_Comm comm, int nprocsx, int nprocsy, int Mx, int My);
static void FreeUserData(UserData data);
static void SetInitialProfiles(UserData data, HYPRE_IJVector Uij);
static void PrintOutput(UserData data, void *arkode_mem, N_Vector u, realtype t);
static void PrintFinalStats(void *arkode_mem);
static void BSend(UserData data, realtype udata[]);
static void BRecvPost(UserData data, MPI_Request request[]);
static void BRecvWait(UserData data, MPI_Request request[]);
static void ucomm(UserData data, realtype t, N_Vector u);

/* RHS calculation and preconditioner functions for SERIAL hypre backend */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
static void fcalc(UserData data, realtype t, realtype udata[], realtype dudata[]);
static int Precond(realtype tn, N_Vector u, N_Vector fu,
                   booleantype jok, booleantype *jcurPtr,
                   realtype gamma, void *user_data);
static int PSolve(realtype tn, N_Vector u, N_Vector fu,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data);
#endif

/* RHS kernel and send buffer population kernel for CUDA/HIP hypre backend */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
__global__ void fcalcKernel(UserData data_dev, realtype t, realtype *udata, realtype *udotdata);
__global__ void FillSendBufferKernel(UserData data_dev, realtype *udata);      
#endif

/* Functions Called by the Solver */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt, int id);


/***************************** Main Program ******************************/
int main(int argc, char *argv[])
{
  int             iout, flag, myproc, nprocs, nprocsx, nprocsy, Mx, My;
  void            *arkode_mem;
  realtype        abstol, reltol, t, tout;

  UserData        data;
  N_Vector        u;
  SUNContext      sunctx;
  SUNLinearSolver LS;
  MPI_Comm        comm;
  HYPRE_ParVector Upar; /* Declare HYPRE parallel vector */
  HYPRE_IJVector  Uij;  /* Declare "IJ" interface to HYPRE vector */

  data       = NULL;
  u          = NULL;
  LS         = NULL;
  arkode_mem = NULL;

  /* Get processor number and total number of proc's */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myproc);

  /* Parse inputs (required: nprocsx, nprocsy, Mx, My) */
  MPI_ASSERT
  (
    argc>=5,
    "ERROR: FOUR (4) inputs required:\n\
     - nprocsx (# of procs in x)\n\
     - nprocsy (# of procs in y)\n\
     - Mx      (# of interior grid points in x)\n\
     - My      (# of interior grid points in y)\n",
    comm,myproc,1
  )
  nprocsx = atoi(argv[1]);
  nprocsy = atoi(argv[2]);
  Mx      = atoi(argv[3]);
  My      = atoi(argv[4]);
  MPI_ASSERT(nprocsx*nprocsy==nprocs,"ERROR: nprocsx*nprocsy =/= N (Requested proc grid size does not equal available procs)",comm,myproc,1)
  MPI_ASSERT(Mx>=nprocsx,"ERROR: Mx < nprocsx (We require at least one interior grid point per process)",comm,myproc,1)
  MPI_ASSERT(My>=nprocsy,"ERROR: My < nprocsy (We require at least one interior grid point per process)",comm,myproc,1)

  /* Create the SUNDIALS context object for this simulation */
  flag = SUNContext_Create(&comm, &sunctx);
  if (check_flag(&flag, "SUNContext_Create", 1, myproc)) return 1;

  /* Allocate and load user data block; allocate preconditioner block */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2, myproc)) MPI_Abort(comm, 1);
  InitUserData(data, comm, nprocsx, nprocsy, Mx, My);

  /* Allocate hypre vector */
  HYPRE_IJVectorCreate(comm, NVARS*data->mybase, NVARS*data->mybase+data->local_N-1, &Uij);
  HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(Uij);

  /* Set initial values and allocate u */
  SetInitialProfiles(data, Uij);
  HYPRE_IJVectorAssemble(Uij);
  HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

  u = N_VMake_ParHyp(Upar, sunctx);  /* Create wrapper u around hypre vector */
  if (check_flag((void *)u, "N_VNew", 0, myproc)) MPI_Abort(comm, 1);

  /* Set tolerances */
  abstol = ATOL; reltol = RTOL;

  /* Create SPGMR solver structure -- use left preconditioning
     and the defauL Krylov dimension maxl */
  // Note: For CUDA/HIP, we could use MAGMA (Matrix Algebra on GPU and MuLi-
  // core Architectures) to write preconditioner routines. Without it, for
  // simplicity, we forego preconditioning in the CUDA/HIP case.
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  LS = SUNLinSol_SPGMR(u, SUN_PREC_LEFT, 0, sunctx);
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, 0, sunctx);
#endif
  if (check_flag((void *)LS, "SUNLinSol_SPGMR", 0, myproc)) MPI_Abort(comm, 1);

  /* Call ARKStepCreate to initialize the integrator memory and specify the
     user's right hand side function in u'=fi(t,u) [here fe is NULL],
     the inital time T0, and the initial dependent variable vector u. */
  arkode_mem = ARKStepCreate(NULL, f, T0, u, sunctx);
  if (check_flag((void *)arkode_mem, "ARKStepCreate", 0, myproc)) MPI_Abort(comm, 1);

  /* Set the pointer to user-defined data */
  flag = ARKStepSetUserData(arkode_mem, data);
  if (check_flag(&flag, "ARKStepSetUserData", 1, myproc)) MPI_Abort(comm, 1);

  /* Call ARKStepSetMaxNumSteps to increase defauL */
  flag = ARKStepSetMaxNumSteps(arkode_mem, 10000);
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1, myproc)) return(1);

  /* Call ARKStepSStolerances to specify the scalar relative tolerance
     and scalar absolute tolerances */
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKStepSStolerances", 1, myproc)) return(1);

  /* Attach SPGMR solver structure to ARKStep interface */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1, myproc)) MPI_Abort(comm, 1);

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  flag = ARKStepSetPreconditioner(arkode_mem, Precond, PSolve);
  if (check_flag(&flag, "ARKStepSetPreconditioner", 1, myproc)) MPI_Abort(comm, 1);
#endif

  /* Print heading */
  printf("Process myproc=%d (myprocx=%d, myprocy=%d) shape: local_Mx=%d by local_My=%d (local_M=%d)\n",
         data->myproc,data->myprocx,data->myprocy,data->local_Mx,data->local_My,data->local_M);
  if (myproc == 0) {
    printf("\n2-species diurnal advection-diffusion problem\n\n");
  }

  /* In loop over output points, call ARKStepEvolve, print resuLs, test for error */
  PrintOutput(data, arkode_mem, u, 0.0);
  for (iout=1, tout=TWOHR; iout<=NOUT; iout++, tout+=TWOHR) {
    flag = ARKStepEvolve(arkode_mem, tout, u, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKStepEvolve", 1, myproc)) break;
    PrintOutput(data, arkode_mem, u, t);
  }

  /* Print final statistics */
  if (myproc == 0) PrintFinalStats(arkode_mem);

  /* Free memory */
  N_VDestroy(u);              /* Free hypre vector wrapper */
  HYPRE_IJVectorDestroy(Uij); /* Free the underlying hypre vector */
  FreeUserData(data);
  ARKStepFree(&arkode_mem);
  SUNLinSolFree(LS);
  SUNContext_Free(&sunctx);   /* Free context */
  MPI_Finalize();

  return(0);
}


/*********************** Private Helper Functions ************************/

/* Load constants in data */
static void InitUserData(UserData data, MPI_Comm comm, int nprocsx, int nprocsy, int Mx, int My)
{
  int nperprocx, nperprocy, nremx, nremy;

  /* Set problem constants */
  data->om       = PI/HALFDAY;
  data->dx       = (XMAX-XMIN)/((realtype)(Mx-1));
  data->dy       = (YMAX-YMIN)/((realtype)(My-1));
  data->hdco     = KH/SUNSQR(data->dx);
  data->haco     = VEL/(RCONST(2.0)*data->dx);
  data->vdco     = (RCONST(1.0)/SUNSQR(data->dy))*KV0;

  /* Set machine-related constants */
  data->comm     = comm;
  data->nprocsx  = nprocsx;
  data->nprocsy  = nprocsy;
  MPI_Comm_size(comm, &data->nprocs);
  MPI_Comm_rank(comm, &data->myproc);

  /* myprocx and myprocy are the proc grid indices corresponding to myproc */
  // Note: grid indexed left to right (in x), bottom to top (in y)
  data->myprocy  = data->myproc / nprocsx;
  data->myprocx  = data->myproc - data->myprocy*nprocsx;

  /* Set boundary booleans (whether proc lies on given boundary) */
  data->isbottom = (data->myprocy==0)         ? 1 : 0;
  data->istop    = (data->myprocy==nprocsy-1) ? 1 : 0;
  data->isleft   = (data->myprocx==0)         ? 1 : 0;
  data->isright  = (data->myprocx==nprocsx-1) ? 1 : 0;

  /* Set partitioning (# of interior grid points M=Mx*My, # eqs N=2*M) */
  data->Mx       = Mx;
  data->My       = My;
  data->M        = Mx*My;
  data->N        = NVARS*data->M;
  nperprocx      = Mx / nprocsx;
  nremx          = Mx - nprocsx*nperprocx;
  nperprocy      = My / nprocsy;
  nremy          = My - nprocsy*nperprocy;
  data->local_Mx = (data->myprocx < nremx) ? nperprocx+1 : nperprocx;
  data->local_My = (data->myprocy < nremy) ? nperprocy+1 : nperprocy;
  data->local_M  = data->local_Mx*data->local_My;
  data->local_N  = NVARS*data->local_M;

  /* Allocate send/receive buffers and get individual buffer indices */
  // Note: bufs[_dev] is combined buffer storage: [sendB,sendT,sendL,sendR,recvB,recvT,recvL,recvR]
  data->bufslen  = 2*NVARS*(2*data->local_Mx+2*data->local_My);
  data->bufs     = (realtype*) malloc(data->bufslen*sizeof(realtype));
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Malloc)(&data->bufs_dev, data->bufslen*sizeof(realtype));
  CHECK_LAST_ERROR();
#endif
  data->sendB    = 0;
  data->sendT    = data->sendB + NVARS*data->local_Mx;
  data->sendL    = data->sendT + NVARS*data->local_Mx;
  data->sendR    = data->sendL + NVARS*data->local_My;
  data->recvB    = data->sendR + NVARS*data->local_My;
  data->recvT    = data->recvB + NVARS*data->local_Mx;
  data->recvL    = data->recvT + NVARS*data->local_Mx;
  data->recvR    = data->recvL + NVARS*data->local_My;
  
  /* Calculate global offset of first local grid point */
  // Note: This is *not* an equation offset. There are NVARS equations per grid point.
  data->mybasex = (data->myprocx < nremx) ? (nperprocx+1)*data->myprocx : nremx+(nperprocx)*data->myprocx;
  data->mybasey = (data->myprocy < nremy) ? (nperprocy+1)*data->myprocy : nremy+(nperprocy)*data->myprocy;
  data->mybase  = (data->mybasey*data->Mx)+data->mybasex;

  /* Set the sizes of a boundary x/y-line in u and uext */
  data->dsizex  = NVARS*data->local_Mx;
  data->dsizey  = NVARS*data->local_My;

  /* Preconditioner-related fields */
  // Note: To each interior grid point, we associate two real matrices and one index array
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  int lx, ly;
  data->P        = (realtype****)    malloc(data->local_Mx*sizeof(realtype***));
  data->Jbd      = (realtype****)    malloc(data->local_Mx*sizeof(realtype***));
  data->pivot    = (sunindextype***) malloc(data->local_Mx*sizeof(sunindextype**));
  for (lx = 0; lx < data->local_Mx; lx++) {
    data->P[lx]     = (realtype***)    malloc(data->local_My*sizeof(realtype**));
    data->Jbd[lx]   = (realtype***)    malloc(data->local_My*sizeof(realtype**));
    data->pivot[lx] = (sunindextype**) malloc(data->local_My*sizeof(sunindextype*));
    for (ly = 0; ly < data->local_My; ly++) {
      (data->P)[lx][ly]     = SUNDlsMat_newDenseMat(NVARS, NVARS);
      (data->Jbd)[lx][ly]   = SUNDlsMat_newDenseMat(NVARS, NVARS);
      (data->pivot)[lx][ly] = SUNDlsMat_newIndexArray(NVARS);
    }
  }
#endif

  /* Give device a copy of UserData */
  // Note: In particular, data_dev will contain all constants and the correct pointer to bufs_dev
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Malloc)((void**)&data->data_dev,sizeof *data);
  CHECK_LAST_ERROR();
  NV_ADD_LANG_PREFIX_PH(Memcpy)(data->data_dev,data,sizeof *data,NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
  CHECK_LAST_ERROR();
#endif
}

/* Free user data memory */
static void FreeUserData(UserData data)
{
  /* Free bufs */
  free(data->bufs);

  /* Free preconditioner-related fields */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  int lx, ly;
  for (lx = 0; lx < data->local_Mx; lx++) {
    for (ly = 0; ly < data->local_My; ly++) {
      SUNDlsMat_destroyMat((data->P)[lx][ly]);
      SUNDlsMat_destroyMat((data->Jbd)[lx][ly]);
      SUNDlsMat_destroyArray((data->pivot)[lx][ly]);
    }
    free(data->P[lx]);
    free(data->Jbd[lx]);
    free(data->pivot[lx]);
  }
  free(data->P);
  free(data->Jbd);
  free(data->pivot);
#endif

  /* Free device UserData and bufs */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Free)(data->data_dev);
  CHECK_LAST_ERROR();
  NV_ADD_LANG_PREFIX_PH(Free)(data->bufs_dev);
  CHECK_LAST_ERROR();
#endif

  /* Free UserData */
  free(data);
}

/* Set initial conditions in u */
static void SetInitialProfiles(UserData data, HYPRE_IJVector Uij)
{
  int          lx, ly, jx, jy;
  sunindextype offset;
  realtype     dx, dy, x, y, cx, cy, xmid, ymid;
  realtype     *udata;
  HYPRE_Int    *iglobal;

  /* Set pointer to data array in vector u */
  udata   = (realtype*) malloc(data->local_N*sizeof(realtype));
  iglobal = (HYPRE_Int*) malloc(data->local_N*sizeof(HYPRE_Int));

  /* Get mesh spacings for this proc */
  dx      = data->dx;
  dy      = data->dy;

  /* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. */
  offset = 0;
  xmid   = RCONST(0.5)*(XMIN + XMAX);
  ymid   = RCONST(0.5)*(YMIN + YMAX);
  for (ly = 0; ly < data->local_My; ly++) {
    jy = ly + data->mybasey;
    y  = YMIN + jy*dy;
    cy = SUNSQR(RCONST(0.1)*(y - ymid));
    cy = RCONST(1.0) - cy + RCONST(0.5)*SUNSQR(cy);
    for (lx = 0; lx < data->local_Mx; lx++) {
      jx = lx + data->mybasex;
      x  = XMIN + jx*dx;
      cx = SUNSQR(RCONST(0.1)*(x - xmid));
      cx = RCONST(1.0) - cx + RCONST(0.5)*SUNSQR(cx);
      iglobal[offset] = data->mybase + offset;
      udata[offset++] = C1_SCALE*cx*cy;
      iglobal[offset] = data->mybase + offset;
      udata[offset++] = C2_SCALE*cx*cy;
    }
  }
  #if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  /* For SERIAL backend, HYPRE_IJVectorSetValues expects host pointers */
  HYPRE_IJVectorSetValues(Uij, data->local_N, iglobal, udata);

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  /* With GPU backend, HYPRE_IJVectorSetValues expects device pointers */
  realtype*  udata_dev;
  HYPRE_Int* iglobal_dev;

  /* Allocate device arrays */
  NV_ADD_LANG_PREFIX_PH(Malloc)(&udata_dev, data->local_N*sizeof(realtype));
  NV_ADD_LANG_PREFIX_PH(Malloc)(&iglobal_dev, data->local_N*sizeof(HYPRE_Int));
  
  /* Copy host data to device */
  NV_ADD_LANG_PREFIX_PH(Memcpy)(udata_dev, udata, data->local_N*sizeof(realtype), NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
  NV_ADD_LANG_PREFIX_PH(Memcpy)(iglobal_dev, iglobal, data->local_N*sizeof(HYPRE_Int), NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));

  /* Set hypre vector values from device arrays, then free the device arrays */
  HYPRE_IJVectorSetValues(Uij, data->local_N, iglobal_dev, udata_dev);
  NV_ADD_LANG_PREFIX_PH(Free)(udata_dev);
  NV_ADD_LANG_PREFIX_PH(Free)(iglobal_dev);
#endif

  free(udata);
  free(iglobal);
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */
static void PrintOutput(UserData data, void *arkode_mem, N_Vector u, realtype t)
{
  int flag;
  realtype hu, *udata, tempu[2];
  int npelast;
  sunindextype i0, i1;
  long int nst;
  MPI_Status status;
  HYPRE_ParVector uhyp;

  npelast = data->nprocs - 1;

  uhyp  = N_VGetVector_ParHyp(u);
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));

  /* Send c1,c2 at top right mesh point to proc 0 */
  if (data->myproc == npelast) {
    i0 = data->local_N - NVARS;
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
    i1 = i0 + 1;
    tempu[0] = udata[i0];
    tempu[1] = udata[i1];
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
    NV_ADD_LANG_PREFIX_PH(Memcpy)(tempu,&udata[i0],NVARS*sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyDeviceToHost));
    CHECK_LAST_ERROR();
#endif      
    if (npelast != 0)
      MPI_Send(tempu, NVARS, MPI_SUNREALTYPE, 0, 0, data->comm);
  }

  /* On proc 0, receive c1,c2 at top right, then print performance data
     and sampled solution values */
  if (data->myproc == 0) {
    if (npelast != 0)
      MPI_Recv(tempu, NVARS, MPI_SUNREALTYPE, npelast, 0, data->comm, &status);
    flag = ARKStepGetNumSteps(arkode_mem, &nst);
    check_flag(&flag, "ARKStepGetNumSteps", 1, data->myproc);
    flag = ARKStepGetLastStep(arkode_mem, &hu);
    check_flag(&flag, "ARKStepGetLastStep", 1, data->myproc);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("t = %.2Le   no. steps = %ld   stepsize = %.2Le\n",
           t, nst, hu);
    printf("At bottom left:  c1, c2 = %12.3Le %12.3Le \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3Le %12.3Le \n\n", tempu[0], tempu[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("t = %.2e   no. steps = %ld   stepsize = %.2e\n",
           t, nst, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#else
    printf("t = %.2e   no. steps = %ld   stepsize = %.2e\n",
           t, nst, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#endif
  }
}

/* Print final statistics contained in iopt */
static void PrintFinalStats(void *arkode_mem)
{
  long int lenrw, leniw;
  long int lenrwLS, leniwLS;
  long int nst, nfe, nfi, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  flag = ARKStepGetWorkSpace(arkode_mem, &lenrw, &leniw);
  check_flag(&flag, "ARKStepGetWorkSpace", 1, 0);
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKStepGetNumSteps", 1, 0);
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1, 0);
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1, 0);
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKStepGetNumErrTestFails", 1, 0);
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1, 0);
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1, 0);

  flag = ARKStepGetLinWorkSpace(arkode_mem, &lenrwLS, &leniwLS);
  check_flag(&flag, "ARKStepGetLinWorkSpace", 1, 0);
  flag = ARKStepGetNumLinIters(arkode_mem, &nli);
  check_flag(&flag, "ARKStepGetNumLinIters", 1, 0);
  flag = ARKStepGetNumPrecEvals(arkode_mem, &npe);
  check_flag(&flag, "ARKStepGetNumPrecEvals", 1, 0);
  flag = ARKStepGetNumPrecSolves(arkode_mem, &nps);
  check_flag(&flag, "ARKStepGetNumPrecSolves", 1, 0);
  flag = ARKStepGetNumLinConvFails(arkode_mem, &ncfl);
  check_flag(&flag, "ARKStepGetNumLinConvFails", 1, 0);
  flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1, 0);

  printf("\nFinal Statistics: \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n", lenrw, leniw);
  printf("lenrwls = %5ld     leniwls = %5ld\n", lenrwLS, leniwLS);
  printf("nst     = %5ld     nfe     = %5ld\n", nst, nfe);
  printf("nfi     = %5ld     nfels   = %5ld\n", nfi, nfeLS);
  printf("nni     = %5ld     nli     = %5ld\n", nni, nli);
  printf("nsetups = %5ld     netf    = %5ld\n", nsetups, netf);
  printf("npe     = %5ld     nps     = %5ld\n", npe, nps);
  printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);
}

/* Routine to send boundary data to neighboring PEs */
static void BSend(UserData data, realtype udata[])
{
  int      dsizex, dsizey;
  int      myproc, procB, procT, procL, procR;
  MPI_Comm comm;
  realtype *bufs;
  
  dsizex = data->dsizex;
  dsizey = data->dsizey;
  myproc = data->myproc;
  procB  = (data->isbottom) ? MPI_PROC_NULL : myproc-data->nprocsx;
  procT  = (data->istop   ) ? MPI_PROC_NULL : myproc+data->nprocsx;
  procL  = (data->isleft  ) ? MPI_PROC_NULL : myproc-1;
  procR  = (data->isright ) ? MPI_PROC_NULL : myproc+1;
  comm   = data->comm;
  bufs   = data->bufs;

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  int i, ly;
  sunindextype offsetu, offsetbuf;

  /* Populate buffers */
  // Note: y-lines (L/R bdry) are not stored contiguously in udata, so we load them into bufs
  if (!data->isleft) {
    for (ly = 0; ly < data->local_My; ly++) {
      offsetbuf = data->sendL + ly*NVARS;
      offsetu   = ly*dsizex;
      for (i = 0; i < NVARS; i++)
        bufs[offsetbuf+i] = udata[offsetu+i];
    }
  }
  if (!data->isright) {
    for (ly = 0; ly < data->local_My; ly++) {
      offsetbuf = data->sendR + ly*NVARS;
      offsetu   = (ly+1)*dsizex - NVARS;
      for (i = 0; i < NVARS; i++)
        bufs[offsetbuf+i] = udata[offsetu+i];
    }
  }

  /* Send x-lines directly from udata */
  offsetu = (data->local_My-1)*dsizex;
  MPI_Send(&udata[0],          dsizex, MPI_SUNREALTYPE, procB, 0, comm);
  MPI_Send(&udata[offsetu],    dsizex, MPI_SUNREALTYPE, procT, 0, comm);

#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  /* Send x-lines from bufs */
  MPI_Send(&bufs[data->sendB], dsizex, MPI_SUNREALTYPE, procB, 0, comm);
  MPI_Send(&bufs[data->sendT], dsizex, MPI_SUNREALTYPE, procT, 0, comm);
#endif

  /* Send y-lines from bufs */
  MPI_Send(&bufs[data->sendL], dsizey, MPI_SUNREALTYPE, procL, 0, comm);
  MPI_Send(&bufs[data->sendR], dsizey, MPI_SUNREALTYPE, procR, 0, comm);
}

/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*local_My realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */
static void BRecvPost(UserData data, MPI_Request request[])
{
  int      dsizex, dsizey;
  int      myproc, procB, procT, procL, procR;
  MPI_Comm comm;
  realtype *bufs;
  
  dsizex = data->dsizex;
  dsizey = data->dsizey;
  myproc = data->myproc;
  procB  = (data->isbottom) ? MPI_PROC_NULL : myproc-data->nprocsx;
  procT  = (data->istop   ) ? MPI_PROC_NULL : myproc+data->nprocsx;
  procL  = (data->isleft  ) ? MPI_PROC_NULL : myproc-1;
  procR  = (data->isright ) ? MPI_PROC_NULL : myproc+1;
  comm   = data->comm;
  bufs   = data->bufs;

  /* Receive x- and y-lines into bufs */
  MPI_Irecv(&bufs[data->recvB], dsizex, MPI_SUNREALTYPE, procB, 0, comm, &request[0]);
  MPI_Irecv(&bufs[data->recvT], dsizex, MPI_SUNREALTYPE, procT, 0, comm, &request[1]);
  MPI_Irecv(&bufs[data->recvL], dsizey, MPI_SUNREALTYPE, procL, 0, comm, &request[2]);
  MPI_Irecv(&bufs[data->recvR], dsizey, MPI_SUNREALTYPE, procR, 0, comm, &request[3]);
}

/* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*local_My realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */
static void BRecvWait(UserData data, MPI_Request request[])
{
  MPI_Status status[4];

  /* Wait for all requests  */
  MPI_Waitall(4,request,status);
}

/* ucomm routine.  This routine performs all communication
   between processors of data needed to calculate f. */
static void ucomm(UserData data, realtype t, N_Vector u)
{
  MPI_Request request[4];
  HYPRE_ParVector uhyp;
  realtype *udata;

  /* Get udata pointer */
  uhyp  = N_VGetVector_ParHyp(u);
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));

  /* Fill device send buffer and copy to host */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  unsigned sendlen  = data->bufslen/2; 
  unsigned block    = 256;
  unsigned grid     = (sendlen/NVARS + block - 1) / block;
  UserData data_dev = (UserData) data->data_dev;
  
  FillSendBufferKernel<<<grid,block>>>(data_dev, udata);

  NV_ADD_LANG_PREFIX_PH(Memcpy)(data->bufs,data->bufs_dev,sendlen*sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyDeviceToHost));
  CHECK_LAST_ERROR();
#endif

  /* Start receiving boundary data from neighboring PEs */
  BRecvPost(data, request);

  /* Send data from boundary of local grid to neighboring PEs */
  BSend(data, udata);

  /* Finish receiving boundary data from neighboring PEs */
  BRecvWait(data, request);

  /* Copy receive buffers to device (offset by sendlen) */
#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  NV_ADD_LANG_PREFIX_PH(Memcpy)(data->bufs_dev+sendlen,data->bufs+sendlen,sendlen*sizeof(realtype),NV_ADD_LANG_PREFIX_PH(MemcpyHostToDevice));
  CHECK_LAST_ERROR();
#endif
}

#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
/* fcalc routine. Compute f(t,y).  This routine assumes that communication
   between processors of data needed to calculate f has already been done,
   and this data is in the work array uext. */
static void fcalc(UserData data, realtype t, realtype udata[], realtype dudata[])
{
  int          i, lx, ly, jy;

  realtype     q3, q4, qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  realtype     c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt, c1rt, c2rt, cydn, cyup;
  realtype     hord1, hord2, horad1, horad2;
  sunindextype offsetbufs, offsetu;//, offsetue;
  
  int          dsizex, mybasey, local_Mx, local_My;
  realtype     dy, vdco, hdco, haco;
  int          recvB, recvT, recvL, recvR;
  realtype     *bufs;

  /* Make local copies of problem variables, for efficiency */
  dsizex   = data->dsizex;   mybasey  = data->mybasey;
  local_Mx = data->local_Mx; local_My = data->local_My;
  dy       = data->dy;       vdco     = data->vdco;
  hdco     = data->hdco;     haco     = data->haco;
  recvB    = data->recvB;    recvT    = data->recvT;
  recvL    = data->recvL;    recvR    = data->recvR;
  bufs     = data->bufs;

  /* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary proc, copy data from the first interior mesh line of u to uext */

  /* If myprocy = 0, copy x-line 2 of u to uext */
  if (data->isbottom) {
    offsetu = dsizex;
    for (i = 0; i < dsizex; i++) bufs[recvB+i] = udata[offsetu+i];
  }
  /* If myprocy = nprocsy-1, copy x-line local_My-1 of u to uext */
  if (data->istop) {
    offsetu = (local_My-2)*dsizex;
    for (i = 0; i < dsizex; i++) bufs[recvT+i] = udata[offsetu+i];
  }
  /* If myprocx = 0, copy y-line 2 of u to uext */
  if (data->isleft) {
    for (ly = 0; ly < local_My; ly++) {
      offsetbufs = recvL + ly*NVARS;
      offsetu    = ly*dsizex + NVARS;
      for (i = 0; i < NVARS; i++) bufs[offsetbufs+i] = udata[offsetu+i];
    }
  }
  /* If myprocx = nprocsx-1, copy y-line local_Mx-1 of u to uext */
  if (data->isright) {
    for (ly = 0; ly < local_My; ly++) {
      offsetbufs = recvR + ly*NVARS;
      offsetu    = (ly+1)*dsizex - 2*NVARS;
      for (i = 0; i < NVARS; i++) bufs[offsetbufs+i] = udata[offsetu+i];
    }
  }

  /* Set diurnal rate coefficients as functions of t, and save q4 in
  data block for use by preconditioner evaluation routine */
  s = sin((data->om)*t);
  if (s > RCONST(0.0)) {
    q3 = SUNRexp(-A3/s);
    q4 = SUNRexp(-A4/s);
  } else {
    q3 = RCONST(0.0);
    q4 = RCONST(0.0);
  }
  data->q4 = q4;

  /* Loop over all grid points in local subgrid */
  for (ly = 0; ly < local_My; ly++) {

    jy   = ly + mybasey;

    /* Set vertical diffusion coefficients at jy +- 1/2 */
    ydn  = YMIN + (jy - RCONST(0.5))*dy;
    yup  = ydn + dy;
    cydn = vdco*SUNRexp(RCONST(0.2)*ydn);
    cyup = vdco*SUNRexp(RCONST(0.2)*yup);
    for (lx = 0; lx < local_Mx; lx++) {

      /* Extract c1 and c2, and set kinetic rate terms */
      // offsetue = (lx+1)*NVARS + (ly+1)*dsizex2;
      offsetu  = ly*dsizex + lx;
      c1       = udata[offsetu];
      c2       = udata[offsetu+1];
      qq1      = Q1*c1*C3;
      qq2      = Q2*c1*c2;
      qq3      = q3*C3;
      qq4      = q4*c2;
      rkin1    = -qq1 - qq2 + RCONST(2.0)*qq3 + qq4;
      rkin2    = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms */
      c1dn     = (ly==0         ) ? bufs[recvB+lx*NVARS  ] : udata[offsetu-dsizex  ];
      c2dn     = (ly==0         ) ? bufs[recvB+lx*NVARS+1] : udata[offsetu-dsizex+1];
      c1up     = (ly==local_My-1) ? bufs[recvT+lx*NVARS  ] : udata[offsetu+dsizex  ];
      c2up     = (ly==local_My-1) ? bufs[recvT+lx*NVARS+1] : udata[offsetu+dsizex+1];
      vertd1   = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      vertd2   = cyup*(c2up - c2) - cydn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms */
      c1lt     = (lx==0         ) ? bufs[recvL+ly*NVARS  ] : udata[offsetu-NVARS  ];
      c2lt     = (lx==0         ) ? bufs[recvL+ly*NVARS+1] : udata[offsetu-NVARS+1];
      c1rt     = (lx==local_Mx-1) ? bufs[recvR+ly*NVARS  ] : udata[offsetu+NVARS  ];
      c2rt     = (lx==local_Mx-1) ? bufs[recvR+ly*NVARS+1] : udata[offsetu+NVARS+1];
      hord1    = hdco*(c1rt - RCONST(2.0)*c1 + c1lt);
      hord2    = hdco*(c2rt - RCONST(2.0)*c2 + c2lt);
      horad1   = haco*(c1rt - c1lt);
      horad2   = haco*(c2rt - c2lt);

      /* Load all terms into dudata */
      offsetu  = lx*NVARS + ly*dsizex;
      dudata[offsetu]   = vertd1 + hord1 + horad1 + rkin1;
      dudata[offsetu+1] = vertd2 + hord2 + horad2 + rkin2;
    }
  }
}
#endif // end SERIAL-only fcalc

/****** Kernels and device helper functions (for CUDA/HIP backend) ******/

#if defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
__global__ void fcalcKernel(UserData data_dev, realtype t, realtype *udata, realtype *udotdata)
{
  sunindextype tid, lx, ly, gy; // lx,ly: local grid point coords; gy: global grid point y coord

  realtype     q3, q4, qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  realtype     c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt, c1rt, c2rt, cydn, cyup;
  realtype     hord1, hord2, horad1, horad2;
  // sunindextype offsetu, offsetue;
  
  int          local_Mx, local_My;
  realtype     dy, vdco, hdco, haco;
  int          recvB, recvT, recvL, recvR;
  realtype     *bufs_dev;

  /* Make local copies of problem variables, for efficiency */
  local_Mx = data_dev->local_Mx; local_My = data_dev->local_My;
  dy       = data_dev->dy;       vdco     = data_dev->vdco;
  hdco     = data_dev->hdco;     haco     = data_dev->haco;
  recvB    = data_dev->recvB;    recvT    = data_dev->recvT;
  recvL    = data_dev->recvL;    recvR    = data_dev->recvR;
  bufs_dev = data_dev->bufs_dev;

  /* Get tid and calculate udot */
  tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < local_Mx*local_My)
  {
    ly   = tid/local_Mx; // Conceptually indexed L->R, B->T: ┌ 2 3 ┐
    lx   = tid%local_Mx; //                                  └ 0 1 ┘
    gy   = ly + data_dev->mybasey;

    /* Get c1&c2 values from udata (take from Recv buffers as needed) */
    // Note: [homogeneous Neumann BC] if a boundary proc, pull from interior instead of buffer
    c1   = udata[NVARS*tid];
    c2   = udata[NVARS*tid+1];
    //   thread on proc bdry?  proc on domain bdry?   [Y,Y] pull from interior        [Y,N] pull from buffer        [N,Y/N] pull directly  
    c1dn = (ly==0         ) ? ((data_dev->isbottom) ? udata[NVARS*(tid+local_Mx)  ] : bufs_dev[recvB+NVARS*lx  ]) : udata[NVARS*(tid-local_Mx)  ];
    c2dn = (ly==0         ) ? ((data_dev->isbottom) ? udata[NVARS*(tid+local_Mx)+1] : bufs_dev[recvB+NVARS*lx+1]) : udata[NVARS*(tid-local_Mx)+1];
    c1up = (ly==local_My-1) ? ((data_dev->istop   ) ? udata[NVARS*(tid-local_Mx)  ] : bufs_dev[recvT+NVARS*lx  ]) : udata[NVARS*(tid+local_Mx)  ];
    c2up = (ly==local_My-1) ? ((data_dev->istop   ) ? udata[NVARS*(tid-local_Mx)+1] : bufs_dev[recvT+NVARS*lx+1]) : udata[NVARS*(tid+local_Mx)+1];
    c1lt = (lx==0         ) ? ((data_dev->isleft  ) ? udata[NVARS*(tid+1)         ] : bufs_dev[recvL+NVARS*ly  ]) : udata[NVARS*(tid-1)         ];
    c2lt = (lx==0         ) ? ((data_dev->isleft  ) ? udata[NVARS*(tid+1)+1       ] : bufs_dev[recvL+NVARS*ly+1]) : udata[NVARS*(tid-1)+1       ];
    c1rt = (lx==local_Mx-1) ? ((data_dev->isright ) ? udata[NVARS*(tid-1)         ] : bufs_dev[recvR+NVARS*ly  ]) : udata[NVARS*(tid+1)         ];
    c2rt = (lx==local_Mx-1) ? ((data_dev->isright ) ? udata[NVARS*(tid-1)+1       ] : bufs_dev[recvR+NVARS*ly+1]) : udata[NVARS*(tid+1)+1       ];

    /* Set diurnal rate coefficients as functions of t */
    s = sin((data_dev->om)*t);
    if (s > RCONST(0.0)) {
      q3 = SUNRexp(-A3/s);
      q4 = SUNRexp(-A4/s);
    } else {
      q3 = RCONST(0.0);
      q4 = RCONST(0.0);
    }
    // data_dev->q4 = q4; // Only need to export q4 if we later implement a preconditioner for GPU

    /* Set vertical diffusion coefficients at jy +- 1/2 */
    ydn    = YMIN + (gy - RCONST(0.5))*dy;
    yup    = ydn + dy;
    cydn   = vdco*SUNRexp(RCONST(0.2)*ydn);
    cyup   = vdco*SUNRexp(RCONST(0.2)*yup);

    /* Set kinetic rate terms */
    qq1    = Q1*c1*C3;
    qq2    = Q2*c1*c2;
    qq3    = q3*C3;
    qq4    = q4*c2;
    rkin1  = -qq1 - qq2 + RCONST(2.0)*qq3 + qq4;
    rkin2  = qq1 - qq2 - qq4;

    /* Set vertical diffusion terms */
    vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
    vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);

    /* Set horizontal diffusion and advection terms */
    hord1  = hdco*(c1rt - RCONST(2.0)*c1 + c1lt);
    hord2  = hdco*(c2rt - RCONST(2.0)*c2 + c2lt);
    horad1 = haco*(c1rt - c1lt);
    horad2 = haco*(c2rt - c2lt);

    /* Load all terms into dudata */
    udotdata[NVARS*tid]   = vertd1 + hord1 + horad1 + rkin1;
    udotdata[NVARS*tid+1] = vertd2 + hord2 + horad2 + rkin2;
  }
  return;
}

__global__ void FillSendBufferKernel(UserData data_dev, realtype *udata)
{
  sunindextype tid, lx, ly, i;
  int offsetB, offsetT, offsetL, offsetR;
  int local_Mx, local_My;
  realtype *bufs_dev;

  offsetB  = data_dev->sendB/NVARS; offsetT  = data_dev->sendT/NVARS;
  offsetL  = data_dev->sendL/NVARS; offsetR  = data_dev->sendR/NVARS;
  local_Mx = data_dev->local_Mx;    local_My = data_dev->local_My;
  bufs_dev = data_dev->bufs_dev;

  tid = blockDim.x * blockIdx.x + threadIdx.x;
  if      (tid < offsetB) { ly = 0;           lx = tid-offsetB; }
  else if (tid < offsetT) { ly = local_My-1;  lx = tid-offsetT; }
  else if (tid < offsetL) { ly = tid-offsetL; lx = 0;           }
  else if (tid < offsetR) { ly = tid-offsetR; lx = local_Mx-1;  }
  else return;

  for (i = 0; i < NVARS; i++) {
    bufs_dev[NVARS*tid+i]   = udata[NVARS*(ly*local_Mx + lx)+i];
  }

  return;
}
#endif // end GPU-only kernels


/***************** Functions Called by the Solver *************************/

/* f routine.  Evaluate f(t,y).  First call ucomm to do communication of
   subgrid boundary data into uext.  Then calculate f by a call to fcalc. */
static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype *udata, *udotdata;
  UserData data;

  HYPRE_ParVector uhyp;
  HYPRE_ParVector udothyp;

  /* Extract hypre vectors */
  uhyp     = N_VGetVector_ParHyp(u);
  udothyp  = N_VGetVector_ParHyp(udot);

  /* Access hypre vectors local data */
  udata    = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));
  udotdata = hypre_VectorData(hypre_ParVectorLocalVector(udothyp));

  /* Extract UserData */
  data     = (UserData) user_data;

  /* Call ucomm to do inter-processor communication */
  ucomm(data, t, u);

  /* Call fcalc to calculate all right-hand sides */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
  fcalc(data, t, udata, udotdata);
#elif defined(SUNDIALS_HYPRE_BACKENDS_CUDA_OR_HIP)
  unsigned block    = 256;
  unsigned grid     = (data->local_M + block - 1) / block;
  UserData data_dev = (UserData) data->data_dev;

  fcalcKernel<<<grid,block>>>(data_dev, t, udata, udotdata);
#endif

  return(0);
}

/* Preconditioner setup routine. Generate and preprocess P. */
#if defined(SUNDIALS_HYPRE_BACKENDS_SERIAL)
static int Precond(realtype tn, N_Vector u, N_Vector fu,
                   booleantype jok, booleantype *jcurPtr,
                   realtype gamma, void *user_data)
{
  int             lx, ly, jy, ier, offset;
  realtype        c1, c2, cydn, cyup, diag, ydn, yup;
  realtype        *udata, **a, **j;
  HYPRE_ParVector uhyp;

  UserData        data;
  realtype        q4, dy, vdco, hdco;

  /* Make local copies of UserData variables for legibility */
  data  = (UserData) user_data;
  q4    = data->q4;   dy   = data->dy;
  vdco  = data->vdco; hdco = data->hdco;

  /* Get hypre vector data pointer */
  uhyp  = N_VGetVector_ParHyp(u);
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));

  /* jok = SUNTRUE: Copy Jbd to P */
  if (jok) {
    for (ly = 0; ly < data->local_My; ly++)
      for (lx = 0; lx < data->local_Mx; lx++)
        SUNDlsMat_denseCopy(data->Jbd[lx][ly], data->P[lx][ly], NVARS, NVARS);

  *jcurPtr = SUNFALSE;

  }

  /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */
  else {
    /* Compute 2x2 diagonal Jacobian blocks (using q4 values
     computed on the last f call).  Load into P. */
    for (ly = 0; ly < data->local_My; ly++) {
      jy   = ly + data->myprocy*data->local_My;
      ydn  = YMIN + (jy - RCONST(0.5))*dy;
      yup  = ydn + dy;
      cydn = vdco*SUNRexp(RCONST(0.2)*ydn);
      cyup = vdco*SUNRexp(RCONST(0.2)*yup);
      diag = -(cydn + cyup + RCONST(2.0)*hdco); 
      for (lx = 0; lx < data->local_Mx; lx++) {
        offset = lx*NVARS + ly*data->dsizex;
        c1     = udata[offset];
        c2     = udata[offset+1];
        j      = data->Jbd[lx][ly];
        a      = data->P[lx][ly];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) =  -Q2*c1 + q4;
        IJth(j,2,1) =   Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4)    + diag;
        SUNDlsMat_denseCopy(j, a, NVARS, NVARS);
      }
    }

    *jcurPtr = SUNTRUE;

  }

  /* Scale by -gamma */
  for (ly = 0; ly < data->local_My; ly++)
    for (lx = 0; lx < data->local_Mx; lx++)
      SUNDlsMat_denseScale(-gamma, data->P[lx][ly], NVARS, NVARS);

  /* Add identity matrix and do LU decompositions on blocks in place */
  for (lx = 0; lx < data->local_Mx; lx++) {
    for (ly = 0; ly < data->local_My; ly++) {
      SUNDlsMat_denseAddIdentity(data->P[lx][ly], NVARS);
      ier = SUNDlsMat_denseGETRF(data->P[lx][ly], NVARS, NVARS, data->pivot[lx][ly]);
      if (ier != 0) return(1);
    }
  }

  return(0);
}

/* Preconditioner solve routine */
static int PSolve(realtype tn, N_Vector u, N_Vector fu,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data)
{
  int lx, ly;
  realtype *zdata, *v;
  HYPRE_ParVector zhyp;

  UserData data;

  /* Get UserData pointer */
  data = (UserData) user_data;

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z.
     First copy vector r to z. */
  N_VScale(RCONST(1.0), r, z);
  zhyp   = N_VGetVector_ParHyp(z); /* extract hypre vector */
  zdata  = hypre_VectorData(hypre_ParVectorLocalVector(zhyp));

  for (lx = 0; lx < data->local_Mx; lx++) {
    for (ly = 0; ly < data->local_My; ly++) {
      v = &(zdata[lx*NVARS + ly*data->dsizex]);
      SUNDlsMat_denseGETRS(data->P[lx][ly], NVARS, data->pivot[lx][ly], v);
    }
  }

  return(0);
}
#endif // end SERIAL-only preconditioning

/*********** Private function to check function return values ***********/

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_flag(void *flagvalue, const char *funcname, int opt, int id)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  } else if (opt == 1) { /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with flag = %d\n\n",
              id, funcname, *errflag);
      return(1);
    }
  } else if (opt == 2 && flagvalue == NULL) { /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  }

  return(0);
}