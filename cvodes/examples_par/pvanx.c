/************************************************************************
 *                                                                      *
 * File       : pvanx.c                                                 *
 * Programmers: Radu Serban @ LLNL                                      *
 * Version of : 30 March 2003                                           *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * The following is a simple example problem, with the program for its  *
 * solution by CVODE.  The problem is the semi-discrete form of the     *
 * advection-diffusion equation in 1-D:                                 *
 *   du/dt = p1 * d^2u / dx^2 + p2 * du / dx                            *
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.      *
 * Homogeneous Dirichlet boundary conditions are posed, and the         *
 * initial condition is                                                 *
 *   u(x,t=0) = x(2-x)exp(2x) .                                         *
 * The nominal values of the two parameters are: p1=1.0, p2=0.5         *
 * The PDE is discretized on a uniform grid of size MX+2 with           *
 * central differencing, and with boundary values eliminated,           *
 * leaving an ODE system of size NEQ = MX.                              *
 * This program solves the problem with the option for nonstiff systems:*
 * ADAMS method and functional iteration.                               *
 * It uses scalar relative and absolute tolerances.                     *
 *                                                                      *
 * In addition to the solution, sensitivities with respect to p1 and p2 *
 * as well as with respect to initial conditions are computed for the   *
 * quantity:                                                            *
 *    g(t, u, p) = int_x u(x,t) at t = 5                                *
 * These sensitivities are obtained by solving the adjoint system:      *
 *    dv/dt = -p1 * d^2 v / dx^2 + p2 * dv / dx                         *
 * with homogeneous Ditrichlet boundary conditions and the final        *
 * condition                                                            *
 *    v(x,t=5) = 1.0                                                    *
 * Then, v(x, t=0) represents the sensitivity of g(5) with respect to   *
 * u(x, t=0) and the gradient of g(5) with respect to p1, p2 is         *
 *    (dg/dp)^T = [  int_t int_x (v * d^2u / dx^2) dx dt ]              *
 *                [  int_t int_x (v * du / dx) dx dt     ]              *
 *                                                                      *
 * This version uses MPI for user routines.                             *
 * Execute with Number of Processors = N,  with 1 <= N <= MX.           *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"
#include "cvodea.h"
#include "nvector_parallel.h"
#include "mpi.h"

/* Problem Constants */
#define XMAX  2.0          /* domain boundary            */
#define MX    20           /* mesh dimension             */
#define NEQ   MX           /* number of equations        */
#define ATOL  1.e-5        /* scalar absolute tolerance  */
#define T0    0.0          /* initial time               */
#define TOUT  2.5          /* output time increment      */

/* Adjoint Problem Constants */
#define NP    2            /* number of parameters       */
#define STEPS 100          /* steps between check points */

/* Type : UserData */
typedef struct {
  realtype p[2];            /* model parameters                         */
  realtype dx;              /* spatial discretization grid              */
  realtype hdcoef, hacoef;  /* diffusion and advection coefficients     */
  integertype npes, my_pe;  /* total number of processes and current ID */
  MPI_Comm comm;        /* MPI communicator                         */
  realtype *z1, *z2;        /* work space                               */
} *UserData;


/* Private Helper Functions */
static void SetIC(N_Vector u, realtype dx, integertype my_length, integertype my_base);
static void SetICback(N_Vector uB, integertype my_base);
static realtype Xintgr(realtype *z, integertype l, realtype dx);
static realtype Compute_g(N_Vector u, UserData data);

/* Functions Called by the CVODES Solver */
static void f(realtype t, N_Vector u, N_Vector udot, void *f_data);
static void fB(realtype t, N_Vector u, 
               N_Vector uB, N_Vector uBdot, void *f_dataB);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  M_Env machEnvF, machEnvB;

  UserData data;

  void *cvadj_mem;
  void *cvode_mem;
  
  long int iopt[OPT_SIZE];
  realtype ropt[OPT_SIZE];
  N_Vector u;
  realtype reltol, abstol;

  N_Vector uB;
  realtype *uBdata;

  realtype dx, t, umax, g_val;
  int flag, my_pe, nprocs, npes, ncheck;
  integertype local_N, nperpe, nrem, my_base, i, iglobal;

  MPI_Comm comm;

  /* Initialize MPI and get total number of pe's, and my_pe. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &my_pe);

  npes = nprocs - 1; /* pe's dedicated to PDE integration */

  /* Allocate and load user data structure */
  data = (UserData) malloc(sizeof *data);
  data->p[0] = 1.0;
  data->p[1] = 0.5;
  dx = data->dx = XMAX/((realtype)(MX+1));
  data->hdcoef = data->p[0]/(dx*dx);
  data->hacoef = data->p[1]/(2.0*dx);
  data->comm = comm;
  data->npes = npes;
  data->my_pe = my_pe;

  /* Set local vector length. */
  if (my_pe < npes) {
    nperpe = NEQ/npes;
    nrem = NEQ - npes*nperpe;
    local_N = (my_pe < nrem) ? nperpe+1 : nperpe;
    my_base = (my_pe < nrem) ? my_pe*local_N : my_pe*nperpe + nrem;
  }

  /*------------------------- 
    FORWARD INTEGRATION PHASE 
  ---------------------------*/

  /* Make last process inactive for forward phase */
  if (my_pe == npes) local_N = 0;

  /* Initialize machine environment for forward phase */
  machEnvF = M_EnvInit_Parallel(comm, local_N, NEQ, &argc, &argv);
  if (machEnvF == NULL) {
    if(my_pe == 0) printf("M_EnvInit_Parallel failed.\n"); 
    return(1);
  }

  /* Set relative and absolute tolerances for forward phase */
  reltol = 0.0;
  abstol = ATOL;

  /* Allocate and initialize forward variables */
  u = N_VNew(machEnvF);
  SetIC(u, dx, local_N, my_base);

  /* Allocate CVODES memory for forward integration */
  cvode_mem = CVodeMalloc(f, T0, u, ADAMS, FUNCTIONAL, SS, &reltol,
                          &abstol, data, NULL, FALSE, iopt, ropt, machEnvF);
  if (cvode_mem == NULL) {
    if(my_pe == 0) printf("CVodeMalloc failed.\n"); 
    return(1);
  }

  /* Allocate combined forward/backward memory */
  cvadj_mem = CVadjMalloc(cvode_mem, STEPS);

  /* Integrate to TOUT and collect check point information */
  flag = CVodeF(cvadj_mem, TOUT, u, &t, NORMAL, &ncheck);
  if (flag != SUCCESS) {
    if(my_pe == 0) printf("CVode failed, flag=%d.\n", flag);
    return(1);
  }

  /* Compute and print value of g(5) */
  g_val = Compute_g(u, data);
  if(my_pe == npes) {
    printf("\n (PE# %d)\n", my_pe);
    printf("    g(t1) = %g\n", g_val);
  }

  /*-------------------------- 
    BACKWARD INTEGRATION PHASE 
  ----------------------------*/

  /* Allocate work space */
  data->z1 = (realtype *)malloc(local_N*sizeof(realtype));
  data->z2 = (realtype *)malloc(local_N*sizeof(realtype));

  /* Activate last process for integration of the quadrature equations */
  if(my_pe == npes) local_N = NP;

  /* Initialize machine environment for backward phase */
  machEnvB = M_EnvInit_Parallel(comm, local_N, NEQ+NP, &argc, &argv);
  if (machEnvB == NULL) {
    if(my_pe == 0) printf("M_EnvInit_Parallel failed.\n"); 
    return(1);
  }

  /* Allocate and initialize backward variables */
  uB = N_VNew(machEnvB);
  SetICback(uB, my_base);

  /* Allocate CVODES memory for the backward integration */
  flag = CVodeMallocB(cvadj_mem, fB, TOUT, uB, ADAMS, FUNCTIONAL, SS, &reltol, 
                      &abstol, data, NULL, FALSE, NULL, NULL, machEnvB);
  if (flag != SUCCESS) { 
    if(my_pe == 0) printf("CVodeMallocB failed, flag=%d.\n", flag);
    return(1);
  }

  /* Integrate to T0 */
  flag = CVodeB(cvadj_mem, uB);
  if (flag < 0) { 
    if(my_pe == 0) printf("CVodeB failed, flag=%d.\n", flag);
    return(1);
  }

  /* Print results (adjoint states and quadrature variables) */
  uBdata = NV_DATA_P(uB);
  printf("\n (PE# %d)\n", my_pe);
  if (my_pe == npes) {
    printf("    dgdp(t1) = [ %g  %g ]\n", -uBdata[0], -uBdata[1]);
  } else {
    for (i=1; i<=local_N; i++) {
      iglobal = my_base + i;
      printf("    mu(t0)[%2d] = %g\n", iglobal, uBdata[i-1]);
    }
  }

  /* Clean-Up */
  N_VFree(u);                                  /* forward variables     */
  N_VFree(uB);                                 /* backward variables    */
  CVodeFree(cvode_mem);                        /* CVODES memory block   */    
  CVadjFree(cvadj_mem);                        /* combined memory block */
  free(data->z1); free(data->z2); free(data);  /* user data structure   */
  M_EnvFree_Parallel(machEnvF);                /* forward M_Env         */
  M_EnvFree_Parallel(machEnvB);                /* backward M_Env        */
  
  /* Finalize MPI */
  MPI_Finalize();

  return(0);
}


/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */
static void SetIC(N_Vector u, realtype dx, integertype my_length, integertype my_base)
{
  int i;
  integertype iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u */
  udata = NV_DATA_P(u);
  my_length = NV_LOCLENGTH_P(u);

  /* Load initial profile into u vector */
  for (i=1; i<=my_length; i++) {
    iglobal = my_base + i;
    x = iglobal*dx;
    udata[i-1] = x*(XMAX - x)*exp(2.0*x);
  }  
}

/* Set final conditions in uB vector */
static void SetICback(N_Vector uB, integertype my_base)
{
  int i;
  realtype *uBdata;
  integertype my_length;

  /* Set pointer to data array and get local length of uB */
  uBdata = NV_DATA_P(uB);
  my_length = NV_LOCLENGTH_P(uB);

  /* Set adjoint states to 1.0 and quadrature variables to 0.0 */
  if(my_base == -1) for (i=0; i<my_length; i++) uBdata[i] = 0.0;
  else              for (i=0; i<my_length; i++) uBdata[i] = 1.0;
}

/* Compute local value of the space integral int_x z(x) dx */
static realtype Xintgr(realtype *z, integertype l, realtype dx)
{
  realtype my_intgr;
  integertype i;

  my_intgr = 0.5*(z[0] + z[l-1]);
  for (i = 1; i < l-1; i++)
    my_intgr += z[i]; 
  my_intgr *= dx;

  return(my_intgr);
}

/* Compute value of g(u) */
static realtype Compute_g(N_Vector u, UserData data)
{
  realtype intgr, my_intgr, dx, *udata;
  integertype my_length;
  int npes, my_pe, i;
  MPI_Status status;
  MPI_Comm comm;

  /* Extract MPI info. from data */
  comm = data->comm;
  npes = data->npes;
  my_pe = data->my_pe;

  dx = data->dx;

  if (my_pe == npes) {  /* Loop over all other processes and sum */
    intgr = 0.0;
    for (i=0; i<npes; i++) {
      MPI_Recv(&my_intgr, 1, PVEC_REAL_MPI_TYPE, i, 0, comm, &status); 
      intgr += my_intgr;
    }
    return(intgr);
  } else {              /* Compute local portion of the integral */
    udata = NV_DATA_P(u);
    my_length = NV_LOCLENGTH_P(u);
    my_intgr = Xintgr(udata, my_length, dx);
    MPI_Send(&my_intgr, 1, PVEC_REAL_MPI_TYPE, npes, 0, comm);
    return(my_intgr);
  }
}

/***************** Function Called by the CVODE Solver ******************/

/* f routine. Compute f(t,u) for forward phase. */
static void f(realtype t, N_Vector u, N_Vector udot, void *f_data)
{
  realtype uLeft, uRight, ui, ult, urt;
  realtype hordc, horac, hdiff, hadv;
  realtype *udata, *dudata;
  integertype i, my_length;
  int npes, my_pe, my_pe_m1, my_pe_p1, last_pe, my_last;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  /* Extract MPI info. from data */
  data = (UserData) f_data;
  comm = data->comm;
  npes = data->npes;
  my_pe = data->my_pe;
  
  /* If this process is inactive, return now */
  if (my_pe == npes) return;

  /* Extract problem constants from data */
  hordc = data->hdcoef;
  horac = data->hacoef;

  /* Find related processes */
  my_pe_m1 = my_pe - 1;
  my_pe_p1 = my_pe + 1;
  last_pe = npes - 1;

  /* Obtain local arrays */
  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);
  my_length = NV_LOCLENGTH_P(u);
  my_last = my_length - 1;

  /* Pass needed data to processes before and after current process. */
   if (my_pe != 0)
     MPI_Send(&udata[0], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
   if (my_pe != last_pe)
     MPI_Send(&udata[my_length-1], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);   

  /* Receive needed data from processes before and after current process. */
   if (my_pe != 0)
     MPI_Recv(&uLeft, 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
   else uLeft = 0.0;
   if (my_pe != last_pe)
     MPI_Recv(&uRight, 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm,
              &status);   
   else uRight = 0.0;

  /* Loop over all grid points in current process. */
  for (i=0; i<my_length; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = udata[i];
    ult = (i==0) ? uLeft: udata[i-1];
    urt = (i==my_length-1) ? uRight : udata[i+1];

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - 2.0*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i] = hdiff + hadv;
  }
}

/* fB routine. Compute right hand side of backward problem */
static void fB(realtype t, N_Vector u, 
               N_Vector uB, N_Vector uBdot, void *f_dataB)
{
  realtype *uBdata, *duBdata, *udata, *zB;
  realtype uBLeft, uBRight, uBi, uBlt, uBrt;
  realtype uLeft, uRight, ui, ult, urt;
  realtype dx, hordc, horac, hdiff, hadv;
  realtype *z1, *z2, intgr1, intgr2;
  integertype i, my_length;
  int npes, my_pe, my_pe_m1, my_pe_p1, last_pe, my_last;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  /* Extract MPI info. from data */
  data = (UserData) f_dataB;
  comm = data->comm;
  npes = data->npes;
  my_pe = data->my_pe;

  if (my_pe == npes) { /* This process performs the quadratures */

    /* Obtain local arrays */
    duBdata = NV_DATA_P(uBdot);
    my_length = NV_LOCLENGTH_P(uB);

    /* Loop over all other processes and load right hand side of quadrature eqs. */
    duBdata[0] = 0.0;
    duBdata[1] = 0.0;
    for (i=0; i<npes; i++) {
      MPI_Recv(&intgr1, 1, PVEC_REAL_MPI_TYPE, i, 0, comm, &status); 
      duBdata[0] += intgr1;
      MPI_Recv(&intgr2, 1, PVEC_REAL_MPI_TYPE, i, 0, comm, &status); 
      duBdata[1] += intgr2;
    }

  } else { /* This process integrates part of the PDE */

    /* Extract problem constants and work arrays from data */
    dx    = data->dx;
    hordc = data->hdcoef;
    horac = data->hacoef;
    z1    = data->z1;
    z2    = data->z2;

    /* Compute related parameters. */
    my_pe_m1 = my_pe - 1;
    my_pe_p1 = my_pe + 1;
    last_pe  = npes - 1;
    my_last  = my_length - 1;
    
    /* Obtain local arrays */
    uBdata = NV_DATA_P(uB);
    duBdata = NV_DATA_P(uBdot);
    udata = NV_DATA_P(u);
    my_length = NV_LOCLENGTH_P(uB);

    /* Pass needed data to processes before and after current process. */
    if (my_pe != 0) {
      MPI_Send(&udata[0], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
      MPI_Send(&uBdata[0], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
    }
    if (my_pe != last_pe) {
       MPI_Send(&udata[my_length-1], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);
       MPI_Send(&uBdata[my_length-1], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);
    }
    
    /* Receive needed data from processes before and after current process. */
    if (my_pe != 0) {
      MPI_Recv(&uLeft, 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
      MPI_Recv(&uBLeft, 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
    } else {
      uLeft = 0.0;
      uBLeft = 0.0;
    }
    if (my_pe != last_pe) {
      MPI_Recv(&uRight, 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm, &status);
      MPI_Recv(&uBRight, 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm, &status);
    } else {
      uRight = 0.0;
      uBRight = 0.0;
    }

    /* Loop over all grid points in current process. */
    for (i=0; i<my_length; i++) {
      
      /* Extract uB at x_i and two neighboring points */
      uBi = uBdata[i];
      uBlt = (i==0) ? uBLeft: uBdata[i-1];
      uBrt = (i==my_length-1) ? uBRight : uBdata[i+1];
      
      /* Set diffusion and advection terms and load into udot */
      hdiff = hordc*(uBlt - 2.0*uBi + uBrt);
      hadv = horac*(uBrt - uBlt);
      duBdata[i] = - hdiff + hadv;

      /* Extract u at x_i and two neighboring points */
      ui = udata[i];
      ult = (i==0) ? uLeft: udata[i-1];
      urt = (i==my_length-1) ? uRight : udata[i+1];

      /* Load integrands of the two space integrals */
      z1[i] = uBdata[i]*(ult - 2.0*ui + urt)/(dx*dx);
      z2[i] = uBdata[i]*(urt - ult)/(2.0*dx);
    }

    /* Compute local integrals */
    intgr1 = Xintgr(z1, my_length, dx);
    intgr2 = Xintgr(z2, my_length, dx);

    /* Send local integrals to 'quadrature' process */
    MPI_Send(&intgr1, 1, PVEC_REAL_MPI_TYPE, npes, 0, comm);
    MPI_Send(&intgr2, 1, PVEC_REAL_MPI_TYPE, npes, 0, comm);

  }

}
