/************************************************************************
 * File       : pvnx.c                                                  *
 * Programmers: Scott D. Cohen, Alan C. Hindmarsh, George Byrne, and    *
 *              Radu Serban @LLNL                                       *
 * Version of : 11 July 2003                                            *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * The following is a simple example problem, with the program for its  *
 * solution by CVODE.  The problem is the semi-discrete form of the     *
 * advection-diffusion equation in 1-D:                                 *
 *   du/dt = d^2 u / dx^2 + .5 du/dx                                    *
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.      *
 * Homogeneous Dirichlet boundary conditions are posed, and the         *
 * initial condition is                                                 *
 *   u(x,t=0) = x(2-x)exp(2x) .                                         *
 * The PDE is discretized on a uniform grid of size MX+2 with           *
 * central differencing, and with boundary values eliminated,           *
 * leaving an ODE system of size NEQ = MX.                              *
 * This program solves the problem with the option for nonstiff systems:*
 * ADAMS method and functional iteration.                               *
 * It uses scalar relative and absolute tolerances.                     *
 * Output is printed at t = .5, 1.0, ..., 5.                            *
 * Run statistics (optional outputs) are printed at the end.            *
 *                                                                      *
 * This version uses MPI for user routines, and the MPI_PVODE solver.   *
 * Execute with Number of Processors = N,  with 1 <= N <= MX.           *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* CVODE header files with a description of contents used here */

#include "sundialstypes.h"     /* definitions of realtype, integertype        */
#include "cvode.h"             /* prototypes for CVodeMalloc, CVode, and      */
                               /* CVodeFree, constants OPT_SIZE, FUNCTIONAL,  */
                               /* ADAMS, SS, SUCCESS, NST,NFE, NNI, NCFN,NETF */
#include "nvector_parallel.h"  /* definitions of type N_Vector and vector     */
                               /* macros, prototypes for N_Vector functions   */
#include "mpi.h"               /* for MPI constants and types                 */


/* Problem Constants */

#define XMAX  2.0          /* domain boundary           */
#define MX    10           /* mesh dimension            */
#define NEQ   MX           /* number of equations       */
#define ATOL  1.e-5        /* scalar absolute tolerance */
#define T0    0.0          /* initial time              */
#define T1    0.5          /* first output time         */
#define DTOUT 0.5          /* output time increment     */
#define NOUT  10           /* number of output times    */


/* Type : UserData 
   contains grid constants, parallel machine parameters, work array. */

typedef struct {
  realtype dx, hdcoef, hacoef;
  integertype npes, my_pe;
  MPI_Comm comm;
  realtype z[100];
} *UserData;


/* Private Helper Functions */

static void SetIC(N_Vector u, realtype dx, integertype my_length,
                  integertype my_base);

static void PrintFinalStats(void *cvode_mem);

/* Functions Called by the CVODE Solver */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  NV_Spec nvSpec;
  realtype dx, reltol, abstol, t, tout, umax;
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, flag, my_pe, npes, nst;
  integertype local_N, nperpe, nrem, my_base;

  MPI_Comm comm;

  /* Get processor number, total number of pe's, and my_pe. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  /* Set local vector length. */
  nperpe = NEQ/npes;
  nrem = NEQ - npes*nperpe;
  local_N = (my_pe < nrem) ? nperpe+1 : nperpe;
  my_base = (my_pe < nrem) ? my_pe*local_N : my_pe*nperpe + nrem;

  data = (UserData) malloc(sizeof *data);  /* Allocate data memory */

  data->comm = comm;
  data->npes = npes;
  data->my_pe = my_pe;

  nvSpec = NV_SpecInit_Parallel(comm, local_N, NEQ, &argc, &argv);

  u = N_VNew(nvSpec);    /* Allocate u vector */

  reltol = 0.0;           /* Set the tolerances */
  abstol = ATOL;

  dx = data->dx = XMAX/((realtype)(MX+1));   /* Set grid coefficients in data */
  data->hdcoef = 1.0/(dx*dx);
  data->hacoef = 0.5/(2.0*dx);

  SetIC(u, dx, local_N, my_base);  /* Initialize u vector */

  /* 
     Call CVodeCreate to create CVODE memory:
     
     ADAMS   specifies the Adams Method
     FUNCTIONAL  specifies functional iteration

     A pointer to CVODE problem memory is returned and stored in cvode_mem.
  */

  cvode_mem = CVodeCreate(ADAMS, FUNCTIONAL);
  if (cvode_mem == NULL) { 
    if (my_pe == 0) printf("CVodeCreate failed.\n"); 
    return(1); 
  }

  flag = CVodeSetFdata(cvode_mem, data);
  if (flag != SUCCESS) { 
    if (my_pe == 0) printf("CVodeSetFdata failed.\n"); 
    return(1); 
  }

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
  if (flag != SUCCESS) { 
    if (my_pe == 0) printf("CVodeMalloc failed.\n"); 
    return(1); 
  }

  if (my_pe == 0) {
    printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", MX);
    printf("\n Number of PEs = %3d \n\n",npes); }

  umax = N_VMaxNorm(u);

  if (my_pe == 0) 
    printf("At t = %4.2f    max.norm(u) =%14.6e \n", T0,umax);

  /* In loop over output points, call CVode, print results, test for error */

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    flag = CVode(cvode_mem, tout, u, &t, NORMAL);
    if (flag != SUCCESS) { 
      if (my_pe == 0) printf("CVode failed, flag=%d.\n", flag); 
      break; 
    }
    umax = N_VMaxNorm(u);
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    if (my_pe == 0) 
      printf("At t = %4.2f  max.norm(u) =%14.6e  nst =%4d \n", t,umax,nst);
  }

  if (my_pe == 0) 
    PrintFinalStats(cvode_mem);     /* Print some final statistics   */

  N_VFree(u);                  /* Free the u vector */
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */
  free(data);                  /* Free user data */
  NV_SpecFree_Parallel(nvSpec);
  MPI_Finalize();

  return(0);
}


/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, realtype dx, integertype my_length,
                  integertype my_base)
{
  int i;
  integertype iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  udata = NV_DATA_P(u);
  my_length = NV_LOCLENGTH_P(u);

  /* Load initial profile into u vector */
  for (i=1; i<=my_length; i++) {
    iglobal = my_base + i;
    x = iglobal*dx;
    udata[i-1] = x*(XMAX - x)*exp(2.0*x);
  }  
}


/* Print some final statistics located in the iopt array */

static void PrintFinalStats(void *cvode_mem)
{
  int nst, nfe, nni, ncfn, netf;
  
  CVodeGetNumSteps(cvode_mem, &nst);
  CVodeGetNumRhsEvals(cvode_mem, &nfe);
  CVodeGetNumErrTestFails(cvode_mem, &netf);
  CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  printf("\nFinal Statistics.. \n\n");
  printf("nst = %-6d  nfe  = %-6d  ", nst, nfe);
  printf("nni = %-6d  ncfn = %-6d  netf = %d\n \n", nni, ncfn, netf);
}


/***************** Function Called by the CVODE Solver ******************/

/* f routine. Compute f(t,u). */

static void f(realtype t, N_Vector u, N_Vector udot, void *f_data)
{
  realtype ui, ult, urt, hordc, horac, hdiff, hadv;
  realtype *udata, *dudata, *z;
  int i;
  int npes, my_pe, my_length, my_pe_m1, my_pe_p1, last_pe, my_last;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  udata = NV_DATA_P(u);
  dudata = NV_DATA_P(udot);

  /* Extract needed problem constants from data */
  data = (UserData) f_data;
  hordc = data->hdcoef;
  horac = data->hacoef;

  /* Extract parameters for parallel computation. */
  comm = data->comm;
  npes = data->npes;           /* Number of processes. */ 
  my_pe = data->my_pe;         /* Current process number. */
  my_length = NV_LOCLENGTH_P(u); /* Number of local elements of u. */ 
  z = data->z;

  /* Compute related parameters. */
  my_pe_m1 = my_pe - 1;
  my_pe_p1 = my_pe + 1;
  last_pe = npes - 1;
  my_last = my_length - 1;

  /* Store local segment of u in the working array z. */
   for (i = 1; i <= my_length; i++)
     z[i] = udata[i - 1];

  /* Pass needed data to processes before and after current process. */
   if (my_pe != 0)
     MPI_Send(&z[1], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
   if (my_pe != last_pe)
     MPI_Send(&z[my_length], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);   

  /* Receive needed data from processes before and after current process. */
   if (my_pe != 0)
     MPI_Recv(&z[0], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
   else z[0] = 0.0;
   if (my_pe != last_pe)
     MPI_Recv(&z[my_length+1], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm,
              &status);   
   else z[my_length + 1] = 0.0;

  /* Loop over all grid points in current process. */
  for (i=1; i<=my_length; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = z[i];
    ult = z[i-1];
    urt = z[i+1];

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - 2.0*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i-1] = hdiff + hadv;
  }
}
