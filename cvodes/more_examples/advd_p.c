#include <stdio.h>
#include <math.h>
#include "llnltyps.h"
#include "cvodes.h"
#include "nvector.h"
#include "mpi.h"

/* Problem Constants */
#define XMAX  2.0          /* domain boundary           */
#define MX    10           /* mesh dimension            */
#define NEQ   MX           /* number of equations       */
#define ATOL  1.e-5        /* scalar absolute tolerance */
#define T0    0.0          /* initial time              */
#define T1    0.5          /* first output time         */
#define DTOUT 0.5          /* output time increment     */
#define NOUT  10           /* number of output times    */

#define NP    2
#define NS    2

#define ZERO  RCONST(0.0)

/* Type : UserData 
   contains problem parameters, grid constants, work array. */

typedef struct {
  real *p;
  real dx;
  integer npes, my_pe;
  MPI_Comm comm;
  real z[100];
} *UserData;


/* Private Helper Functions */

static void SetIC(N_Vector u, real dx, integer my_length, integer my_base);

static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[]);

/* Functions Called by the CVODE Solver */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data);


/***************************** Main Program ******************************/

main(int argc, char *argv[])
{
  real ropt[OPT_SIZE], dx, reltol, abstol, t, tout, umax;
  long int iopt[OPT_SIZE];
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, flag;
  integer local_N, my_pe, npes, nperpe, nrem, my_base;
  machEnvType machEnv;

  real *pbar, rhomax;
  integer is, *plist;
  N_Vector *uS;
  int sensi, sensi_meth, err_con, ifS, sensi_flags[3];

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

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data); /* Allocate data memory */
  data->comm = comm;
  data->npes = npes;
  data->my_pe = my_pe;
  data->p = (real *) malloc(NP * sizeof(real));
  dx = data->dx = XMAX/((real)(MX+1));
  data->p[0] = 1.0;
  data->p[1] = 0.5;

  /* SET machEnv BLOCK */
  machEnv = PVecInitMPI(comm, local_N, NEQ, &argc, &argv);
  if (machEnv == NULL) return(1);

  /* INITIAL STATES */
  u = N_VNew(NEQ, machEnv);          /* Allocate u vector */
  SetIC(u, dx, local_N, my_base);    /* Initialize u vector */

  /* TOLERANCES */
  reltol = 0.0;                /* Set the tolerances */
  abstol = ATOL;

  if(my_pe == 0) {
    printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", MX);
    printf("\n Number of PEs = %3d \n\n",npes);
  }

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(NEQ, f, T0, u, ADAMS, FUNCTIONAL, SS, &reltol,
                          &abstol, data, NULL, FALSE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed.\n"); return(1); }

  /* SENSITIVTY DATA */
  if (my_pe == 0) {
    printf("\nPerform sensitivity analysis? (0:NO , 1:YES): ");scanf("%d",&sensi);
    if(sensi) {
      printf("\nSensitivity method (%d:SIMULTANEOUS , %d:STAGGERED, %d:STAGGERED1): ",
             SIMULTANEOUS,STAGGERED,STAGGERED1);
      scanf("%d",&sensi_meth);
      printf("\nError control (%d:FULL , %d:PARTIAL): ",FULL,PARTIAL);
      scanf("%d",&err_con);
    }
    sensi_flags[0] = sensi;
    sensi_flags[1] = sensi_meth;
    sensi_flags[2] = err_con;
  }
  
  MPI_Bcast(sensi_flags,3,MPI_INT,0,comm);
  sensi      = sensi_flags[0];
  sensi_meth = sensi_flags[1];
  err_con    = sensi_flags[2];

  if(sensi) {
    pbar  = (real *) malloc(NP * sizeof(real));
    pbar[0] = 1.0;
    pbar[1] = 0.5;
    plist = (integer *) malloc(NS * sizeof(integer));
    for(is=0; is<NS; is++)
      plist[is] = is+1; /* sensitivity w.r.t. i-th parameter */

    uS = N_VNew_S(NS,NEQ,machEnv);
    for(is=0;is<NS;is++)
      N_VConst(0.0,uS[is]);

    rhomax = ZERO;
    ifS = ALLSENS;
    if(sensi_meth==STAGGERED1) ifS = ONESENS;
    flag = CVodeSensMalloc(cvode_mem,NS,sensi_meth,data->p,pbar,plist,
                           ifS,NULL,err_con,rhomax,uS,NULL,NULL);
    if (flag != SUCCESS) {
      printf("CVodeSensMalloc failed, flag=%d\n",flag);
      return(1);
    }
  }

  umax = N_VMaxNorm(u);
  if (my_pe == 0)
    printf("At t = %4.2f    max.norm(u) =%14.6e \n", T0, umax);

  if (sensi) {
    for (is=0;is<NS;is++) {
      umax = N_VMaxNorm(uS[is]);
      if (my_pe == 0)
        printf("sensitivity s_%d:  max.norm =%14.6e \n", is, umax);
    }
    if (my_pe == 0) printf("\n");
  }

  /* In loop over output points, call CVode, print results, test for error */

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {

    flag = CVode(cvode_mem, tout, u, &t, NORMAL);

    umax = N_VMaxNorm(u);
    if (my_pe == 0)
      printf("At t = %4.2f    max.norm(u) =%14.6e   nst =%4d \n", t,umax,iopt[NST]);
    
    if (flag != SUCCESS) { 
      printf("CVode failed, flag=%d.\n", flag); 
      break; 
    }

    if (sensi) {
      flag = CVodeSensExtract(cvode_mem, t, uS);
      for (is=0; is<NS; is++) {
        umax = N_VMaxNorm(uS[is]);
        if (my_pe == 0)
          printf("sensitivity s_%d:  max.norm =%14.6e \n", is, umax);
      }
      if (my_pe == 0) printf("\n");
    }

  }

  if (my_pe == 0) 
    PrintFinalStats(sensi,sensi_meth,err_con,iopt);

  N_VFree(u);                  /* Free the u vector */
  if(sensi) N_VFree_S(NS, uS); /* Free the uS vectors */
  CVodeFree(cvode_mem);        /* Free the CVODE problem memory */

  free(data->p);  
  free(data);                  /* Free block of UserData */

  PVecFreeMPI(machEnv);
  MPI_Finalize();

  return(0);
}


/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, real dx, integer my_length, integer my_base)
{
  int i;
  integer iglobal;
  real x;
  real *udata;

  /* Set pointer to data array and get local length of u. */
  udata = N_VDATA(u);
  my_length = N_VLOCLENGTH(u);

  /* Load initial profile into u vector */
  for (i=1; i<=my_length; i++) {
    iglobal = my_base + i;
    x = iglobal*dx;
    udata[i-1] = x*(XMAX - x)*exp(2.0*x);
  }  
}

/* Print some final statistics located in the iopt array */

static void PrintFinalStats(int sensi, int sensi_meth, int err_con, long int iopt[])
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

  printf("========================================================\n");

}

/***************** Function Called by the CVODE Solver ******************/

/* f routine. Compute f(t,u). */

static void f(integer N, real t, N_Vector u, N_Vector udot, void *f_data)
{
  real ui, ult, urt, hordc, horac, hdiff, hadv;
  real *udata, *dudata, *z;
  real dx;
  int i, j;
  int npes, my_pe, my_length, my_pe_m1, my_pe_p1, last_pe, my_last;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  udata = N_VDATA(u);
  dudata = N_VDATA(udot);

  /* Extract needed problem constants from data */
  data  = (UserData) f_data;
  dx    = data->dx; 
  hordc = data->p[0]/(dx*dx);
  horac = data->p[1]/(2.0*dx);

  /* Extract parameters for parallel computation. */
  comm = data->comm;
  npes = data->npes;           /* Number of processes. */ 
  my_pe = data->my_pe;         /* Current process number. */
  my_length = N_VLOCLENGTH(u); /* Number of local elements of u. */ 
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
