/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2004-10-15 20:26:52 $
 * -----------------------------------------------------------------
 * Programmer(s): Lukas Jager and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Parallel Krylov adjoint sensitivity example problem.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mpi.h"
#include "cvodes.h"
#include "cvodea.h"
#include "cvspgmr.h" 
#include "cvbbdpre.h"
#include "nvector_parallel.h"
#include "sundialstypes.h"
#include "sundialsmath.h"

/*
 *------------------------------------------------------------------
 * Constants
 *------------------------------------------------------------------
 */

#ifdef USE3D
#define DIM 3
#else
#define DIM 2
#endif

/* Domain definition */

#define XMIN  0.0
#define XMAX 20.0
#define MX   80
#define NPX  2

#define YMIN  0.0
#define YMAX 20.0
#define MY   80
#define NPY  4

#define ZMIN  0.0
#define ZMAX 16.0
#define MZ   16
#define NPZ  1

/* Parameters for source Gaussians */

#define G1_AMPL   1.0
#define G1_SIGMA  1.7  
#define G1_X      4.0
#define G1_Y      8.0
#define G1_Z      8.0

#define G2_AMPL   0.8
#define G2_SIGMA  3.0  
#define G2_X     16.0
#define G2_Y     12.0
#define G2_Z     12.0

#define G_MIN    1.0e-5

#define DIFF_COEF 1.0
#define V_MAX     1.0
#define L         (YMAX-YMIN)/2.0
#define V_COEFF   V_MAX/L/L

#define ti    0.0
#define tf   10.0

/* Integration tolerances */

#define RTOL    1.0e-8
#define ATOL    1.0e-6

#define RTOL_Q  1.0e-8
#define ATOL_Q  1.0e-6

#define RTOL_B  1.0e-8
#define ATOL_B  1.0e-6

#define RTOL_QB 1.0e-8
#define ATOL_QB 1.0e-6

/* Steps between check points */

#define STEPS 200

/*
 *------------------------------------------------------------------
 * Macros
 *------------------------------------------------------------------
 */

#define FOR_DIM for(dim=0; dim<DIM; dim++)

#ifdef USE3D
#define IJth(y,i) ( y[(i[0])+(l_m[0]*((i[1])+(i[2])*l_m[1]))] )
#define IJth_ext(y,i)     ( y[(i[0]+1)+((l_m[0]+2)*((i[1]+1)+(i[2]+1)*(l_m[1]+2)))] )
#else
#define IJth(y,i) (y[i[0]+(i[1])*l_m[0]])
#define IJth_ext(y,i)     (y[ (i[0]+1) + (i[1]+1) * (l_m[0]+2)])
#endif

/*
 *------------------------------------------------------------------
 * Type definition: ProblemData 
 *------------------------------------------------------------------
 */

typedef struct {
  /* Domain */
  realtype xmin[DIM];  /* "left" boundaries */  
  realtype xmax[DIM];  /* "right" boundaries */
  int m[DIM];          /* number of grid points */
  realtype dx[DIM];    /* grid spacing */
  realtype dOmega;     /* differential volume */

  /* Parallel stuff */
  MPI_Comm comm;       /* MPI communicator */
  int myId;            /* process id */ 
  int npes;            /* total number of processes */
  int num_procs[DIM];  /* number of processes in each direction */
  int nbr_left[DIM];   /* MPI ID of "left" neighbor */
  int nbr_right[DIM];  /* MPI ID of "right" neighbor */
  int m_start[DIM];    /* "left" index in the global domain */
  int l_m[DIM];        /* number of local grid points */ 
  realtype *y_ext;     /* extended data array */
  realtype *buf_send;  /* Send buffer */
  realtype *buf_recv;  /* Receive buffer */
  int buf_size;        /* Buffer size */

  /* Source */
  N_Vector p;          /* Source parameters */ 

} *ProblemData;

/*
 *------------------------------------------------------------------
 * Interface functions to CVODES
 *------------------------------------------------------------------
 */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static void f_local(long int Nlocal, realtype t, N_Vector y, 
                    N_Vector ydot, void *f_data);
static void ucomm(long int Nlocal, realtype t, N_Vector y, void *f_data);

static void fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data);


static void fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, 
		      void *f_dataB);
static void fB_local(long int NlocalB, realtype t, 
                     N_Vector y, N_Vector yB, N_Vector yBdot, 
                     void *f_dataB);
static void ucommB(long int NlocalB, realtype t, 
		   N_Vector y, N_Vector yB,
		   void *f_dataB);

static void fQB(realtype t, N_Vector y, N_Vector yB, 
                N_Vector qBdot, void *fQ_dataB);

/*
 *------------------------------------------------------------------
 * Private functions
 *------------------------------------------------------------------
 */

static void SetData(ProblemData d, MPI_Comm comm, int npes, int myId,
                    long int *neq, long int *l_neq);
static void SetSource(ProblemData d);
static void f_comm( long int Nlocal, realtype t, N_Vector y, void *f_data);
static void Load_yext(realtype *src, ProblemData d);
static void PrintFinalStats(void *cvode_mem);
static void OutputGradient(int myId, N_Vector qB, ProblemData d);

/*
 *------------------------------------------------------------------
 * Main program
 *------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  ProblemData d;

  MPI_Comm comm;
  int npes, npes_needed;
  int myId;
 
  long int neq, l_neq;

  void *cvode_mem;
  N_Vector y, q;
  realtype abstol, reltol, abstolQ, reltolQ;
  void *bbdp_data;
  int mudq, mldq, mukeep, mlkeep;

  void *cvadj_mem;
  void *cvode_memB;
  N_Vector yB, qB;
  realtype abstolB, reltolB, abstolQB, reltolQB;
  int mudqB, mldqB, mukeepB, mlkeepB;

  realtype tret, *qdata, G;

  int ncheckpnt, flag;

  /* Initialize MPI and set Ids */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myId);

  /* Check number of processes */
  npes_needed = NPX * NPY;
#ifdef USE3D
  npes_needed *= NPZ;
#endif
  MPI_Comm_size(comm, &npes);
  if (npes_needed != npes) {
    if (myId == 0)
      fprintf(stderr,"I need %d processes but I only got %d\n",
              npes_needed, npes);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  /* Allocate and set problem data structure */
  d = (ProblemData) malloc(sizeof *d);
  SetData(d, comm, npes, myId, &neq, &l_neq);
  
  if (myId == 0) {
    printf("\nParallel Krylov adjoint sensitivity analysis example\n");
    printf("%1dD Advection diffusion PDE with homogeneous Neumann B.C.\n",DIM);
    printf("Computes gradient of G = int_t_Omega ( c_i^2 ) dt dOmega\n");
    printf("with respect to the source values at each grid point.\n\n");

    printf("Domain:\n");
    printf("   %f < x < %f   mx = %d  npe_x = %d \n",XMIN,XMAX,MX,NPX);
    printf("   %f < y < %f   my = %d  npe_y = %d \n",YMIN,YMAX,MY,NPY);
#ifdef USE3D
    printf("   %f < z < %f   mz = %d  npe_z = %d \n",ZMIN,ZMAX,MZ,NPZ);
#endif

    printf("\n");
  }

  /*-------------------------- 
    Forward integration phase
    --------------------------*/

  /* Allocate space for y and set it with the I.C. */
  y = N_VNew_Parallel(comm, l_neq, neq);
  N_VConst(0.0, y);
  
  /* Allocate and initialize qB (local contributin to cost) */
  q = N_VNew_Parallel(comm, 1, npes); 
  N_VConst(0.0, q);

  /* Create CVODES object, attach user data, and allocate space */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeSetFdata(cvode_mem, d);
  abstol = ATOL;  
  reltol = RTOL; 
  flag = CVodeMalloc(cvode_mem, f, ti, y, CV_SS, &reltol, &abstol);
  
  /* Attach preconditioner and linear solver modules */
  mudq = mldq = d->l_m[0]+1;
  mukeep = mlkeep = 2;  
  bbdp_data = (void *) CVBBDPrecAlloc(cvode_mem, l_neq, mudq, mldq, 
                                      mukeep, mlkeep, 0.0, f_local, ucomm);
  flag = CVBBDSpgmr(cvode_mem, PREC_LEFT, 0, bbdp_data);

  /* Initialize quadrature calculations */
  abstolQ = ATOL_Q;
  reltolQ = RTOL_Q;
  flag = CVodeSetQuadFdata(cvode_mem, d);
  flag = CVodeSetQuadErrCon(cvode_mem, TRUE);
  flag = CVodeSetQuadTolerances(cvode_mem, CV_SS, &reltolQ, &abstolQ); 
  flag = CVodeQuadMalloc(cvode_mem, fQ, q);

  /* Allocate space for the adjoint calculation */
  cvadj_mem = CVadjMalloc(cvode_mem, STEPS);

  /* Integrate forward in time while storing check points */
  if (myId == 0) printf("Begin forward integration... ");
  flag = CVodeF(cvadj_mem, tf, y, &tret, CV_NORMAL, &ncheckpnt);
  if (myId == 0) printf("done. ");

   /* Extract quadratures */
  flag = CVodeGetQuad(cvode_mem, tf, q);
  qdata = NV_DATA_P(q);
  MPI_Allreduce(&qdata[0], &G, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
  if (myId == 0) printf("  G = %e\n",G);

  /* Print statistics for forward run */
  if (myId == 0) PrintFinalStats(cvode_mem);

  /*-------------------------- 
    Backward integration phase
    --------------------------*/
 
  /* Allocate and initialize yB */
  yB = N_VNew_Parallel(comm, l_neq, neq); 
  N_VConst(0.0, yB);

  /* Allocate and initialize qB (gradient) */
  qB = N_VNew_Parallel(comm, l_neq, neq); 
  N_VConst(0.0, qB);

  /* Create and allocate backward CVODE memory */
  flag = CVodeCreateB(cvadj_mem, CV_BDF, CV_NEWTON);
  flag = CVodeSetFdataB(cvadj_mem, d);
  flag = CVodeSetMaxNumStepsB(cvadj_mem, 10000);
  abstolB = ATOL_B;  
  reltolB = RTOL_B; 
  flag = CVodeMallocB(cvadj_mem, fB, tf, yB, CV_SS, &reltolB, &abstolB);

  /* Attach preconditioner and linear solver modules */
  mudqB = mldqB = d->l_m[0]+1;
  mukeepB = mlkeepB = 2;  
  flag = CVBBDPrecAllocB(cvadj_mem, l_neq, mudqB, mldqB, 
                         mukeepB, mlkeepB, 0.0, fB_local, ucommB);
  flag = CVBBDSpgmrB(cvadj_mem, PREC_LEFT, 0); 

  /* Initialize quadrature calculations */
  abstolQB = ATOL_QB;
  reltolQB = RTOL_QB;
  flag = CVodeSetQuadFdataB(cvadj_mem, d);
  flag = CVodeSetQuadErrConB(cvadj_mem, TRUE);
  flag = CVodeSetQuadTolerancesB(cvadj_mem, CV_SS, &reltolQB, &abstolQB); 
  flag = CVodeQuadMallocB(cvadj_mem, fQB, qB);

  /* Integrate backwards */
  if (myId == 0) printf("Begin backward integration... ");
  flag = CVodeB(cvadj_mem, ti, yB, &tret, CV_NORMAL);
  if (myId == 0) printf("done.\n");
  
  /* Print statistics for backward run */
  if (myId == 0) {
    cvode_memB = CVadjGetCVodeBmem(cvadj_mem);
    PrintFinalStats(cvode_memB);
  }

   /* Extract quadratures */
  flag = CVodeGetQuadB(cvadj_mem, qB);

  /* Process 0 collects the gradient components and prints them */
  OutputGradient(myId, qB, d);
  if (myId == 0) printf("Wrote matlab file 'grad.m'.\n");

  /* Free memory */
  N_VDestroy_Parallel(y);
  N_VDestroy_Parallel(q);
  N_VDestroy_Parallel(qB);
  N_VDestroy_Parallel(yB);

  CVBBDPrecFree(bbdp_data);
  CVadjFree(cvadj_mem);
  CVodeFree(cvode_mem);

  MPI_Finalize();

}

/*
 *------------------------------------------------------------------
 * SetData:
 * Allocate space for the ProblemData structure.
 * Set fields in the ProblemData structure.
 * Return local and global problem dimensions.
 *
 * SetSource:
 * Instantiates the source parameters for a combination of two
 * Gaussian sources.
 *------------------------------------------------------------------
 */

static void SetData(ProblemData d, MPI_Comm comm, int npes, int myId,
                    long int *neq, long int *l_neq)
{
  int n[DIM], nd[DIM];
  int dim, size;

  /* Set MPI communicator, id,  and total number of processes */
  d->comm = comm;
  d->myId = myId;
  d->npes = npes;

  /* Set domain boundaries */

  d->xmin[0] = XMIN;
  d->xmax[0] = XMAX;
  d->m[0]    = MX;

  d->xmin[1] = YMIN;
  d->xmax[1] = YMAX;
  d->m[1]    = MY;
#ifdef USE3D
  d->xmin[2] = ZMIN;
  d->xmax[2] = ZMAX;
  d->m[2]    = MZ;
#endif

  /* Calculate grid spacing and differential volume */
  d->dOmega = 1.0;
  FOR_DIM {
    d->dx[dim] = (d->xmax[dim] - d->xmin[dim]) / d->m[dim];
    d->m[dim] +=1;
    d->dOmega *= d->dx[dim];
  }

  /* Set partitioning */

  d->num_procs[0] = NPX;
  n[0] = NPX; 
  nd[0] = d->m[0] / NPX;

  d->num_procs[1] = NPY;
  n[1] = NPY; 
  nd[1] = d->m[1] / NPY;
#ifdef USE3D
  d->num_procs[2] = NPZ;
  n[2] = NPZ; 
  nd[2] = d->m[2] / NPZ;
#endif
  
  /* Compute the neighbors */

  d->nbr_left[0]  = (myId%n[0]) == 0                ? myId : myId-1;
  d->nbr_right[0] = (myId%n[0]) == n[0]-1           ? myId : myId+1;

  d->nbr_left[1]  = (myId/n[0])%n[1] == 0           ? myId : myId-n[0];
  d->nbr_right[1] = (myId/n[0])%n[1] == n[1]-1      ? myId : myId+n[0];
#ifdef USE3D
  d->nbr_left[2]  = (myId/n[0]/n[1])%n[2] == 0      ? myId : myId-n[0]*n[1];
  d->nbr_right[2] = (myId/n[0]/n[1])%n[2] == n[2]-1 ? myId : myId+n[0]*n[1];
#endif
 
  /* Compute the local subdomains 
     m_start: left border in global index space 
     l_m:     length of the subdomain */
  
  d->m_start[0] = (myId%n[0])*nd[0];
  d->l_m[0]     = d->nbr_right[0] == myId ? d->m[0] - d->m_start[0] : nd[0];

  d->m_start[1] = ((myId/n[0])%n[1])*nd[1];
  d->l_m[1]     = d->nbr_right[1] == myId ? d->m[1] - d->m_start[1] : nd[1];

#ifdef USE3D
  d->m_start[2] = (myId/n[0]/n[1])*nd[2];
  d->l_m[2]     = d->nbr_right[2] == myId ? d->m[2] - d->m_start[2] : nd[2];
#endif

  /* Allocate memory for the y_ext array 
     (local solution + ghost cells) */
  size = 1;
  FOR_DIM size *= d->l_m[dim]+2;
  d->y_ext = (realtype *) malloc( size*sizeof(realtype));

  /* Initialize Buffer field.
     Size of buffer is checked when needed */
  d->buf_send = NULL;
  d->buf_recv = NULL;
  d->buf_size = 0;   


  /* Allocate space for the source parameters */
  *neq = 1; *l_neq = 1;
  FOR_DIM {*neq *= d->m[dim];  *l_neq *= d->l_m[dim];}
  d->p = N_VNew_Parallel(comm, *l_neq, *neq);

  /* Initialize the parameters for a source with Gaussian profile */
  SetSource(d);

}

static void SetSource(ProblemData d)
{
  int *l_m, *m_start;
  realtype *xmin, *xmax, *dx;
  realtype x[DIM], g, *pdata;
  int i[DIM];

  l_m  = d->l_m;
  m_start = d->m_start;
  xmin = d->xmin;
  xmax = d->xmax;
  dx = d->dx;


  pdata = NV_DATA_P(d->p);

  for(i[0]=0; i[0]<l_m[0]; i[0]++) {
    x[0] = xmin[0] + (m_start[0]+i[0]) * dx[0];
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      x[1] = xmin[1] + (m_start[1]+i[1]) * dx[1];
#ifdef USE3D
      for(i[2]=0; i[2]<l_m[2]; i[2]++) {
        x[2] = xmin[2] + (m_start[2]+i[2]) * dx[2];
        
        g = G1_AMPL 
          * exp( -SQR(G1_X-x[0])/SQR(G1_SIGMA) ) 
          * exp( -SQR(G1_Y-x[1])/SQR(G1_SIGMA) )
          * exp( -SQR(G1_Z-x[2])/SQR(G1_SIGMA) ); 
        
        g += G2_AMPL 
          * exp( -SQR(G2_X-x[0])/SQR(G2_SIGMA) ) 
          * exp( -SQR(G2_Y-x[1])/SQR(G2_SIGMA) )
          * exp( -SQR(G2_Z-x[2])/SQR(G2_SIGMA) ); 
        
        if( g < G_MIN ) g = 0.0;

        IJth(pdata, i) = g;
      }
#else
      g = G1_AMPL 
        * exp( -SQR(G1_X-x[0])/SQR(G1_SIGMA) ) 
        * exp( -SQR(G1_Y-x[1])/SQR(G1_SIGMA) ); 

      g += G2_AMPL 
        * exp( -SQR(G2_X-x[0])/SQR(G2_SIGMA) ) 
        * exp( -SQR(G2_Y-x[1])/SQR(G2_SIGMA) ); 
      
      if( g < G_MIN ) g = 0.0;

      IJth(pdata, i) = g;
#endif 
    }
  }
}

/*
 *------------------------------------------------------------------
 * f_comm: 
 * Function for inter-process communication
 * Used both for the forward and backward phase.
 *------------------------------------------------------------------
 */

static void f_comm(long int N_local, realtype t, N_Vector y, void *f_data)
{
  int id, n[DIM], proc_cond[DIM], nbr[DIM][2];
  ProblemData d;
  realtype *y_ext, *y_data;
  int l_m[DIM], dim;
  int c, i[DIM], l[DIM-1];
  realtype *buf_send, *buf_recv;
  MPI_Status stat;
  MPI_Comm comm;
  int dir, size = 1, small = INT_MAX;

  d  = (ProblemData) f_data;
  comm = d->comm;
  id = d->myId;
  
  /* extract data from domain*/
  FOR_DIM {
    n[dim] = d->num_procs[dim];
    l_m[dim] = d->l_m[dim];
  }
  y_ext = d->y_ext;
  y_data = NV_DATA_P(y);
  
  /* Calculate required buffer size */
  FOR_DIM {
    size *= l_m[dim];
    if( l_m[dim] < small) small = l_m[dim];
  }
  size /= small;
  
  /* Adjust buffer size if necessary */
  if( d->buf_size < size ) {
    d->buf_send = (realtype*) realloc( d->buf_send, size * sizeof(realtype));
    d->buf_recv = (realtype*) realloc( d->buf_recv, size * sizeof(realtype));
    d->buf_size = size;
  }

  buf_send = d->buf_send;
  buf_recv = d->buf_recv;
  
  /* Compute the communication pattern; who sends first? */
  /* if proc_cond==1 , process sends first in this dimension */
  proc_cond[0] = (id%n[0])%2;
  proc_cond[1] = ((id/n[0])%n[1])%2;
#ifdef USE3D
  proc_cond[2] = (id/n[0]/n[1])%2;
#endif

  /* Compute the actual communication pattern */
  /* nbr[dim][0] is first proc to communicate with in dimension dim */
  /* nbr[dim][1] the second one */
  FOR_DIM {
    nbr[dim][proc_cond[dim]]  = d->nbr_left[dim];
    nbr[dim][!proc_cond[dim]] = d->nbr_right[dim];
  }
  
  /* Communication: loop over dimension and direction (left/right) */
  FOR_DIM {

    for (dir=0; dir<=1; dir++) {

      /* If subdomain at boundary, no communication in this direction */

      if (id != nbr[dim][dir]) {
        c=0;
        /* Compute the index of the boundary (right or left) */
        i[dim] = (dir ^ proc_cond[dim]) ? (l_m[dim]-1) : 0;
        /* Loop over all other dimensions and copy data into buf_send */
        l[0]=(dim+1)%DIM;
#ifdef USE3D
        l[1]=(dim+2)%DIM;
        for(i[l[1]]=0; i[l[1]]<l_m[l[1]]; i[l[1]]++) 
#endif
          for(i[l[0]]=0; i[l[0]]<l_m[l[0]]; i[l[0]]++) 
            buf_send[c++] = IJth(y_data, i);
	  
        if ( proc_cond[dim] ) {
          /* Send buf_send and receive into buf_recv */
          MPI_Send(buf_send, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm);
          MPI_Recv(buf_recv, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm, &stat);
        } else {
          /* Receive into buf_recv and send buf_send*/
          MPI_Recv(buf_recv, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm, &stat);
          MPI_Send(buf_send, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm);
        }

        c=0;

        /* Compute the index of the boundary (right or left) in y_ext */
        i[dim] = (dir ^ proc_cond[dim]) ? l_m[dim] : -1;

        /* Loop over all other dimensions and copy data into y_ext */
#ifdef USE3D
        for(i[l[1]]=0; i[l[1]]<l_m[l[1]]; i[l[1]]++)
#endif
          for(i[l[0]]=0; i[l[0]]<l_m[l[0]]; i[l[0]]++)
            IJth_ext(y_ext, i) = buf_recv[c++];
      }
    } /* end loop over direction */
  } /* end loop over dimension */ 
}

/*
 *------------------------------------------------------------------
 * f and f_local:
 * Forward phase ODE right-hand side
 *------------------------------------------------------------------
 */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  ProblemData d;
  int l_neq=1;
  int dim;

  d = (ProblemData) f_data;
  FOR_DIM l_neq *= d->l_m[dim];
  
  /* Do all inter-processor communication */
  f_comm(l_neq, t, y, f_data);

  /* Compute right-hand side locally */
  f_local(l_neq, t, y, ydot, f_data);
}


static void f_local(long int Nlocal, realtype t, N_Vector y, 
                    N_Vector ydot, void *f_data)
{
  realtype *yextdata, *dydata, *ydata, *pdata;
  realtype dx[DIM], c, v[DIM], cl[DIM], cr[DIM];
  realtype adv[DIM], diff[DIM];
  realtype xmin[DIM], xmax[DIM], x[DIM], x1;
  int i[DIM], l_m[DIM], m_start[DIM];
  ProblemData d;
  int dim;

  d = (ProblemData) f_data;

  FOR_DIM {
    xmin[dim]    = d->xmin[dim];
    xmax[dim]    = d->xmax[dim];
    l_m[dim]     = d->l_m[dim];
    m_start[dim] = d->m_start[dim];
    dx[dim]      = d->dx[dim];
  } 

  yextdata = d->y_ext;
  dydata   = NV_DATA_P(ydot);
  ydata    = NV_DATA_P(y); 
  pdata    = NV_DATA_P(d->p);

  /* Copy local segment of y to y_ext, boundary */
  Load_yext(ydata, d);

  /* Velocity components in x1 and x2 directions (Poiseuille profile) */
  v[1] = 0.0;
#ifdef USE3D
  v[2] = 0.0;
#endif

  /* Local domain is [xmin+(m_start+1)*dx, xmin+(m_start+1+l_m-1)*dx] */
#ifdef USE3D
  for(i[2]=0; i[2]<l_m[2]; i[2]++) {

    x[2] = xmin[2] + (m_start[2]+i[2])*dx[2];
#endif    
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {

      x[1] = xmin[1] + (m_start[1]+i[1])*dx[1];

      /* Velocity component in x0 direction (Poiseuille profile) */
      x1 = x[1] - xmin[1] - L;
      v[0] = V_COEFF * (L + x1) * (L - x1);

      for(i[0]=0; i[0]<l_m[0]; i[0]++) {

        x[0] = xmin[0] + (m_start[0]+i[0])*dx[0];

        c  = IJth_ext(yextdata, i);	       

        /* Source term*/
        IJth(dydata, i) = IJth(pdata, i);

        FOR_DIM {
          i[dim]+=1;
          cr[dim] = IJth_ext(yextdata, i);
          i[dim]-=2;
          cl[dim] = IJth_ext(yextdata, i);
          i[dim]+=1;
          adv[dim]  = v[dim] * (cr[dim]-cl[dim]) / (2.0*dx[dim]);
          diff[dim] = DIFF_COEF * (cr[dim]-2.*c+cl[dim]) / SQR(dx[dim]);

          IJth(dydata, i) += (diff[dim] - adv[dim]);
        } 
      }
    }
#ifdef USE3D
  }
#endif
}

/*
 *------------------------------------------------------------------
 * Communication function for CVBBDPRE for forwardintegration.
 * Nothing is done here, as all necessary communication was already
 * done in f_comm.
 *------------------------------------------------------------------
 */

static void ucomm(long int Nlocal, realtype t, N_Vector y, 
                  void *f_data){}

/*
 *------------------------------------------------------------------
 * fQ:
 * Right-hand side of quadrature equations on forward integration.
 * The only quadrature on this phase computes the local contribution
 * to the function G.
 *------------------------------------------------------------------
 */

static void fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data)
{
  ProblemData d;
  realtype *dqdata;

  d = (ProblemData) fQ_data;

  dqdata = NV_DATA_P(qdot);

  dqdata[0] = N_VDotProd_Parallel(y,y);
  dqdata[0] *= 0.5 * (d->dOmega);
}

/*
 *------------------------------------------------------------------
 * fB and fB_local:
 * Backward phase ODE right-hand side (the discretized adjoint PDE)
 *------------------------------------------------------------------
 */

static void fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, 
               void *f_dataB)
{
  ProblemData d;
  int l_neq=1;
  int dim;

  d = (ProblemData) f_dataB;
  FOR_DIM l_neq *= d->l_m[dim];
  
  /* Do all inter-processor communication */
  f_comm(l_neq, t, yB, f_dataB);

  /* Compute right-hand side locally */
  fB_local(l_neq, t, y, yB, yBdot, f_dataB);
}

static void fB_local(long int NlocalB, realtype t, 
                     N_Vector y, N_Vector yB, N_Vector gB, 
                     void *f_dataB)
{
  realtype *yextdata, *dydata, *ydata, *old_y;
  realtype dx[DIM], c, v[DIM], cl[DIM], cr[DIM];
  realtype adv[DIM], diff[DIM];
  realtype xmin[DIM], xmax[DIM], x[DIM], x1;
  int i[DIM], l_m[DIM], m_start[DIM], nbr_left[DIM], nbr_right[DIM], id;
  ProblemData d;
  int dim;
  
  d = (ProblemData) f_dataB;

  /* Extract stuff from data */
  id = d->myId;
  FOR_DIM {
    xmin[dim]      = d->xmin[dim];
    xmax[dim]      = d->xmax[dim];
    l_m[dim]       = d->l_m[dim];
    m_start[dim]   = d->m_start[dim];
    dx[dim]        = d->dx[dim];
    nbr_left[dim]  = d->nbr_left[dim];
    nbr_right[dim] = d->nbr_right[dim];
  }
 
  yextdata = d->y_ext;
  dydata   = NV_DATA_P(gB);
  ydata    = NV_DATA_P(yB); 
  old_y    = NV_DATA_P(y);

  /* Copy local segment of y to y_ext, boundary */
  Load_yext(ydata, d);

  /* Velocity components in x1 and x2 directions (Poiseuille profile) */
  v[1] = 0.0;
#ifdef USE3D
  v[2] = 0.0;
#endif
 
  /* local domain is [xmin+(m_start)*dx, xmin+(m_start+l_m-1)*dx] */
#ifdef USE3D
  for(i[2]=0; i[2]<l_m[2]; i[2]++) {

    x[2] = xmin[2] + (m_start[2]+i[2])*dx[2];
#endif
    
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      
      x[1] = xmin[1] + (m_start[1]+i[1])*dx[1];
	  
      /* Velocity component in x0 direction (Poiseuille profile) */
      x1 = x[1] - xmin[1] - L;
      v[0] = V_COEFF * (L + x1) * (L - x1);

      for(i[0]=0; i[0]<l_m[0]; i[0]++) {

        x[0] = xmin[0] + (m_start[0]+i[0])*dx[0];
        
        c  = IJth_ext(yextdata, i);	       
        
        /* Source term for adjoint PDE */
        IJth(dydata, i) = -IJth(old_y, i);
        
        FOR_DIM {
          
          i[dim]+=1;
          cr[dim] = IJth_ext(yextdata, i);
          i[dim]-=2;
          cl[dim] = IJth_ext(yextdata, i);
          i[dim]+=1;

          /* Boundary conditions for the adjoint variables */
          if( i[dim]==l_m[dim]-1 && nbr_right[dim]==id)
            cr[dim] = cl[dim]-(2.0*dx[dim]*v[dim]/DIFF_COEF)*c;
          else if( i[dim]==0 && nbr_left[dim]==id )
            cl[dim] = cr[dim]+(2.0*dx[dim]*v[dim]/DIFF_COEF)*c;
		  
          adv[dim]  = v[dim] * (cr[dim]-cl[dim]) / (2.0*dx[dim]);
          diff[dim] = DIFF_COEF * (cr[dim]-2.*c+cl[dim]) / SQR(dx[dim]);
          
          IJth(dydata, i) -= (diff[dim] + adv[dim]);
        } 
      }
    }
#ifdef USE3D
  }
#endif
}

/*
 *------------------------------------------------------------------
 * Communication function for CVBBDPRE for backward integration.
 * Nothing is done here, as all necessary communication was already
 * done in f_comm.
 *------------------------------------------------------------------
 */

static void ucommB(long int NlocalB, realtype t, 
		   N_Vector y, N_Vector yB,
		   void *f_dataB){}
  
/*
 *------------------------------------------------------------------
 * fQB:
 * Right-hand side of quadrature equations on backward integration
 * The i-th component of the gradient is nothing but int_t yB_i dt
 *------------------------------------------------------------------
 */

static void fQB(realtype t, N_Vector y, 
		N_Vector yB, N_Vector qBdot, 
		void *fQ_dataB)
{
  ProblemData d;

  d = (ProblemData) fQ_dataB;

  N_VScale_Parallel(-(d->dOmega), yB, qBdot);
}

/*
 *------------------------------------------------------------------
 * Load_yext: 
 * copies data from y into y_ext, the array with boundaries
 *------------------------------------------------------------------
 */

static void Load_yext(realtype *src, ProblemData d)
{
  int i[DIM], j[DIM];
  int id, l_m[DIM];
  int dim, l1;
#ifdef USE3D
  int  l2;
#endif

  id = d->myId;
  FOR_DIM l_m[dim] = d->l_m[dim];
     
  /* copy local segment */
#ifdef USE3D
  for  (i[2]=0; i[2]<l_m[2]; i[2]++)
#endif
    for(i[1]=0; i[1]<l_m[1]; i[1]++)
      for(i[0]=0; i[0]<l_m[0]; i[0]++)
	IJth_ext(d->y_ext, i) = IJth(src, i);

  /* copy borders if subdomain at boundary */
  FOR_DIM {
    /* subdomain at boundary ? */
    if( id == d->nbr_left[dim]) {
      /* left boundary*/
      i[dim] = 1;
      j[dim] = -1;
      l1=(dim+1)%DIM;
#ifdef USE3D
      l2=(dim+2)%DIM;
      for(i[l2]=0; i[l2]<l_m[l2]; i[l2]++) {
        j[l2]=i[l2];
#endif
        for(i[l1]=0; i[l1]<l_m[l1]; i[l1]++) {
          j[l1]=i[l1]; 
          IJth_ext(d->y_ext,j) = IJth(src, i);
        }
#ifdef USE3D
      }
#endif
    }
    /* subdomain at boundary ? */
    if( id == d->nbr_right[dim]) {
      /* right boundary*/
      i[dim] = l_m[dim]-2;
      j[dim] = l_m[dim];
      l1=(dim+1)%DIM;
#ifdef USE3D
      l2=(dim+2)%DIM;
      for(i[l2]=0; i[l2]<l_m[l2]; i[l2]++) {
        j[l2]=i[l2];
#endif
        for(i[l1]=0; i[l1]<l_m[l1]; i[l1]++) {
          j[l1]=i[l1]; 
          IJth_ext(d->y_ext,j) = IJth(src, i);
        }
#ifdef USE3D
      }
#endif
    }
  }
}

/*
 *------------------------------------------------------------------
 * PrintFinalStats:
 * Print final statistics contained in cvode_mem
 *------------------------------------------------------------------
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int lenrw, leniw ;
  long int lenrwSPGMR, leniwSPGMR;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeSPGMR;
  int flag;

  flag = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  flag = CVSpgmrGetWorkSpace(cvode_mem, &lenrwSPGMR, &leniwSPGMR);
  flag = CVSpgmrGetNumLinIters(cvode_mem, &nli);
  flag = CVSpgmrGetNumPrecEvals(cvode_mem, &npe);
  flag = CVSpgmrGetNumPrecSolves(cvode_mem, &nps);
  flag = CVSpgmrGetNumConvFails(cvode_mem, &ncfl);
  flag = CVSpgmrGetNumRhsEvals(cvode_mem, &nfeSPGMR);

  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %6ld     leniw = %6ld\n", lenrw, leniw);
  printf("llrw    = %6ld     lliw  = %6ld\n", lenrwSPGMR, leniwSPGMR);
  printf("nst     = %6ld\n"                  , nst);
  printf("nfe     = %6ld     nfel  = %6ld\n"  , nfe, nfeSPGMR);
  printf("nni     = %6ld     nli   = %6ld\n"  , nni, nli);
  printf("nsetups = %6ld     netf  = %6ld\n"  , nsetups, netf);
  printf("npe     = %6ld     nps   = %6ld\n"  , npe, nps);
  printf("ncfn    = %6ld     ncfl  = %6ld\n\n", ncfn, ncfl); 
}

/*
 *------------------------------------------------------------------
 * OutputGradient:
 * Generate matlab m files for visualization
 * One file gradXXXX.m from each process + a driver grad.m
 *------------------------------------------------------------------
 */

static void OutputGradient(int myId, N_Vector qB, ProblemData d)
{
  FILE *fid;
  char filename[20];
  int *l_m, *m_start, i[DIM],ip;
  realtype *xmin, *xmax, *dx;
  realtype x[DIM], *pdata, p, *qBdata, g;

  sprintf(filename,"grad%03d.m",myId);
  fid = fopen(filename,"w");

  l_m  = d->l_m;
  m_start = d->m_start;
  xmin = d->xmin;
  xmax = d->xmax;
  dx = d->dx;

  qBdata = NV_DATA_P(qB);
  pdata  = NV_DATA_P(d->p);

  for(i[0]=0; i[0]<l_m[0]; i[0]++) {
    x[0] = xmin[0] + (m_start[0]+i[0]) * dx[0];
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      x[1] = xmin[1] + (m_start[1]+i[1]) * dx[1];
#ifdef USE3D
      for(i[2]=0; i[2]<l_m[2]; i[2]++) {
        x[2] = xmin[2] + (m_start[2]+i[2]) * dx[2];
        
        g = IJth(qBdata, i);
      }
#else
      g = IJth(qBdata, i);;
      p = IJth(pdata, i);;
      fprintf(fid,"x%d(%d,1) = %e; \n",  myId, i[0]+1,         x[0]);
      fprintf(fid,"y%d(%d,1) = %e; \n",  myId, i[1]+1,         x[1]);
      fprintf(fid,"p%d(%d,%d) = %e; \n", myId, i[1]+1, i[0]+1, p);
      fprintf(fid,"g%d(%d,%d) = %e; \n", myId, i[1]+1, i[0]+1, g);
      
#endif 
    }
  }
  fclose(fid);

  if (myId == 0) {

    fid = fopen("grad.m","w");

    fprintf(fid,"clear;\nfigure;\n");
    fprintf(fid,"trans = 0.7;\n");
    fprintf(fid,"ecol  = 'none';\n");

    for (ip=0; ip<d->npes; ip++) {

      fprintf(fid,"\ngrad%03d;\n",ip);

      fprintf(fid,"\nsubplot(1,2,1)\n");
      fprintf(fid,"s=surf(x%d,y%d,g%d);\n",ip,ip,ip);
      fprintf(fid,"set(s,'FaceAlpha',trans);\n");
      fprintf(fid,"set(s,'EdgeColor',ecol);\n");
      fprintf(fid,"hold on\n");
      
      fprintf(fid,"\nsubplot(1,2,2)\n");
      fprintf(fid,"s=surf(x%d,y%d,p%d);\n",ip,ip,ip);
      fprintf(fid,"set(s,'CData',g%d);\n",ip);
      fprintf(fid,"set(s,'FaceAlpha',trans);\n");
      fprintf(fid,"set(s,'EdgeColor',ecol);\n");
      fprintf(fid,"hold on\n");
    }

    fprintf(fid,"\nsubplot(1,2,1)\n");
    fprintf(fid,"axis tight\n");
    fprintf(fid,"box on\n");

    fprintf(fid,"\nsubplot(1,2,2)\n");
    fprintf(fid,"axis tight\n");
    fprintf(fid,"box on\n");

    fclose(fid);

  }
    
}

