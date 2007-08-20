/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2007-08-20 20:56:23 $
 * -----------------------------------------------------------------
 * Programmer(s): Lukas Jager and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * 3D parallel Krylov adjoint sensitivity example problem.
 * Use for scalability tests.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_spgmr.h> 
#include <cvodes/cvodes_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <mpi.h>

/*
 *------------------------------------------------------------------
 * Constants
 *------------------------------------------------------------------
 */


/* Domain definition */

#define XMIN RCONST(0.0)
#define XMAX RCONST(60.0)
#define MX   60    /* no. of divisions in x dir. */
#define NPX  2     /* no. of procs. in x dir.    */

#define YMIN RCONST(0.0)
#define YMAX RCONST(20.0)
#define MY   20    /* no. of divisions in y dir. */
#define NPY  2     /* no. of procs. in y dir.    */

#define ZMIN RCONST(0.0)
#define ZMAX RCONST(20.0)
#define MZ   20    /* no. of divisions in z dir. */
#define NPZ  1     /* no. of procs. in z dir.    */

/* Parameters for source Gaussians */

#define G1_AMPL   RCONST(1.0)
#define G1_SIGMA  RCONST(1.7) 
#define G1_X      RCONST(4.0)
#define G1_Y      RCONST(8.0)
#define G1_Z      RCONST(8.0)

#define G2_AMPL   RCONST(0.8)
#define G2_SIGMA  RCONST(3.0)
#define G2_X      RCONST(16.0)
#define G2_Y      RCONST(12.0)
#define G2_Z      RCONST(12.0)

#define G_MIN     RCONST(1.0e-5)

/* Diffusion coeff., max. velocity, domain width in y dir. */

#define DIFF_COEF RCONST(1.0)
#define V_MAX     RCONST(1.0)
#define L         (YMAX-YMIN)/RCONST(2.0)
#define V_COEFF   V_MAX/L/L

/* Initial and final times */

#define ti    RCONST(0.0)
#define tf    RCONST(10.0)

/* Integration tolerances */

#define RTOL    RCONST(1.0e-8) /* states */
#define ATOL    RCONST(1.0e-6)

#define RTOL_Q  RCONST(1.0e-8) /* forward quadrature */
#define ATOL_Q  RCONST(1.0e-6)

#define RTOL_B  RCONST(1.0e-8) /* adjoint variables */
#define ATOL_B  RCONST(1.0e-6)

#define RTOL_QB RCONST(1.0e-8) /* backward quadratures */
#define ATOL_QB RCONST(1.0e-6)

/* Steps between check points */

#define STEPS 200

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/*
 *------------------------------------------------------------------
 * Macros
 *------------------------------------------------------------------
 */

#define DIM 3
#define FOR_DIM for(dim=0; dim<DIM; dim++)

/* IJth:     (i[0],i[1],i[2])-th vector component                       */
/* IJth_ext: (i[0],i[1],i[2])-th vector component in the extended array */

#define IJth(y,i)     ( y[(i[0])+(l_m[0]*((i[1])+(i[2])*l_m[1]))] )
#define IJth_ext(y,i) ( y[(i[0]+1)+((l_m[0]+2)*((i[1]+1)+(i[2]+1)*(l_m[1]+2)))] )

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

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
static int f_local(int Nlocal, realtype t, N_Vector y, 
                   N_Vector ydot, void *f_data);

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data);


static int fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, 
              void *f_dataB);
static int fB_local(int NlocalB, realtype t, 
                    N_Vector y, N_Vector yB, N_Vector yBdot, 
                    void *f_dataB);

static int fQB(realtype t, N_Vector y, N_Vector yB, 
               N_Vector qBdot, void *fQ_dataB);

/*
 *------------------------------------------------------------------
 * Private functions
 *------------------------------------------------------------------
 */

static void SetData(ProblemData d, MPI_Comm comm, int npes, int myId,
                    int *neq, int *l_neq);
static void SetSource(ProblemData d);
static void f_comm( int Nlocal, realtype t, N_Vector y, void *f_data);
static void Load_yext(realtype *src, ProblemData d);
static void PrintHeader();
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
 
  int neq, l_neq;

  void *cvode_mem;
  N_Vector y, q;
  realtype abstol, reltol, abstolQ, reltolQ;
  int mudq, mldq, mukeep, mlkeep;

  void *cvode_memB;
  N_Vector yB, qB;
  realtype abstolB, reltolB, abstolQB, reltolQB;
  int mudqB, mldqB, mukeepB, mlkeepB;

  realtype tret, *qdata, G;

  int ncheckpnt, indexB, flag;

  booleantype output;

  /* Initialize MPI and set Ids */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myId);

  /* Check number of processes */
  npes_needed = NPX * NPY * NPZ;
  MPI_Comm_size(comm, &npes);
  if (npes_needed != npes) {
    if (myId == 0)
      fprintf(stderr,"I need %d processes but I only got %d\n",
              npes_needed, npes);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  /* Test if matlab output is requested */
  if (argc > 1) output = TRUE;
  else          output = FALSE;

  /* Allocate and set problem data structure */
  d = (ProblemData) malloc(sizeof *d);
  SetData(d, comm, npes, myId, &neq, &l_neq);
  
  if (myId == 0) PrintHeader();

  /*-------------------------- 
    Forward integration phase
    --------------------------*/

  /* Allocate space for y and set it with the I.C. */
  y = N_VNew_Parallel(comm, l_neq, neq);
  N_VConst(ZERO, y);
  
  /* Allocate and initialize qB (local contributin to cost) */
  q = N_VNew_Parallel(comm, 1, npes); 
  N_VConst(ZERO, q);

  /* Create CVODES object, attach user data, and allocate space */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeSetUserData(cvode_mem, d);
  abstol = ATOL;  
  reltol = RTOL; 
  flag = CVodeInit(cvode_mem, f, ti, y);
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);

  /* Attach preconditioner and linear solver modules */
  flag = CVSpgmr(cvode_mem, PREC_LEFT, 0);
  mudq = mldq = d->l_m[0]+1;
  mukeep = mlkeep = 2;  
  flag = CVBBDPrecInit(cvode_mem, l_neq, mudq, mldq, 
                       mukeep, mlkeep, ZERO,
                       f_local, NULL);

  /* Initialize quadrature calculations */
  abstolQ = ATOL_Q;
  reltolQ = RTOL_Q;
  flag = CVodeQuadInit(cvode_mem, fQ, q);
  flag = CVodeQuadSStolerances(cvode_mem, reltolQ, abstolQ);
  flag = CVodeSetQuadErrCon(cvode_mem, TRUE);

  /* Allocate space for the adjoint calculation */
  flag = CVodeAdjInit(cvode_mem, STEPS, CV_HERMITE);

  /* Integrate forward in time while storing check points */
  if (myId == 0) printf("Begin forward integration... ");
  flag = CVodeF(cvode_mem, tf, y, &tret, CV_NORMAL, &ncheckpnt);
  if (myId == 0) printf("done. ");

   /* Extract quadratures */
  flag = CVodeGetQuad(cvode_mem, &tret, q);
  qdata = NV_DATA_P(q);
  MPI_Allreduce(&qdata[0], &G, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
  if (myId == 0) printf("  G = %le\n",G);

  /* Print statistics for forward run */
  if (myId == 0) PrintFinalStats(cvode_mem);

  /*-------------------------- 
    Backward integration phase
    --------------------------*/
 
  /* Allocate and initialize yB */
  yB = N_VNew_Parallel(comm, l_neq, neq); 
  N_VConst(ZERO, yB);

  /* Allocate and initialize qB (gradient) */
  qB = N_VNew_Parallel(comm, l_neq, neq); 
  N_VConst(ZERO, qB);

  /* Create and allocate backward CVODE memory */
  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB);
  flag = CVodeSetUserDataB(cvode_mem, indexB, d);
  flag = CVodeInitB(cvode_mem, indexB, fB, tf, yB);
  abstolB = ATOL_B;  
  reltolB = RTOL_B; 
  flag = CVodeSStolerancesB(cvode_mem, indexB, reltolB, abstolB);
  

  /* Attach preconditioner and linear solver modules */
  flag = CVSpgmrB(cvode_mem, indexB, PREC_LEFT, 0); 
  mudqB = mldqB = d->l_m[0]+1;
  mukeepB = mlkeepB = 2;  
  flag = CVBBDPrecInitB(cvode_mem, indexB, l_neq, mudqB, mldqB, 
                        mukeepB, mlkeepB, ZERO, fB_local, NULL);

  /* Initialize quadrature calculations */
  abstolQB = ATOL_QB;
  reltolQB = RTOL_QB;
  flag = CVodeQuadInitB(cvode_mem, indexB, fQB, qB);
  flag = CVodeQuadSStolerancesB(cvode_mem, indexB, reltolQB, abstolQB);
  flag = CVodeSetQuadErrConB(cvode_mem, indexB, TRUE);

  /* Integrate backwards */
  if (myId == 0) printf("Begin backward integration... ");
  flag = CVodeB(cvode_mem, ti, CV_NORMAL);
  if (myId == 0) printf("done.\n");
  
  flag = CVodeGetB(cvode_mem, indexB, &tret, yB);

  /* Print statistics for backward run */
  if (myId == 0) {
    cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, indexB);
    PrintFinalStats(cvode_memB);
  }

   /* Extract quadratures */
  flag = CVodeGetQuadB(cvode_mem, indexB, &tret, qB);

  /* Process 0 collects the gradient components and prints them */
  if (output) {
    OutputGradient(myId, qB, d);
    if (myId == 0) printf("Wrote matlab file 'grad.m'.\n");
  }

  /* Free memory */
  N_VDestroy_Parallel(y);
  N_VDestroy_Parallel(q);
  N_VDestroy_Parallel(qB);
  N_VDestroy_Parallel(yB);

  CVodeFree(&cvode_mem);

  MPI_Finalize();

  return(0);
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
                    int *neq, int *l_neq)
{
  int n[DIM], nd[DIM];
  int dim, size;

  /* Set MPI communicator, id, and total number of processes */

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

  d->xmin[2] = ZMIN;
  d->xmax[2] = ZMAX;
  d->m[2]    = MZ;

  /* Calculate grid spacing and differential volume */

  d->dOmega = ONE;
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

  d->num_procs[2] = NPZ;
  n[2] = NPZ; 
  nd[2] = d->m[2] / NPZ;
  
  /* Compute the neighbors */

  d->nbr_left[0]  = (myId%n[0]) == 0                ? myId : myId-1;
  d->nbr_right[0] = (myId%n[0]) == n[0]-1           ? myId : myId+1;

  d->nbr_left[1]  = (myId/n[0])%n[1] == 0           ? myId : myId-n[0];
  d->nbr_right[1] = (myId/n[0])%n[1] == n[1]-1      ? myId : myId+n[0];

  d->nbr_left[2]  = (myId/n[0]/n[1])%n[2] == 0      ? myId : myId-n[0]*n[1];
  d->nbr_right[2] = (myId/n[0]/n[1])%n[2] == n[2]-1 ? myId : myId+n[0]*n[1];
 
  /* Compute the local subdomains 
     m_start: left border in global index space 
     l_m:     length of the subdomain */
  
  d->m_start[0] = (myId%n[0])*nd[0];
  d->l_m[0]     = d->nbr_right[0] == myId ? d->m[0] - d->m_start[0] : nd[0];

  d->m_start[1] = ((myId/n[0])%n[1])*nd[1];
  d->l_m[1]     = d->nbr_right[1] == myId ? d->m[1] - d->m_start[1] : nd[1];

  d->m_start[2] = (myId/n[0]/n[1])*nd[2];
  d->l_m[2]     = d->nbr_right[2] == myId ? d->m[2] - d->m_start[2] : nd[2];

  /* Allocate memory for the y_ext array 
     (local solution + data from neighbors) */

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
      for(i[2]=0; i[2]<l_m[2]; i[2]++) {
        x[2] = xmin[2] + (m_start[2]+i[2]) * dx[2];
        
        g = G1_AMPL 
          * EXP( -SQR(G1_X-x[0])/SQR(G1_SIGMA) ) 
          * EXP( -SQR(G1_Y-x[1])/SQR(G1_SIGMA) )
          * EXP( -SQR(G1_Z-x[2])/SQR(G1_SIGMA) ); 
        
        g += G2_AMPL 
          * EXP( -SQR(G2_X-x[0])/SQR(G2_SIGMA) ) 
          * EXP( -SQR(G2_Y-x[1])/SQR(G2_SIGMA) )
          * EXP( -SQR(G2_Z-x[2])/SQR(G2_SIGMA) ); 
        
        if( g < G_MIN ) g = ZERO;

        IJth(pdata, i) = g;
      }
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

static void f_comm(int N_local, realtype t, N_Vector y, void *f_data)
{
  int id, n[DIM], proc_cond[DIM], nbr[DIM][2];
  ProblemData d;
  realtype *yextdata, *ydata;
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
  yextdata = d->y_ext;
  ydata    = NV_DATA_P(y);
  
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
  proc_cond[2] = (id/n[0]/n[1])%2;

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
        l[1]=(dim+2)%DIM;
        for(i[l[1]]=0; i[l[1]]<l_m[l[1]]; i[l[1]]++) 
          for(i[l[0]]=0; i[l[0]]<l_m[l[0]]; i[l[0]]++) 
            buf_send[c++] = IJth(ydata, i);
	  
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

        /* Compute the index of the boundary (right or left) in yextdata */
        i[dim] = (dir ^ proc_cond[dim]) ? l_m[dim] : -1;

        /* Loop over all other dimensions and copy data into yextdata */
        for(i[l[1]]=0; i[l[1]]<l_m[l[1]]; i[l[1]]++)
          for(i[l[0]]=0; i[l[0]]<l_m[l[0]]; i[l[0]]++)
            IJth_ext(yextdata, i) = buf_recv[c++];
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

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
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

  return(0);
}

static int f_local(int Nlocal, realtype t, N_Vector y, 
                   N_Vector ydot, void *f_data)
{
  realtype *Ydata, *dydata, *pdata;
  realtype dx[DIM], c, v[DIM], cl[DIM], cr[DIM];
  realtype adv[DIM], diff[DIM];
  realtype xmin[DIM], xmax[DIM], x[DIM], x1;
  int i[DIM], l_m[DIM], m_start[DIM], nbr_left[DIM], nbr_right[DIM], id;
  ProblemData d;
  int dim;

  d = (ProblemData) f_data;

  /* Extract stuff from data structure */
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

  /* Get pointers to vector data */
  dydata = NV_DATA_P(ydot);
  pdata  = NV_DATA_P(d->p);

  /* Copy local segment of y to y_ext */
  Load_yext(NV_DATA_P(y), d);
  Ydata = d->y_ext;

  /* Velocity components in x1 and x2 directions (Poiseuille profile) */
  v[1] = ZERO;
  v[2] = ZERO;

  /* Local domain is [xmin+(m_start+1)*dx, xmin+(m_start+1+l_m-1)*dx] */

  for(i[2]=0; i[2]<l_m[2]; i[2]++) {

    x[2] = xmin[2] + (m_start[2]+i[2])*dx[2];

    for(i[1]=0; i[1]<l_m[1]; i[1]++) {

      x[1] = xmin[1] + (m_start[1]+i[1])*dx[1];

      /* Velocity component in x0 direction (Poiseuille profile) */
      x1 = x[1] - xmin[1] - L;
      v[0] = V_COEFF * (L + x1) * (L - x1);

      for(i[0]=0; i[0]<l_m[0]; i[0]++) {

        x[0] = xmin[0] + (m_start[0]+i[0])*dx[0];

        c  = IJth_ext(Ydata, i);	       

        /* Source term*/
        IJth(dydata, i) = IJth(pdata, i);

        FOR_DIM {
          i[dim]+=1;
          cr[dim] = IJth_ext(Ydata, i);
          i[dim]-=2;
          cl[dim] = IJth_ext(Ydata, i);
          i[dim]+=1;

          /* Boundary conditions for the state variables */
          if( i[dim]==l_m[dim]-1 && nbr_right[dim]==id)
            cr[dim] = cl[dim];
          else if( i[dim]==0 && nbr_left[dim]==id )
            cl[dim] = cr[dim];

          adv[dim]  = v[dim] * (cr[dim]-cl[dim]) / (TWO*dx[dim]);
          diff[dim] = DIFF_COEF * (cr[dim]-TWO*c+cl[dim]) / SQR(dx[dim]);

          IJth(dydata, i) += (diff[dim] - adv[dim]);
        } 
      }
    }
  }

  return(0);
}

/*
 *------------------------------------------------------------------
 * fQ:
 * Right-hand side of quadrature equations on forward integration.
 * The only quadrature on this phase computes the local contribution
 * to the function G.
 *------------------------------------------------------------------
 */

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *fQ_data)
{
  ProblemData d;
  realtype *dqdata;

  d = (ProblemData) fQ_data;

  dqdata = NV_DATA_P(qdot);

  dqdata[0] = N_VDotProd_Parallel(y,y);
  dqdata[0] *= RCONST(0.5) * (d->dOmega);

  return(0);
}

/*
 *------------------------------------------------------------------
 * fB and fB_local:
 * Backward phase ODE right-hand side (the discretized adjoint PDE)
 *------------------------------------------------------------------
 */

static int fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, 
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

  return(0);
}

static int fB_local(int NlocalB, realtype t, 
                    N_Vector y, N_Vector yB, N_Vector dyB, 
                    void *f_dataB)
{
  realtype *YBdata, *dyBdata, *ydata;
  realtype dx[DIM], c, v[DIM], cl[DIM], cr[DIM];
  realtype adv[DIM], diff[DIM];
  realtype xmin[DIM], xmax[DIM], x[DIM], x1;
  int i[DIM], l_m[DIM], m_start[DIM], nbr_left[DIM], nbr_right[DIM], id;
  ProblemData d;
  int dim;
  
  d = (ProblemData) f_dataB;

  /* Extract stuff from data structure */
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
 
  dyBdata = NV_DATA_P(dyB);
  ydata   = NV_DATA_P(y);

  /* Copy local segment of yB to y_ext */
  Load_yext(NV_DATA_P(yB), d);
  YBdata = d->y_ext;

  /* Velocity components in x1 and x2 directions (Poiseuille profile) */
  v[1] = ZERO;
  v[2] = ZERO;
 
  /* local domain is [xmin+(m_start)*dx, xmin+(m_start+l_m-1)*dx] */

  for(i[2]=0; i[2]<l_m[2]; i[2]++) {

    x[2] = xmin[2] + (m_start[2]+i[2])*dx[2];
    
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      
      x[1] = xmin[1] + (m_start[1]+i[1])*dx[1];
	  
      /* Velocity component in x0 direction (Poiseuille profile) */
      x1 = x[1] - xmin[1] - L;
      v[0] = V_COEFF * (L + x1) * (L - x1);

      for(i[0]=0; i[0]<l_m[0]; i[0]++) {

        x[0] = xmin[0] + (m_start[0]+i[0])*dx[0];
        
        c  = IJth_ext(YBdata, i);	       
        
        /* Source term for adjoint PDE */
        IJth(dyBdata, i) = -IJth(ydata, i);
        
        FOR_DIM {
          
          i[dim]+=1;
          cr[dim] = IJth_ext(YBdata, i);
          i[dim]-=2;
          cl[dim] = IJth_ext(YBdata, i);
          i[dim]+=1;

          /* Boundary conditions for the adjoint variables */
          if( i[dim]==l_m[dim]-1 && nbr_right[dim]==id)
	    cr[dim] = cl[dim]-(TWO*dx[dim]*v[dim]/DIFF_COEF)*c;
          else if( i[dim]==0 && nbr_left[dim]==id )
	      cl[dim] = cr[dim]+(TWO*dx[dim]*v[dim]/DIFF_COEF)*c;
		  
          adv[dim]  = v[dim] * (cr[dim]-cl[dim]) / (TWO*dx[dim]);
          diff[dim] = DIFF_COEF * (cr[dim]-TWO*c+cl[dim]) / SQR(dx[dim]);
          
          IJth(dyBdata, i) -= (diff[dim] + adv[dim]);
        } 
      }
    }
  }

  return(0);
}

/*
 *------------------------------------------------------------------
 * fQB:
 * Right-hand side of quadrature equations on backward integration
 * The i-th component of the gradient is nothing but int_t yB_i dt
 *------------------------------------------------------------------
 */

static int fQB(realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, 
               void *fQ_dataB)
{
  ProblemData d;

  d = (ProblemData) fQ_dataB;

  N_VScale_Parallel(-(d->dOmega), yB, qBdot);

  return(0);
}

/*
 *------------------------------------------------------------------
 * Load_yext: 
 * copies data from src (y or yB) into y_ext, which already contains
 * data from neighboring processes.
 *------------------------------------------------------------------
 */

static void Load_yext(realtype *src, ProblemData d)
{
  int i[DIM], l_m[DIM], dim;

  FOR_DIM l_m[dim] = d->l_m[dim];
     
  /* copy local segment */
  for  (i[2]=0; i[2]<l_m[2]; i[2]++)
    for(i[1]=0; i[1]<l_m[1]; i[1]++)
      for(i[0]=0; i[0]<l_m[0]; i[0]++)
	IJth_ext(d->y_ext, i) = IJth(src, i);
}

/*
 *------------------------------------------------------------------
 * PrintHeader:
 * Print first lins of output (problem description)
 *------------------------------------------------------------------
 */

static void PrintHeader()
{
    printf("\nParallel Krylov adjoint sensitivity analysis example\n");
    printf("%1dD Advection diffusion PDE with homogeneous Neumann B.C.\n",DIM);
    printf("Computes gradient of G = int_t_Omega ( c_i^2 ) dt dOmega\n");
    printf("with respect to the source values at each grid point.\n\n");

    printf("Domain:\n");

    printf("   %f < x < %f   mx = %d  npe_x = %d \n",XMIN,XMAX,MX,NPX);
    printf("   %f < y < %f   my = %d  npe_y = %d \n",YMIN,YMAX,MY,NPY);
    printf("   %f < z < %f   mz = %d  npe_z = %d \n",ZMIN,ZMAX,MZ,NPZ);

    printf("\n");
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

  flag = CVSpilsGetWorkSpace(cvode_mem, &lenrwSPGMR, &leniwSPGMR);
  flag = CVSpilsGetNumLinIters(cvode_mem, &nli);
  flag = CVSpilsGetNumPrecEvals(cvode_mem, &npe);
  flag = CVSpilsGetNumPrecSolves(cvode_mem, &nps);
  flag = CVSpilsGetNumConvFails(cvode_mem, &ncfl);
  flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeSPGMR);

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

  /* Write matlab files with solutions from each process */

  for(i[0]=0; i[0]<l_m[0]; i[0]++) {
    x[0] = xmin[0] + (m_start[0]+i[0]) * dx[0];
    fprintf(fid,"x%d(%d,1) = %le; \n",  myId, i[0]+1, x[0]);
  }
  for(i[1]=0; i[1]<l_m[1]; i[1]++) {
    x[1] = xmin[1] + (m_start[1]+i[1]) * dx[1];
    fprintf(fid,"y%d(%d,1) = %le; \n",  myId, i[1]+1, x[1]);
  }
  for(i[2]=0; i[2]<l_m[2]; i[2]++) {
    x[2] = xmin[2] + (m_start[2]+i[2]) * dx[2];
    fprintf(fid,"z%d(%d,1) = %le; \n",  myId, i[2]+1, x[2]);
  }

  for(i[0]=0; i[0]<l_m[0]; i[0]++) {
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      for(i[2]=0; i[2]<l_m[2]; i[2]++) {
        g = IJth(qBdata, i);
        p = IJth(pdata, i);
        fprintf(fid,"p%d(%d,%d,%d) = %le; \n", myId, i[1]+1, i[0]+1, i[2]+1, p);
        fprintf(fid,"g%d(%d,%d,%d) = %le; \n", myId, i[1]+1, i[0]+1, i[2]+1, g);
      }
    }
  }
  fclose(fid);

  /* Write matlab driver */

  if (myId == 0) {

    fid = fopen("grad.m","w");

    fprintf(fid,"clear;\n\n");
    for (ip=0; ip<d->npes; ip++) {
      fprintf(fid,"grad%03d;\n",ip);
    }
    fprintf(fid,"\n\n");    

    fprintf(fid,"figure;\nhold on\n");
    fprintf(fid,"trans = 0.7;\n");
    fprintf(fid,"ecol  = 'none';\n");
    fprintf(fid,"xp=[%f %f];\n",G1_X,G2_X);
    fprintf(fid,"yp=[%f %f];\n",G1_Y,G2_Y);
    fprintf(fid,"zp=[%f %f];\n",G1_Z,G2_Z);
    fprintf(fid,"\n\n");

    for (ip=0; ip<d->npes; ip++) {
      fprintf(fid,"[X,Y,Z]=meshgrid(x%d,y%d,z%d);\n",ip,ip,ip);
      fprintf(fid,"s%d=slice(X,Y,Z,g%d,xp,yp,zp);\n",ip,ip);
      fprintf(fid,"for i = 1:length(s%d)\n",ip);
      fprintf(fid,"  set(s%d(i),'FaceAlpha',trans);\n",ip);
      fprintf(fid,"  set(s%d(i),'EdgeColor',ecol);\n",ip);
      fprintf(fid,"end\n");
    }
    
    fprintf(fid,"view(3)\n");
    fprintf(fid,"\nshading interp\naxis equal\n\n\n");

    fprintf(fid,"figure;\nhold on\n");
    fprintf(fid,"iso=[16.0 8.0 4.0 2.0];\n");
    fprintf(fid,"trans=[1.0 0.8 0.4 0.2];\n");
    fprintf(fid,"fcol={'red','blue','green','yellow'};\n");
    fprintf(fid,"\n\n");

    fprintf(fid,"for i=1:4\n");
    for (ip=0; ip<d->npes; ip++) {
      fprintf(fid,"  [X,Y,Z]=meshgrid(x%d,y%d,z%d);\n",ip,ip,ip);
      fprintf(fid,"  p = patch(isosurface(X,Y,Z,g%d,iso(i)));\n",ip);
      fprintf(fid,"  isonormals(X,Y,Z,g%d,p);\n",ip);
      fprintf(fid,"  set(p, 'FaceColor', fcol{i}, 'EdgeColor', 'none', 'FaceAlpha', trans(i));\n");
      fprintf(fid,"  daspect([1 1 1]);\n");
    }
    fprintf(fid,"end\n");

    fprintf(fid,"view(3)\n");
    fprintf(fid,"camlight; lighting phong;\n\n");

    fclose(fid);
  }
}
