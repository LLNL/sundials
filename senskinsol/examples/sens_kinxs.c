 /**********************************************************************
  *                                                                    *
  * File: sens_kinxp.c                                                 *
  *                                                                    *
  * Programmers: Allan G. Taylor, Alan C. Hindmarsh, and               *
  *              Keith E. Grant @ LLNL                                 *
  *                                                                    *
  * Version of 03 Jan 2001                                             *
  *--------------------------------------------------------------------*
  *                                                                    *
  *                                                                    *
  * Example problem for Sens_KINSOL, serial and parallel machine       *
  * cases. This example solves a nonlinear system that arises from     *
  * a system of partial differential equations. The PDE system is a    *
  * food web population model, with predator-prey interaction and      *
  * diffusion on the unit square in two dimensions. The dependent      *
  * variable vector is                                                 *
  *                                                                    *
  *       1   2         ns                                             *
  * c = (c , c ,  ..., c  )              (denoted by the variable cc)  *
  *                                                                    *
  * and the pde's are as follows:                                      *
  *                                                                    *
  *                    i       i                                       *
  *         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)     *
  *                    xx      yy         i                            *
  *                                                                    *
  *   where                                                            *
  *                                                                    *
  *                   i             ns         j                       *
  *   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )                      *
  *    i                           j=1                                 *
  *                                                                    *
  * The number of species is ns = 2 * np, with the first np being      *
  * prey and the last np being predators. The number np is both the    *
  * number of prey and predator species. The coefficients a(i,j),      *
  * b(i), d(i) are                                                     *
  *                                                                    *
  *   a(i,i) = -P0*AA  (all i)                                         *
  *   a(i,j) = -P1*GG  (i <= np , j >  np)                             *
  *   a(i,j) =  P2*EE  (i >  np,  j <= np)                             *
  *   b(i) = BB * (1 + alpha * x * y)  (i <= np)                       *
  *   b(i) =-BB * (1 + alpha * x * y)  (i >  np)                       *
  *   d(i) = DPREY  (i <= np)                                          *
  *   d(i) = DPRED  ( i > np)                                          *
  *                                                                    *
  *  The various scalar parameters are set using define's              *
  *  or in routine InitUserData.                                       * 
  *                                                                    *
  *  The boundary conditions are .. normal derivative  =  0.           *
  *                                                                    *
  *  The initial guess is constant in x and y, although the final      *
  *  solution is not.                                                  *
  *                                                                    *
  *  The PDEs are discretized by central differencing on a mx by       *
  *  my mesh.                                                          *
  *                                                                    *
  *  The nonlinear system is solved by KINSOL using the method         *
  *  specified in local variable globalstrat .                         *
  *                                                                    *
  *  The preconditioner matrix is a block-diagonal matrix based on the *
  *  partial derivatives of the interaction terms f (in the above      *
  *  equation) only                                                    *
  *                                                                    *
  *                                                                    *
  *  Execution:                                                        *
  *                                                                    *
  *     Parallel:                                                      *
  *        mpirun -np N -machinefile machines sens_kinxp               *
  *        {with N = NPEX*NPEY, total number of processors, see below} *
  *                                                                    *
  *     Serial:                                                        *
  *        sens_kinxs                                                  *
  *                                                                    *
  *                                                                    *
  *  References...                                                     *
  *                                                                    *
  * 1.                                                                 *
  *  Peter N Brown and Youcef Saad, Hybrid Krylov Methods for          *
  *  Nonlinear Systems of Equations LLNL report UCRL-97645,            *
  *  November 1987.                                                    *
  *                                                                    *
  * 2.                                                                 *
  *  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage Matrix      *
  *  Methods in Stiff ODE systems, Lawrence Livermore National         *
  *  Laboratory Report  UCRL-95088, Rev. 1, June 1987, and Journal of  *
  *  Applied Mathematics and Computation, Vol. 31 (May 1989),          *
  *  pp. 40-91. ( for a description of the time-dependent version of   *
  *  this test problem.)                                               *
  *                                                                    *
  *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"    /* definitions of real, integer, boole, TRUE, FALSE */
#include "llnlmath.h"    /* contains RSqrt and UnitRoundoff routines         */
#include "nvector.h"     /* definitions of type N_Vector, macro N_VDATA      */
#include "iterativ.h"    /* contains the enum for types of preconditioning   */
#include "smalldense.h"  /* use generic DENSE solver for preconditioning     */
#include "kinspgmr.h"    /* use KINSpgmr linear solver                       */
#include "kinsol.h"      /* main KINSOL header file                          */
#include "sens_kinsol.h" /* Sensitivity parameter definitions                */
#include "sens_kinspgmr.h" /* Definitions for the sens. linear solver API    */

#ifdef PARALLEL
#include "mpi.h"         /* MPI include file */
#else
typedef int  MPI_Comm;   /* Defined for call compatibility */
#endif

/***********************************************************************
 *
 * Problem Constants
 *
 * The following definitions are made to keep the maximum similarity
 * of structure between the serial and parallel cases. For the serial
 * case MXSUB and MYSUB are the full sizes of the grid. Using them
 * makes several later definitions identical for the serial and
 * parallel cases.
 *
 **********************************************************************/

#define NUM_SPECIES  6     /* Must equal 2*(number of prey or predators)
                              number of prey = number of predators         */

#ifdef PARALLEL

#define NPEX      2        /* number of processors in the x-direction      */
#define NPEY      2        /* number of processors in the y-direction      */
#define MXSUB     10       /* MXSUB = number of x mesh points per subgrid  */
#define MYSUB     10       /* MYSUB = number of y mesh points per subgrid  */
#define MX        (NPEX*MXSUB) /* Number of grid points in the x-direction */
#define MY        (NPEY*MYSUB) /* Number of grid points in the y-direction */

#else

#define NPEX      1        /* For the serial case, this must always be one */
#define NPEY      1        /* For the serial case, this must always be one */
#define MXSUB     4        /* MXSUB = number of x mesh points (full grid)  */
#define MYSUB     4        /* MYSUB = number of y mesh points (full grid)  */

#endif

#define MX        (NPEX*MXSUB) /* Number of grid points in the x-direction */
#define MY        (NPEY*MYSUB) /* Number of grid points in the y-direction */

#define NEQ       (NUM_SPECIES * MX * MY)    /* Number of equations */
#define NSMXSUB   (NUM_SPECIES * MXSUB)
#define NSMXSUB2  (NUM_SPECIES * (MXSUB+2))

#define NPSENS    3              /* Number of sensitivity parameters */

#define PI        RCONST(3.1415926535898)   /* pi */
#define AA        RCONST(1.0)    /* value of coefficient a, above eqns */
#define EE        RCONST(10000.) /* value of coefficient e, above eqns */
#define GG        RCONST(0.5e-6) /* value of coefficient g, above eqns */
#define BB        RCONST(1.0)    /* value of coefficient b, above eqns */
#define DPREY     RCONST(1.0)    /* value of coefficient dprey, above eqns */
#define DPRED     RCONST(0.5)    /* value of coefficient dpred, above eqns */
#define ALPHA     RCONST(1.0)    /* value of coefficient alpha, above eqns */
#define AX        RCONST(1.0)    /* total range of x variable */
#define AY        RCONST(1.0)    /* total range of y variable */
#define FTOL      RCONST(1.e-7)  /*  ftol tolerance */
#define STOL      RCONST(1.e-13) /*  stol tolerance */
#define THOUSAND  RCONST(1000.0) /* one thousand */
#define ZERO      RCONST(0.)     /* 0. */
#define ONE       RCONST(1.0)    /* 1. */


/* Sensitivity Parameter (SP) access enumerations */
enum {
   SP_SELFSELF,         /* A Species interacting with itself */
   SP_PREYPRED,         /* Prey interacting with predators   */
   SP_PREDPREY          /* Predators interacting with prey   */
};


/***********************************************************************
 *
 * User-defined vector accessor macro: IJ_Vptr
 *
 * IJ_Vptr is defined in order to isolate the underlying 3-d structure
 * of the dependent variable vector from its underlying 1-d storage
 * (an N_Vector). IJ_Vptr returns a pointer to the location in vv
 * corresponding to ns = 0 ,  jx = i,  jy = j .
 *
 **********************************************************************/

#define IJ_Vptr(vv,i,j)   (&(((vv)->data)[(i)*NUM_SPECIES + (j)*NSMXSUB]))


#ifndef Parallel

/***********************************************************************
 *
 * User-defined vector and matrix accessor macros: IJKth, IJth
 *

   IJKth(vdata,i,j,k) references the element in the vdata array for
 * species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
 * 0 <= j <= MX-1, 0 <= k <= MY-1. The vdata array is obtained via
 * the macro call vdata = N_VDATA(v), where v is an N_Vector.
 * For each mesh point (j,k), the elements for species i and i+1 are
 * contiguous within vdata.
 *
 * IJth(a,i,j) references the (i,j)th entry of the small matrix
 * real **a, where 1 <= i,j <= NUM_SPECIES. The small matrix routines
 * in smalldense.h work with matrices stored by column in a
 * 2-dimensional array. In C, arrays are indexed starting at 0, not 1.
 *
 **********************************************************************/

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMXSUB])
#define IJth(a,i,j)        (a[j-1][i-1])

#endif


/* Type : UserData
   contains preconditioner blocks, pivot arrays, and problem constants */

typedef struct {

  real      ax, ay, dx, dy;
  real      uround, sqruround;
  integer   Neq, mx, my, ns, np;
  real    **acoef, *bcoef;
  real     *cox, *coy;
  N_Vector  rates;
  boole     sens_active;      /* Sensitivity active flag */
  real     *sens_p;           /* Sensitivity parameters  */
  real    **P[MXSUB][MYSUB];
  integer  *pivot[MXSUB][MYSUB];
  real      cext[NUM_SPECIES * (MXSUB+2)*(MYSUB+2)];

  integer   my_pe;
  integer   isubx;
  integer   isuby;
  integer   nsmxsub;
  integer   nsmxsub2;
  MPI_Comm  comm;

} *UserData;


/* Private Helper Functions */

static UserData AllocUserData(void);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector cc, N_Vector sc);
static void PrintFinalStats(long int *iopt);
static void WebRate(real xx, real yy, real *cxy, real *ratesxy, void *f_data);
static real DotProd(integer size, real *x1, real *x2);
static void InitUserData(integer my_pe, MPI_Comm comm, real *sens_p,
               UserData data);
static void PrintOutput(integer my_pe, MPI_Comm comm, N_Vector cc);
static void PrintSensOut(integer my_pe, MPI_Comm comm, N_Vector ww);
static void fcalcprpr(integer Neq, N_Vector cc, N_Vector fval,
               void *f_data);

#ifdef PARALLEL

static void BSend(MPI_Comm comm, integer my_pe, integer isubx, integer isuby,
                  integer dsizex, integer dsizey, real *cdata);
static void BRecvPost(MPI_Comm comm, MPI_Request request[], integer my_pe,
		      integer isubx, integer isuby,
		      integer dsizex, integer dsizey,
		      real *cext, real *buffer);
static void BRecvWait(MPI_Request request[], integer isubx, integer isuby,
		      integer dsizex, real *cext, real *buffer);
static void ccomm(integer Neq, real *cdata, UserData data);

#endif



/* Functions Called by the KINSOL Solver */

static void funcprpr(integer Neq, N_Vector cc, N_Vector fval,
               void *f_data);


static int Precondbd(integer Neq, N_Vector cc, N_Vector cscale,
         N_Vector fval, N_Vector fscale, N_Vector vtemp1,N_Vector vtemp2,
         SysFn func, real uround, long int *nfePtr, void *P_data);


static int PSolvebd(integer Neq, N_Vector cc, N_Vector cscale,
        N_Vector fval, N_Vector fscale, N_Vector vtem, N_Vector ftem,
        SysFn func, real uround, long int *nfePtr, void *P_data);


/***********************************************************************
 *
 * Main Program of Preditor - Prey Example
 *
 **********************************************************************/

int  main(int argc, char *argv[])

{
  FILE *msgfile;
  integer Neq=NEQ;
  integer globalstrategy, i;
  real fnormtol, scsteptol, ropt[OPT_SIZE];
  long int iopt[OPT_SIZE];
  N_Vector cc, sc, constraints;
  N_Vector ww;
  int flag;
  boole optIn;

  integer  ip;
#if ( NPSENS > 0 )
  real psens[NPSENS];         /* Sensitivity parameters */
  real psbar[NPSENS];         /* Sens parameter nominal magnitudes */
#else
  real *psens = NULL;
  real *psbar = NULL;
#endif

  void *sens_params;
  SensKINParams sens;
  UserData data;

  machEnvType machEnv;
  int my_pe;
  int npes;
  int npelast = NPEX*NPEY-1;
  MPI_Comm  comm;

  const int kryldim =  16; /* Maximum Krylov dimension      */
  const int maxlrst = 128; /* Maximum no of linear restarts */

#ifdef PARALLEL
  integer local_N;
  const char mymode[]="parallel";
#else
  const char mymode[]="serial";
#endif

  /* Allocate memory, and set problem data, initial values, tolerances */


#ifdef PARALLEL

  msgfile = fopen("sens_kinxp.errout","w");


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

#else

  msgfile = fopen("sens_kinxs.errout","w");

  /* Definitions for compatibility with parallel coding */
  my_pe = 0;
  npes  = 1;
  comm  = 0;

#endif


#ifdef PARALLEL
  local_N = NUM_SPECIES*MXSUB*MYSUB;   /* Set local length */
  machEnv = PVecInitMPI(comm, local_N, Neq, &argc, &argv);
  if(machEnv==NULL) return(1);
#else
  machEnv = NULL;
#endif

  /* Allocate and initialze sensitivity data block */

  for ( ip=0; ip<NPSENS; ip++) {
    psens[ip] = ONE;
    psbar[ip] = ONE;
  }
  sens_params = SensKINMalloc(Neq, NPSENS, DIFF_FRWRD, ONE, psens, psbar,
         msgfile, machEnv);
  sens = (SensKINParams) sens_params;
  if (my_pe==0 && sens == NULL) {
     printf("\nSensKINMalloc failed.\n\n");
     return(1);
  }

  /* allocate and initialize user data block */

  data=(UserData)AllocUserData();
  InitUserData(my_pe, comm, psens, data);

  /* example of changing defaults using iopt */
  optIn = TRUE;
  for(i=0;i<KINSOL_IOPT_SIZE;i++)iopt[i]=0;
  for(i=0;i<KINSOL_ROPT_SIZE;i++)ropt[i]=ZERO;
  iopt[MXITER]=250;
  iopt[PRINTFL]=3;

  /* choose global strategy */
  globalstrategy = INEXACT_NEWTON;

  /* allocate (initialize) vectors */
  cc = N_VNew(Neq, machEnv);
  sc = N_VNew(Neq, machEnv);
  data->rates=N_VNew(Neq,machEnv);

  constraints = N_VNew(Neq, machEnv);
  N_VConst(0.,constraints);

  SetInitialProfiles(cc, sc);

  fnormtol=FTOL; scsteptol=STOL;



  /*
   * Call KINSpgmr to specify the KINSOL linear solver KINSpgmr with
   * solve routines Precondbd and PSolvebd, and the pointer to the
   * user-defined  block data.
   */

  flag = SensKINSpgmr(sens_params,
      kryldim, /* a zero in this position forces use of the KINSpgmr
                  default for maxl, dimension of the Krylov space */
      maxlrst, /* if zero in this position forces use of the KINSpgmr
                  default for maxlrst, the max number of linear solver
                  restarts allowed */
      0,    /* a zero in this position forces use of the KINSpgmr
               default for msbpre, the number of calls to the
               preconditioner allowed without a call to the
               preconditioner setup routine */
      Precondbd, /* user-supplied preconditioner setup routine */
      PSolvebd,  /* user-supplied preconditioner solve routine */
	   NULL,   /* user-supplied ATimes routine -- Null chosen here */
      data);

  if(flag!=0){
    if(my_pe==0)printf("KINSpgmr returned nonzero (failed). %d\n",flag);
    return(1);
  }

  if(my_pe==0) printf(" \n Predator-prey test problem --  "
      "KINSol (%s version)\n\n",mymode);

  /* first,print out the problem size and then the
     initial concentration profile */

  if(my_pe==0){
    printf("Mesh dimensions %d X %d\n",MX,MY);
    printf("Number of species %d \n",NUM_SPECIES);
    printf("Total system size %d\n",Neq);
    printf("Preconditioning uses interaction-only block-diagonal matrix\n");
    printf("tolerance parameters:  fnormtol = %g   scsteptol = %g\n",
	           fnormtol,scsteptol);
    printf("Max Krylov dimension = %d, Max linear restarts = %d\n",
       kryldim, maxlrst);

    printf("\nInitial profile of concentration\n");
  }

  if(my_pe==0 || my_pe==npelast)  PrintOutput(my_pe, comm, cc);


  /* call KINSol and print output concentration profile */

  flag = SensKINSol(sens_params,  /* Sensitivity memory block */
		Neq,            /* system size -- number of equations  */
		cc,             /* solution cc of funcprpr(cc)=0       */
		funcprpr,       /* function for the system equations   */
		globalstrategy, /* global stragegy choice              */
		sc,             /* scaling vector, for the variable cc */
		sc,             /* scaling vector for function values  */
		fnormtol,       /* tolerance on fnorm funcprpr(cc)     */
		scsteptol,      /* step size tolerance                 */
		constraints,    /* constraints vector                  */
		optIn,          /* optional inputs flat: TRUE or FALSE */
		iopt,           /* integer optional input array        */
		ropt,           /* real optional input array           */
		data);          /* pointer to user data                */

  if(my_pe==0) {
    if (flag != KINSOL_SUCCESS) {
      printf("\nKINSol failed, flag=%d.\n\n", flag);
      return(flag); }

      printf("\n\n\nComputed equilibrium species concentrations:\n\n");
  }

  if(my_pe==0 || my_pe==npelast) PrintOutput(my_pe, comm, cc);

  /* cc values are available on each processor */
  if(my_pe==0) PrintFinalStats(iopt);

#if ( NPSENS > 0 )

 /**********************************************************************
  *
  *             *****  Start of Sensitivity Solutions *****
  *
  *********************************************************************/
  data->sens_active = TRUE;
  flag = SensKINLinInit (sens_params);
  if (my_pe==0 && flag != 0) {
     printf("\nSensKINInit failed and returned %d\n\n",flag);
     return (flag);
  }

  ww = N_VNew(Neq, machEnv);

  for (ip=0; ip<NPSENS; ip++) {

    if ( my_pe == 0 )
       printf ("\nSolution for sensitivity parameter %d\n", ip);


    /* Optional call to reset difference option of scale factor */

     flag = SensKINDiff (sens_params, DIFF_FRWRD, ONE);
     if (my_pe==0 && flag != 0) {
        printf("\nSensKINDiff failed and returned %d\n\n",flag);
        return (flag);
     }

     flag = SensKINLinSolve (sens_params, ip, ww);

     if (my_pe==0 && flag != 0) {
        printf("\nSensKINLinSolve failed and returned %d\n\n",flag);
        return (flag);
     }

     if(my_pe==0 || my_pe==npelast) PrintSensOut(my_pe, comm, ww);
     if(my_pe==0) PrintFinalStats(iopt);

  }

  data->sens_active = FALSE;

#endif

  /* Free memory and print final statistics */
  N_VFree(cc);
  N_VFree(sc);
  N_VFree(constraints);
  N_VFree(ww);
  FreeUserData(data);
  SensKINFree(sens);

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return(0);
}


/***********************************************************************
 *                                                                     *
 ********************** Private Helper Functions ***********************
 *                                                                     *
 **********************************************************************/


/***********************************************************************
 *
 * Function: AllocUserData
 *
 * Allocate memory for data structure of type UserData
 *
 **********************************************************************/

static UserData AllocUserData(void)
{
  int jx, jy;
  UserData data;

  data = (UserData) malloc(sizeof *data);

  for (jx=0; jx < MXSUB; jx++) {
    for (jy=0; jy < MYSUB; jy++) {
      (data->P)[jx][jy] = denalloc(NUM_SPECIES);
      (data->pivot)[jx][jy] = denallocpiv(NUM_SPECIES);
    }
  }
 (data->acoef) = denalloc(NUM_SPECIES);
 (data->bcoef) = (real *)malloc(NUM_SPECIES * sizeof(real));
 (data->cox)   = (real *)malloc(NUM_SPECIES * sizeof(real));
 (data->coy)   = (real *)malloc(NUM_SPECIES * sizeof(real));


  return(data);
}


/* readability constants defined */

#define acoef  (data->acoef)
#define bcoef  (data->bcoef)
#define cox    (data->cox)
#define coy    (data->coy)


/***********************************************************************
 *
 * Function: InitUserData
 *
 * Load problem constants into user data structure
 *
 **********************************************************************/

/* Load problem constants in data */

static void InitUserData(integer my_pe, MPI_Comm comm, real *sens_p,
   UserData data)
{
  int i, j, np;
  real *a1,*a2, *a3, *a4, dx2, dy2;

  data->mx  = MX;
  data->my  = MY;
  data->ns  = NUM_SPECIES;
  data->np  = NUM_SPECIES / 2;
  data->ax  = AX;
  data->ay  = AY;
  data->dx  = (data->ax)/(MX-1);
  data->dy  = (data->ay)/(MY-1);
  data->Neq = NEQ;

  data->my_pe    = my_pe;
  data->isuby    = my_pe / NPEX;
  data->isubx    = my_pe - data->isuby*NPEX;
  data->nsmxsub  = NSMXSUB;
  data->nsmxsub2 = NSMXSUB2;

  data->uround = UnitRoundoff();
  data->sqruround = RSqrt(data->uround);
  data->sens_active = FALSE;  /* Sensitivity analysis active flag */
  data->sens_p = sens_p;      /* Sensitivity analysis parameters  */

  data->comm = comm;

  /* set up the coefficients a and b plus others found in the equations */
  np = data->np;

  dx2=(data->dx)*(data->dx); dy2=(data->dy)*(data->dy);

  for(i=0;i<np;i++){
    a1= &(acoef[i][np]);
    a2= &(acoef[i+np][0]);
    a3= &(acoef[i][0]);
    a4= &(acoef[i+np][np]);
    /*  fill in the portion of acoef in the four quadrants, row by row */
    for(j=0;j<np;j++){
      *a1++ =  -GG;
      *a2++ =   EE;
      *a3++ = ZERO;
      *a4++ = ZERO;
    }

    /* and then change the diagonal elements of acoef to -AA */
    acoef[i][i]=-AA;
    acoef[i+np][i+np] = -AA;

    bcoef[i] = BB;
    bcoef[i+np] = -BB;

    cox[i]=DPREY/(dx2);
    cox[i+np]=DPRED/(dx2);

    coy[i]=DPREY/(dy2);
    coy[i+np]=DPRED/(dy2);
  }

}


/***********************************************************************
 *
 * Function: FreeUserData
 *
 * Free user data structure memory
 *
 **********************************************************************/

static void FreeUserData(UserData data)
{
  int jx, jy;

  for (jx=0; jx < MXSUB; jx++) {
    for (jy=0; jy < MYSUB; jy++) {
      denfree((data->P)[jx][jy]);
      denfreepiv((data->pivot)[jx][jy]);
    }
  }

  denfree(acoef);
  free(bcoef);
  free(cox);
  free(coy);
  N_VFree(data->rates);

  free(data);

}


/***********************************************************************
 *
 * Function: SetInitialProfiles
 *
 * Set initial conditions in cc
 *
 **********************************************************************/

static void SetInitialProfiles(N_Vector cc, N_Vector sc)
{
  int i, jx, jy;
  real *ct1, *st1, *ct2, *st2;
  real  ctemp[NUM_SPECIES], stemp[NUM_SPECIES];

  /* Initialize temporary arrays ctemp and stemp to be used
          in the loading process */

  for(i=0;i<NUM_SPECIES;i++)
    if(i<NUM_SPECIES/2){
      ctemp[i]=RCONST(1.16347);
      stemp[i]=ONE;}
    else {
      ctemp[i]=RCONST(34903.1);
      stemp[i]=RCONST(0.00001);}

  /* Load initial profiles into cc and sc vector from temporary arrays */

  for (jy=0; jy < MYSUB; jy++) {
    for (jx=0; jx < MXSUB; jx++) {
      ct1 = IJ_Vptr(cc,jx,jy);
      ct2 = ctemp;
      st1 = IJ_Vptr(sc,jx,jy);
      st2 = stemp;
      for(i=0;i<NUM_SPECIES;i++){
	*ct1++=*ct2++;
	*st1++=*st2++;
      }
    }
  }

}  /* end SetInitialProfiles */


/***********************************************************************
 *
 * Function: PrintOutput
 *
 * Print sample of current cc value
 *
 **********************************************************************/
static void PrintOutput(integer my_pe, MPI_Comm comm, N_Vector cc)
{
  int   is;
  int   i0;
  int   npelast;
  real  *ct;
  real  tempc[NUM_SPECIES];

#ifdef PARALLEL
  MPI_Status status;
#endif

  npelast = NPEX*NPEY - 1;  /* Number of last processor */

  ct = N_VDATA(cc);


  /*
   * Send the cc values (for all species) at the top right mesh point
   * to PE 0
   *
  */

  i0 = NUM_SPECIES*(MXSUB*MYSUB-1);

#ifdef PARALLEL
  if (my_pe == npelast) {
    if(npelast!=0)
       MPI_Send(&ct[i0],NUM_SPECIES,PVEC_REAL_MPI_TYPE,0,0,comm);
    else  /* single processor case */
       for (is=0; is<NUM_SPECIES; is++) tempc[is] = ct[i0+is];
  }
#else
  if (my_pe == npelast) {
    for (is=0; is<NUM_SPECIES; is++) tempc[is] = ct[i0+is];
  }
#endif

  /*
   * On PE 0, receive the cc values at top right, then print performance
   * data and sampled solution values
   *
  */

  if(my_pe == 0) {

#ifdef PARALLEL
    if(npelast != 0) {
      MPI_Recv (&tempc[0], NUM_SPECIES, PVEC_REAL_MPI_TYPE,
         npelast, 0, comm, &status);
    }
#endif

    printf("\n");
    printf("At bottom left::\n");
    for(is=0; is<NUM_SPECIES; is++) {
      if ((is%4) == 0) printf("\n");
      printf(" %16.8e",ct[is]);
    }

    printf("\n\n");
    printf("At top right:\n");
    for(is=0; is<NUM_SPECIES; is++) {
      if ((is%4) == 0) printf("\n");
      printf(" %16.8e",tempc[is]);
    }
    printf("\n\n");
  }
}


/***********************************************************************
 *
 * Function: PrintSensOut
 *
 * Print sample of current ww value
 *
 **********************************************************************/
static void PrintSensOut(integer my_pe, MPI_Comm comm, N_Vector ww)
{
  int   is;
  int   i0;
  int   npelast;
  real  *wt;
  real  tempw[NUM_SPECIES];

#ifdef PARALLEL
  MPI_Status status;
#endif

  npelast = NPEX*NPEY - 1;  /* Number of last processor */

  wt = N_VDATA(ww);


  /*
   * Send the ww values (for all species) at the top right mesh point
   * to PE 0
   *
  */

  i0 = NUM_SPECIES*(MXSUB*MYSUB-1);

#ifdef PARALLEL
  if (my_pe == npelast) {
    if(npelast!=0)
       MPI_Send(&wt[i0],NUM_SPECIES,PVEC_REAL_MPI_TYPE,0,0,comm);
    else  /* single processor case */
       for (is=0; is<NUM_SPECIES; is++) tempw[is] = wt[i0+is];
  }
#else
  if (my_pe == npelast) {
    for (is=0; is<NUM_SPECIES; is++) tempw[is] = wt[i0+is];
  }
#endif

  /*
   * On PE 0, receive the ww values at top right, then print performance
   * data and sampled sensitivity values
   *
  */

  if(my_pe == 0) {

#ifdef PARALLEL
    if(npelast != 0) {
      MPI_Recv (&tempw[0], NUM_SPECIES, PVEC_REAL_MPI_TYPE,
         npelast, 0, comm, &status);
    }
#endif

    printf("\n");
    printf("At bottom left::\n");
    for(is=0; is<NUM_SPECIES; is++) {
      if ((is%6)*6 == is) printf("\n");
      printf(" %g",wt[is]);
    }

    printf("\n\n");
    printf("At top right:\n");
    for(is=0; is<NUM_SPECIES; is++) {
      if ((is%6)*6 == is) printf("\n");
      printf(" %g",tempw[is]);
    }
    printf("\n\n");
  }
}


/***********************************************************************
 *
 * Function: PrintFinalStats
 *
 * Print final statistics contained in iopt
 *
 **********************************************************************/

static void PrintFinalStats(long int *iopt)
{
  printf("\nFinal Statistics.. \n\n");
  printf("nni    = %5ld    nli   = %5ld\n", iopt[NNI], iopt[SPGMR_NLI]);
  printf("nfe    = %5ld    npe   = %5ld\n", iopt[NFE], iopt[SPGMR_NPE]);
  printf("nps    = %5ld    ncfl  = %5ld\n", iopt[SPGMR_NPS], iopt[SPGMR_NCFL]);
}


/***********************************************************************
 *
 * Function: fcalcprpr
 *
 * System function for predator - prey system  calculation part
 *
 **********************************************************************/

static void fcalcprpr(integer Neq, N_Vector cc, N_Vector fval,
   void *f_data)
{
  real xx, yy, *cxy, *rxy, *fxy;

  real    dcydi;     /* y down side index  */
  real    dcyui;     /* y up side index    */
  real    dcxli;     /* x left side index  */
  real    dcxri;     /* x right side index */

  integer isubx=0;   /* Processor offset in x direction */
  integer isuby=0;   /* Processor offset in x direction */
  integer nsmxsub;
  integer nsmxsub2;

  integer i, j, is, ly;
  real *cext, *cdata;
  integer shifty;
  integer offsetc, offsetce, offsetcl, offsetcr, offsetcd, offsetcu;
  UserData data;


  data  = (UserData) f_data;
  cdata = N_VDATA(cc);


  /* Get subgrid indices, data sizes, extended work array cext */

  isubx   = data->isubx;   isuby = data->isuby;
  nsmxsub = data->nsmxsub; nsmxsub2 = data->nsmxsub2;
  cext    = data->cext;

  /* Copy local segment of cc vector into the working extended array cext */

  offsetc = 0;
  offsetce = nsmxsub2 + NUM_SPECIES;
  for (ly = 0; ly < MYSUB; ly++) {
    for (i = 0; i < nsmxsub; i++) cext[offsetce+i] = cdata[offsetc+i];
    offsetc = offsetc + nsmxsub;
    offsetce = offsetce + nsmxsub2;
  }

  /* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of cc to cext */

  /* If isuby = 0, copy x-line 2 of cc to cext */
  if (isuby == 0) {
    for (i = 0; i < nsmxsub; i++) cext[NUM_SPECIES+i] = cdata[nsmxsub+i];
  }

  /* If isuby = NPEY-1, copy x-line MYSUB-1 of cc to cext */
  if (isuby == NPEY-1) {
    offsetc = (MYSUB-2)*nsmxsub;
    offsetce = (MYSUB+1)*nsmxsub2 + NUM_SPECIES;
    for (i = 0; i < nsmxsub; i++) cext[offsetce+i] = cdata[offsetc+i];
  }

  /* If isubx = 0, copy y-line 2 of cc to cext */
  if (isubx == 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetc = ly*nsmxsub + NUM_SPECIES;
      offsetce = (ly+1)*nsmxsub2;
      for (i = 0; i < NUM_SPECIES; i++) cext[offsetce+i] = cdata[offsetc+i];
    }
  }

  /* If isubx = NPEX-1, copy y-line MXSUB-1 of cc to cext */
  if (isubx == NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetc = (ly+1)*nsmxsub - 2*NUM_SPECIES;
      offsetce = (ly+2)*nsmxsub2 - NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++) cext[offsetce+i] = cdata[offsetc+i];
    }
  }


  /* loop over all grid points, evaluating for each species at each */

  shifty = (MXSUB+2)*NUM_SPECIES;

  for(j=0; j<MYSUB; j++) {
    yy = (data->dy) *(j  + isuby * MYSUB);

    for(i=0; i<MXSUB; i++){

      xx = (data->dx) * ( i + isubx * MXSUB);
      cxy = IJ_Vptr(cc,i,j);
      rxy = IJ_Vptr(data->rates,i,j);
      fxy = IJ_Vptr(fval,i,j);

      WebRate(xx, yy, cxy, rxy, f_data);

      offsetc = (i+1)*NUM_SPECIES + (j+1)*NSMXSUB2;
      offsetcd = offsetc - shifty;
      offsetcu = offsetc + shifty;
      offsetcl = offsetc - NUM_SPECIES;
      offsetcr = offsetc + NUM_SPECIES;


      for(is=0; is<NUM_SPECIES; is++){

	/* differencing in x */

	dcydi = cext[offsetc+is]  - cext[offsetcd+is];
        dcyui = cext[offsetcu+is] - cext[offsetc+is];

	/* differencing in y */

	dcxli = cext[offsetc+is]  - cext[offsetcl+is];
	dcxri = cext[offsetcr+is] - cext[offsetc+is];

	/* compute the value at xx , yy */

	fxy[is] = (coy)[is] * (dcyui - dcydi) +
	           (cox)[is] * (dcxri - dcxli) + rxy[is];

      } /* end is loop */

    } /* end of i or x  loop */

  } /* end of j or y loop */

}  /* end of routine fcalcprpr */


/***********************************************************************
 *                                                                     *
 ***************** Functions Called by the KINSol Solver ***************
 *                                                                     *
 **********************************************************************/


/***********************************************************************
 *
 * Function: funcprpr
 *
 * System function routine. Evaluate preditor-prey function.  First
 * call ccomm to do communication of subgrid boundary data into cext.
 * Then calculate funcprpr(cc) by a call to fcalcprpr
 *
 **********************************************************************/

static void funcprpr(integer Neq, N_Vector cc, N_Vector fval, void *f_data)
{
#ifdef PARALLEL
  real *cdata;
#endif
  UserData data;

  data = (UserData) f_data;

#ifdef PARALLEL
  /* Call ccomm to do inter-processor communicaiton */

  cdata = N_VDATA(cc);
  ccomm (Neq, cdata, data);
#endif

  /* Call fcalc to calculate all right-hand sides */

  fcalcprpr (Neq, cc, fval, data);

}


/***********************************************************************
 *
 * Function: Precondbd
 *
 * Preconditioner setup routine. Generate and preprocess P.
 *
 **********************************************************************/

static int Precondbd(integer Neq, N_Vector cc, N_Vector cscale,
		   N_Vector fval, N_Vector fscale,
		   N_Vector vtem, N_Vector vtemp1, SysFn func, real uround,
		   long int *nfePtr, void *P_data)
{
  real r, r0, sqruround;
  real xx, yy, *cxy, *scxy, cctemp, **Pxy, *ratesxy, *Pxycol;
  real fac, perturb_rates[NUM_SPECIES];


  integer i, j, jx, jy, ret;


  UserData data;

  data = (UserData)P_data;

  sqruround = data->sqruround;
  fac = N_VWL2Norm(fval, fscale);
  r0 = THOUSAND * uround * fac * Neq;

  if(r0 == ZERO) r0 = ONE;


  for(jy=0; jy<MYSUB; jy++){

    yy =data->dy *(jy + data->isuby * MYSUB);

    for(jx=0; jx<MXSUB; jx++){

      xx = data->dx * (jx + data->isubx * MXSUB);
      Pxy = (data->P)[jx][jy];
      cxy = IJ_Vptr(cc,jx,jy);
      scxy= IJ_Vptr(cscale,jx,jy);
      ratesxy = IJ_Vptr((data->rates),jx,jy);

      for(j=0; j<NUM_SPECIES; j++){

	cctemp=cxy[j];  /* save the j,jx,jy element of cc */
	r=MAX(sqruround * ABS(cctemp),r0/scxy[j]);
	cxy[j] += r; /* perturb the j,jx,jy element of cc */
	fac = ONE/r;

	WebRate(xx, yy, cxy, perturb_rates,data);

	Pxycol = Pxy[j];

	for(i=0; i<NUM_SPECIES; i++) {
	  Pxycol[i]=(perturb_rates[i]-ratesxy[i]) * fac;
	}

	/* restore j,jx,jy element of cc */
	cxy[j] = cctemp;

      } /* end of j loop */


      /*  lu decomposition of each block */

      ret = gefa(Pxy, NUM_SPECIES, (data->pivot)[jx][jy]);


      if(ret!=0)return(1);

    } /* end jx loop */

  } /* end jy loop */
  return(0);

}  /* end of routine Precondbd */


/***********************************************************************
 *
 * Function: PSolvebd
 *
 * Preconditioner solve routine
 *
 **********************************************************************/

static int PSolvebd(integer Neq, N_Vector cc, N_Vector cscale,
		  N_Vector fval, N_Vector fscale, N_Vector vv, N_Vector ftem,
		  SysFn func, real uround,
		  long int *nfePtr, void *P_data)
{
 real **Pxy, *vxy;
 integer *pivot, jx, jy;
 UserData data;

 data = (UserData)P_data;

 for(  jx=0; jx<MXSUB; jx++) {
   for(jy=0; jy<MYSUB; jy++){
     /* for a given jx,jy block, do the inversion process */
     /* vvxy is the address of the portion of the vector to which the
	inversion process is applied, and Pxy is the first address for the
	jx,jy block of P */
     pivot=(data->pivot)[jx][jy];
     Pxy = (data->P)[jx][jy];
     vxy = IJ_Vptr(vv,jx,jy);
     gesl(Pxy, NUM_SPECIES, pivot, vxy);

   } /* end of jy loop */

 } /* end of jx loop */

 return(0);

} /*  end of PSolvebd  */


/***********************************************************************
 *
 * Function: WebRate
 *
 * The routine evaluates the predator-prey birth/death function for
 * all species in a single grid cell. Variation of sensitivity
 * analysis parameters has been included for the species interaction
 * coefficients:
 *
 *    SP_SELFSELF    A species interacting with itself
 *
 *    SP_PREYPRED    Effects on prey of interacting with predators
 *
 *    SP_PREDPREY    Effects on predators of interacting with prey
 *
 **********************************************************************/

static void WebRate(real xx, real yy, real *cxy, real *ratesxy, void *f_data)
{
  integer i, ip, npr;
  real fac;
  UserData data;
  real *sens_p;

  real pintr, pself, ptmp[NUM_SPECIES];

  data = (UserData) f_data;

  npr = data->np;  /*The number of predators or prey separately  */
  sens_p = data->sens_p;      /* Sensitivity analysis parameters */


   if ( data->sens_active == FALSE ) {

      /* Loop over the base species without including */
      /* a sensitivity parameter variation            */

      for(i=0;i<NUM_SPECIES;i++)
        ratesxy[i]= DotProd(NUM_SPECIES, cxy, acoef[i]);

   } else {


      /* Loop over the base species with sensitivity */
      /* parameter variation included.               */

      pintr = ( SP_PREYPRED < NPSENS ) ? sens_p[SP_PREYPRED] : ONE;
      pself = ( SP_SELFSELF < NPSENS ) ? sens_p[SP_SELFSELF] : ONE;

      for ( i=0; i<npr; i++) {            /* Base species is prey */

         for ( ip=0; ip<npr; ip++) {
            ptmp[ip] = ONE;
            ptmp[ip+npr] = pintr;
         }
         ptmp[i] = pself;
         for ( ip=0; ip<NUM_SPECIES; ip++) ptmp[ip] *= acoef[i][ip];

         ratesxy[i]= DotProd(NUM_SPECIES, cxy, ptmp);
      }


      pintr = (SP_PREDPREY < NPSENS) ? sens_p[SP_PREDPREY] : ONE;

      for ( i=npr; i<NUM_SPECIES; i++) {  /* Base species is predator */

         for ( ip=0; ip<npr; ip++) {
            ptmp[ip] = pintr;
            ptmp[ip+npr] = ONE;
         }
         ptmp[i] = pself;
         for ( ip=0; ip<NUM_SPECIES; ip++) ptmp[ip] *= acoef[i][ip];

         ratesxy[i]= DotProd(NUM_SPECIES, cxy, ptmp);
      }

   }

   fac = ONE + ALPHA * xx * yy;

   for(i=0; i<NUM_SPECIES; i++) {
      ratesxy[i] = cxy[i] * ( bcoef[i] * fac + ratesxy[i] );
   }

} /* end WebRate */


/***********************************************************************
 *
 * Function: DotProd
 *
 **********************************************************************/

static real DotProd(integer size, real *x1, real *x2)
{
  integer i;
  real *xx1, *xx2, temp = ZERO;

  xx1 = x1; xx2 = x2;
  for(i=0; i<size; i++) temp += *xx1++ * *xx2++;
  return(temp);

}


#ifdef PARALLEL
/***********************************************************************
 *                                                                     *
 **** Communication routines for parallel solution from here to EOF ****
 *                                                                     *
 **********************************************************************/


/***********************************************************************
 *
 * Function: BSend
 *
 * Routine to send boundary data to neighboring PEs
 *
 **********************************************************************/

static void BSend(MPI_Comm comm, integer my_pe, integer isubx, integer isuby,
                  integer dsizex, integer dsizey, real *cdata)
{
  int i, ly;
  integer offsetc, offsetbuf;
  real bufleft[NUM_SPECIES*MYSUB], bufright[NUM_SPECIES*MYSUB];

  /* If isuby > 0, send data from bottom x-line of u */

  if (isuby != 0)
    MPI_Send(&cdata[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If isuby < NPEY-1, send data from top x-line of u */

  if (isuby != NPEY-1) {
    offsetc = (MYSUB-1)*dsizex;
    MPI_Send(&cdata[offsetc], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If isubx > 0, send data from left y-line of u (via bufleft) */

  if (isubx != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetc = ly*dsizex;
      for (i = 0; i < NUM_SPECIES; i++)
        bufleft[offsetbuf+i] = cdata[offsetc+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);
  }

  /* If isubx < NPEX-1, send data from right y-line of u (via bufright) */

  if (isubx != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetc = offsetbuf*MXSUB + (MXSUB-1)*NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++)
        bufright[offsetbuf+i] = cdata[offsetc+i];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe+1, 0, comm);
  }

}


/***********************************************************************
 *
 * Function: BRecvPost
 *
 * Routine to start receiving boundary data from neighboring PEs.
 *
 * 1) buffer should be able to hold 2*NUM_SPECIES*MYSUB real entries,
 *    should be passed to both the BRecvPost and BRecvWait functions,
 *    and should not be manipulated between the two calls.
 *
 * 2) request should have 4 entries, and should be passed in both
 *    calls also.
 *
 *
 **********************************************************************/

static void BRecvPost(MPI_Comm comm, MPI_Request request[], integer my_pe,
		      integer isubx, integer isuby,
		      integer dsizex, integer dsizey,
		      real *cext, real *buffer)
{
  integer offsetce;
  /* Have bufleft and bufright use the same buffer */
  real *bufleft = buffer, *bufright = buffer+NUM_SPECIES*MYSUB;

  /* If isuby > 0, receive data for bottom x-line of cext */
  if (isuby != 0)
    MPI_Irecv(&cext[NUM_SPECIES], dsizex, PVEC_REAL_MPI_TYPE,
    					 my_pe-NPEX, 0, comm, &request[0]);

  /* If isuby < NPEY-1, receive data for top x-line of cext */
  if (isuby != NPEY-1) {
    offsetce = NUM_SPECIES*(1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&cext[offsetce], dsizex, PVEC_REAL_MPI_TYPE,
                                         my_pe+NPEX, 0, comm, &request[1]);
  }

  /* If isubx > 0, receive data for left y-line of cext (via bufleft) */
  if (isubx != 0) {
    MPI_Irecv(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         my_pe-1, 0, comm, &request[2]);
  }

  /* If isubx < NPEX-1, receive data for right y-line of cext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Irecv(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         my_pe+1, 0, comm, &request[3]);
  }

}


/***********************************************************************
 *
 * Function: BRecvWait
 *
 * Routine to finish receiving boundary data from neighboring PE
 *
 * 1) buffer should be able to hold 2*NUM_SPECIES*MYSUB real entries,
 *    should be passed to both the BRecvPost and BRecvWait functions,
 *    and should not be manipulated between the two calls.
 *
 * 2) request should have 4 entries, and should be passed in both
 *    calls also.
 *
 * **********************************************************************/

static void BRecvWait(MPI_Request request[], integer isubx, integer isuby,
		      integer dsizex, real *cext, real *buffer)
{
  int i, ly;
  integer dsizex2, offsetce, offsetbuf;
  real *bufleft = buffer, *bufright = buffer+NUM_SPECIES*MYSUB;
  MPI_Status status;

  dsizex2 = dsizex + 2*NUM_SPECIES;

  /* If isuby > 0, receive data for bottom x-line of cext */
  if (isuby != 0)
    MPI_Wait(&request[0],&status);

  /* If isuby < NPEY-1, receive data for top x-line of cext */
  if (isuby != NPEY-1)
    MPI_Wait(&request[1],&status);

  /* If isubx > 0, receive data for left y-line of cext (via bufleft) */
  if (isubx != 0) {
    MPI_Wait(&request[2],&status);

    /* Copy the buffer to cext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetce = (ly+1)*dsizex2;
      for (i = 0; i < NUM_SPECIES; i++)
        cext[offsetce+i] = bufleft[offsetbuf+i];
    }
  }

  /* If isubx < NPEX-1, receive data for right y-line of cext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Wait(&request[3],&status);

    /* Copy the buffer to cext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetce = (ly+2)*dsizex2 - NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++)
	cext[offsetce+i] = bufright[offsetbuf+i];
    }
  }

}


/***********************************************************************
 *
 * Function: ccomm
 *
 * This routine performs all communication between processors of data
 * needed to calculate f
 *
 **********************************************************************/

static void ccomm(integer Neq,real *cdata, UserData data)
{

  real *cext, buffer[2*NUM_SPECIES*MYSUB];
  MPI_Comm comm;
  integer my_pe, isubx, isuby, nsmxsub, nsmysub;
  MPI_Request request[4];


  /* Get comm, my_pe, subgrid indices, data sizes, extended array cext */

  comm = data->comm;  my_pe = data->my_pe;
  isubx = data->isubx;   isuby = data->isuby;
  nsmxsub = data->nsmxsub;
  nsmysub = NUM_SPECIES*MYSUB;
  cext = data->cext;

  /* Start receiving boundary data from neighboring PEs */

  BRecvPost(comm, request, my_pe, isubx, isuby, nsmxsub, nsmysub, cext, buffer);

  /* Send data from boundary of local grid to neighboring PEs */

  BSend(comm, my_pe, isubx, isuby, nsmxsub, nsmysub, cdata);

  /* Finish receiving boundary data from neighboring PEs */

  BRecvWait(request, isubx, isuby, nsmxsub, cext, buffer);

}

/********** End of Parallel only communication routines ***************/

#endif
