/************************************************************************
 *                                                                      *
 * File       : cvfkx.c                                                 *
 * Programmers: Scott D. Cohen and Alan C. Hindmarsh and                *
 *              Radu Serban @ LLNL                                      *
 * Version of : 30 March 2003                                           *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * An ODE system is generated from the following 2-species diurnal      *
 * kinetics advection-diffusion PDE system in 2 space dimensions:       *
 *                                                                      *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)    *
 *                 + Ri(c1,c2,t)      for i = 1,2,   where              *
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,       *
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,                    *
 *   Kv(z) = Kv0*exp(z/5) ,                                             *
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)        *
 * vary diurnally.   The problem is posed on the square                 *
 *   0 <= x <= 20,    30 <= z <= 50   (all in km),                      *
 * with homogeneous Neumann boundary conditions, and for time t in      *
 *   0 <= t <= 86400 sec (1 day).                                       *
 * The PDE system is treated by central differences on a uniform        *
 * 10 x 10 mesh, with simple polynomial initial profiles.               *
 * The problem is solved with CVODES, with the BDF/GMRES method (i.e.   *
 * using the CVSPGMR linear solver) and the block-diagonal part of the  *
 * Newton matrix as a left preconditioner. A copy of the block-diagonal *
 * part of the Jacobian is saved and conditionally reused within the    *
 * Precond routine.                                                     *
 *                                                                      *
 * Optionally, CVODES can compute sensitivities with respect to the     *
 * problem parameters q1 and q2.                                        *
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and       *
 * STAGGERED1) can be used and sensitivities may be included in the     *
 * error test or not (error control set on FULL or PARTIAL,             *
 * respectively).                                                       *
 *                                                                      *
 * Execution:                                                           *
 *                                                                      *
 * If no sensitivities are desired:                                     *
 *    % cvskx -nosensi                                                  *
 * If sensitivities are to be computed:                                 *
 *    % cvskx -sensi sensi_meth err_con                                 *
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of    *
 * {full, partial}.                                                     *
 *                                                                      *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sundialstypes.h"     /* definitions of realtype, integertype            */
#include "cvodes.h"            /* main CVODES header file                         */
#include "iterativ.h"          /* contains the enum for types of preconditioning  */
#include "cvsspgmr.h"          /* use CVSPGMR linear solver each internal step    */
#include "smalldense.h"        /* use generic DENSE solver for preconditioning    */
#include "nvector_serial.h"    /* definitions of type N_Vector, macro NV_DATA_S   */
#include "sundialsmath.h"      /* contains SQR macro                              */

/* Problem Constants */

#define NUM_SPECIES  2             /* number of species */
#define C1_SCALE     1.0e6         /* coefficients in initial profiles */
#define C2_SCALE     1.0e12

#define T0           0.0           /* initial time */
#define NOUT         12            /* number of output times */
#define TWOHR        7200.0        /* number of seconds in two hours  */
#define HALFDAY      4.32e4        /* number of seconds in a half day */
#define PI       3.1415926535898   /* pi */ 

#define XMIN          0.0          /* grid boundaries in x  */
#define XMAX         20.0           
#define ZMIN         30.0          /* grid boundaries in z  */
#define ZMAX         50.0
#define XMID         10.0          /* grid midpoints in x,z */          
#define ZMID         40.0

#define MX           15             /* MX = number of x mesh points */
#define MZ           15             /* MZ = number of z mesh points */
#define NSMX         NUM_SPECIES*MX /* NSMX = NUM_SPECIES*MX */
#define MM           (MX*MZ)        /* MM = MX*MZ */

/* CVodeMalloc Constants */
#define RTOL    1.0e-5            /* scalar relative tolerance */
#define FLOOR   100.0             /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */
#define NEQ     (NUM_SPECIES*MM)  /* NEQ = number of equations */

/* Sensitivity Constants */
#define NP    8
#define NS    2

#define ZERO  RCONST(0.0)

/* User-defined vector and matrix accessor macros: IJKth, IJth */

/* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MZ-1. The vdata array is obtained via
   the macro call vdata = NV_DATA_S(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])


/* Type : UserData 
   contains preconditioner blocks, pivot arrays, 
   problem parameters, and problem constants     */

typedef struct {
  realtype *p;
  realtype **P[MX][MZ], **Jbd[MX][MZ];
  integertype *pivot[MX][MZ];
  realtype q4, om, dx, dz, hdco, haco, vdco;
} *UserData;


/* Private Helper Functions */

static void WrongArgs(char *argv[]);
static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector y, realtype dx, realtype dz);
static void PrintOutput(long int iopt[], realtype ropt[], realtype t, N_Vector y);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(booleantype sensi, int sensi_meth, int err_con, long int iopt[]);

/* Functions Called by the CVODES Solver */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

static int Precond(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, N_Vector ewt, realtype h,
                   realtype uround, long int *nfePtr, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int PSolve(realtype tn, N_Vector y, N_Vector fy, N_Vector vtemp,
                  realtype gamma, N_Vector ewt, realtype delta, long int *nfePtr,
                  N_Vector r, int lr, void *P_data, N_Vector z);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  M_Env machEnv;  
  realtype abstol, reltol, t, tout, ropt[OPT_SIZE];
  long int iopt[OPT_SIZE];
  N_Vector y;
  UserData data;
  void *cvode_mem;
  int iout, flag, i;

  realtype *pbar, rhomax;
  integertype is, *plist;
  N_Vector *uS;
  booleantype sensi;
  int sensi_meth, err_con, ifS;

  /* Process arguments */

  if (argc < 2)
    WrongArgs(argv);

  if (strcmp(argv[1],"-nosensi") == 0)
    sensi = FALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    sensi = TRUE;
  else
    WrongArgs(argv);

  if (sensi) {

    if (argc != 4)
      WrongArgs(argv);

    if (strcmp(argv[2],"sim") == 0)
      sensi_meth = SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      sensi_meth = STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      sensi_meth = STAGGERED1;
    else 
      WrongArgs(argv);

    if (strcmp(argv[3],"full") == 0)
      err_con = FULL;
    else if (strcmp(argv[3],"partial") == 0)
      err_con = PARTIAL;
    else
      WrongArgs(argv);

  }

  /* Initialize serial machine environment */
  machEnv = M_EnvInit_Serial(NEQ);

  /* PROBLEM PARAMETERS */
  data = AllocUserData();
  InitUserData(data);

  /* INITIAL STATES */
  y = N_VNew(machEnv);
  SetInitialProfiles(y, data->dx, data->dz);
  
  /* TOLERANCES */
  abstol=ATOL; 
  reltol=RTOL;

  /* OPTIONAL INPUT */
  for (i = 0; i < OPT_SIZE; i++) {
    iopt[i] = 0;
    ropt[i] = 0.0;
  }
  iopt[MXSTEP] = 2000;

  /* CVODE_MALLOC */
  cvode_mem = CVodeMalloc(f, T0, y, BDF, NEWTON, SS, &reltol,
                          &abstol, data, NULL, TRUE, iopt, ropt, machEnv);
  if (cvode_mem == NULL) { 
    printf("CVodeMalloc failed."); 
    return(1); 
  }

  /* CVSPGMR */
  flag = CVSpgmr(cvode_mem, LEFT, MODIFIED_GS, 0, 0.0, Precond, PSolve, data, NULL, NULL);
  if (flag != SUCCESS) { printf("CVSpgmr failed.\n"); return(1); }

  /* SENSITIVTY */
  if(sensi) {
    pbar = (realtype *) malloc(NP*sizeof(realtype));
    for(is=0; is<NP; is++) pbar[is] = data->p[is];
    plist = (integertype *) malloc(NS * sizeof(integertype));
    for(is=0; is<NS; is++) plist[is] = is+1;

    uS = N_VNew_S(NS, machEnv);
    for(is=0;is<NS;is++)
      N_VConst(ZERO,uS[is]);

    rhomax = ZERO;

    ifS = ALLSENS;
    if(sensi_meth==STAGGERED1) ifS = ONESENS;

    flag = CVodeSensMalloc(cvode_mem, NS, sensi_meth, data->p, pbar, plist,
                           ifS, NULL, err_con, rhomax, uS, NULL, NULL, NULL);
    if (flag != SUCCESS) {printf("CVodeSensMalloc failed, flag=%d\n",flag);return(1);}

  }

  /* In loop over output points, call CVode, print results, test for error */
  
  printf("\n2-species diurnal advection-diffusion problem\n\n");

  printf("========================================================================\n");
  printf("     T     Q       H      NST                    Bottom left  Top right \n");
  printf("========================================================================\n");

  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    if (flag != SUCCESS) { 
      printf("CVode failed, flag=%d.\n", flag); 
      break; 
    }
    PrintOutput(iopt, ropt, t, y);
    if (sensi) {
      flag = CVodeSensExtract(cvode_mem, t, uS);
      if (flag != SUCCESS) { 
        printf("CVodeSensExtract failed, flag=%d.\n", flag); 
        break; 
      }
      PrintOutputS(uS);
    }
    
    printf("------------------------------------------------------------------------\n");
    
  }

  /* Print final statistics */
  PrintFinalStats(sensi, sensi_meth, err_con, iopt);

  /* Free memory */
  N_VFree(y);
  if(sensi) N_VFree_S(NS, uS);
  FreeUserData(data);
  CVodeFree(cvode_mem);
  M_EnvFree_Serial(machEnv);

  return(0);
}


/*********************** Private Helper Functions ************************/

/* ======================================================================= */
/* Exit if arguments are incorrect */

static void WrongArgs(char *argv[])
{
  printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",argv[0]);
  printf("         sensi_meth = sim, stg, or stg1\n");
  printf("         err_con    = full or partial\n");
  
  exit(0);
}

/* ======================================================================= */
/* Allocate memory for data structure of type UserData */

static UserData AllocUserData(void)
{
  int jx, jz;
  UserData data;

  data = (UserData) malloc(sizeof *data);

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      (data->P)[jx][jz] = denalloc(NUM_SPECIES);
      (data->Jbd)[jx][jz] = denalloc(NUM_SPECIES);
      (data->pivot)[jx][jz] = denallocpiv(NUM_SPECIES);
    }
  }

  data->p = (realtype *) malloc(NP*sizeof(realtype));

  return(data);
}

/* ======================================================================= */
/* Load problem constants in data */

static void InitUserData(UserData data)
{
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Set problem parameters */
  Q1 = 1.63e-16; /* Q1  coefficients q1, q2, c3             */
  Q2 = 4.66e-16; /* Q2                                      */
  C3 = 3.7e16;   /* C3                                      */
  A3 = 22.62;    /* A3  coefficient in expression for q3(t) */
  A4 = 7.601;    /* A4  coefficient in expression for q4(t) */
  KH = 4.0e-6;   /* KH  horizontal diffusivity Kh           */ 
  VEL = 0.001;   /* VEL advection velocity V                */
  KV0 = 1.0e-8;  /* KV0 coefficient in Kv(z)                */  

  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/(MX-1);
  data->dz = (ZMAX-ZMIN)/(MZ-1);
  data->hdco = KH/SQR(data->dx);
  data->haco = VEL/(2.0*data->dx);
  data->vdco = (1.0/SQR(data->dz))*KV0;

  data->p[0] = Q1;
  data->p[1] = Q2;
  data->p[2] = C3;
  data->p[3] = A3;
  data->p[4] = A4;
  data->p[5] = KH;
  data->p[6] = VEL;
  data->p[7] = KV0;
}

/* ======================================================================= */
/* Free data memory */

static void FreeUserData(UserData data)
{
  int jx, jz;

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      denfree((data->P)[jx][jz]);
      denfree((data->Jbd)[jx][jz]);
      denfreepiv((data->pivot)[jx][jz]);
    }
  }

  free(data->p);

  free(data);
}

/* ======================================================================= */
/* Set initial conditions in y */

static void SetInitialProfiles(N_Vector y, realtype dx, realtype dz)
{
  int jx, jz;
  realtype x, z, cx, cz;
  realtype *ydata;

  /* Set pointer to data array in vector y. */

  ydata = NV_DATA_S(y);

  /* Load initial profiles of c1 and c2 into y vector */

  for (jz=0; jz < MZ; jz++) {
    z = ZMIN + jz*dz;
    cz = SQR(0.1*(z - ZMID));
    cz = 1.0 - cz + 0.5*SQR(cz);
    for (jx=0; jx < MX; jx++) {
      x = XMIN + jx*dx;
      cx = SQR(0.1*(x - XMID));
      cx = 1.0 - cx + 0.5*SQR(cx);
      IJKth(ydata,1,jx,jz) = C1_SCALE*cx*cz; 
      IJKth(ydata,2,jx,jz) = C2_SCALE*cx*cz;
    }
  }
}

/* ======================================================================= */
/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(long int iopt[], realtype ropt[], realtype t, N_Vector y)
{
  realtype *ydata;

  ydata = NV_DATA_S(y);
  
  printf("%8.3e %2ld  %8.3e %5ld\n", t,iopt[QU],ropt[HU],iopt[NST]);
  printf("                                Solution       ");
  printf("%12.4e %12.4e \n", IJKth(ydata,1,0,0), IJKth(ydata,1,MX-1,MZ-1)); 
  printf("                                               ");
  printf("%12.4e %12.4e \n", IJKth(ydata,2,0,0), IJKth(ydata,2,MX-1,MZ-1));
  
}

/* ======================================================================= */
/* Print sampled sensitivities */

static void PrintOutputS(N_Vector *uS)
{
  realtype *sdata;

  sdata = NV_DATA_S(uS[0]);

  printf("                                ----------------------------------------\n"); 
  printf("                                Sensitivity 1  ");
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
  printf("                                               ");
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));

  sdata = NV_DATA_S(uS[1]);

  printf("                                ----------------------------------------\n"); 
  printf("                                Sensitivity 2  ");
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
  printf("                                               ");
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));

}

/* ======================================================================= */
/* Print final statistics contained in iopt */

static void PrintFinalStats(booleantype sensi, int sensi_meth, int err_con, long int iopt[])
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


/***************** Functions Called by the CVODES Solver ******************/

/* ======================================================================= */
/* f routine. Compute f(t,y). */

static void f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, czdn, czup, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, zdn, zup;
  realtype q4coef, delz, verdco, hordco, horaco;
  realtype *ydata, *dydata;
  int jx, jz, idn, iup, ileft, iright;
  UserData data;
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  data = (UserData) f_data;
  ydata = NV_DATA_S(y);
  dydata = NV_DATA_S(ydot);

  /* Load problem coefficients and parameters */

  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];
  A3 = data->p[3];
  A4 = data->p[4];
  KH = data->p[5];
  VEL = data->p[6];
  KV0 = data->p[7];

  /* Set diurnal rate coefficients. */

  s = sin(data->om*t);
  if (s > 0.0) {
    q3 = exp(-A3/s);
    data->q4 = exp(-A4/s);
  } else {
    q3 = 0.0;
    data->q4 = 0.0;
  }

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  delz = data->dz;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Loop over all grid points. */

  for (jz=0; jz < MZ; jz++) {

    /* Set vertical diffusion coefficients at jz +- 1/2 */

    zdn = ZMIN + (jz - .5)*delz;
    zup = zdn + delz;
    czdn = verdco*exp(0.2*zdn);
    czup = verdco*exp(0.2*zup);
    idn = (jz == 0) ? 1 : -1;
    iup = (jz == MZ-1) ? -1 : 1;
    for (jx=0; jx < MX; jx++) {

      /* Extract c1 and c2, and set kinetic rate terms. */

      c1 = IJKth(ydata,1,jx,jz); 
      c2 = IJKth(ydata,2,jx,jz);
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms. */

      c1dn = IJKth(ydata,1,jx,jz+idn);
      c2dn = IJKth(ydata,2,jx,jz+idn);
      c1up = IJKth(ydata,1,jx,jz+iup);
      c2up = IJKth(ydata,2,jx,jz+iup);
      vertd1 = czup*(c1up - c1) - czdn*(c1 - c1dn);
      vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms. */

      ileft = (jx == 0) ? 1 : -1;
      iright =(jx == MX-1) ? -1 : 1;
      c1lt = IJKth(ydata,1,jx+ileft,jz); 
      c2lt = IJKth(ydata,2,jx+ileft,jz);
      c1rt = IJKth(ydata,1,jx+iright,jz);
      c2rt = IJKth(ydata,2,jx+iright,jz);
      hord1 = hordco*(c1rt - 2.0*c1 + c1lt);
      hord2 = hordco*(c2rt - 2.0*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);

      /* Load all terms into ydot. */

      IJKth(dydata, 1, jx, jz) = vertd1 + hord1 + horad1 + rkin1; 
      IJKth(dydata, 2, jx, jz) = vertd2 + hord2 + horad2 + rkin2;
    }
  }
}

/* ======================================================================= */
/* Preconditioner setup routine. Generate and preprocess P. */

static int Precond(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, N_Vector ewt, realtype h,
                   realtype uround, long int *nfePtr, void *P_data,
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  realtype c1, c2, czdn, czup, diag, zdn, zup, q4coef, delz, verdco, hordco;
  realtype **(*P)[MZ], **(*Jbd)[MZ];
  integertype *(*pivot)[MZ], ier;
  int jx, jz;
  realtype *ydata, **a, **j;
  UserData data;
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Make local copies of pointers in P_data, and of pointer to y's data */
  data = (UserData) P_data;
  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  ydata = NV_DATA_S(y);

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

    for (jz=0; jz < MZ; jz++)
      for (jx=0; jx < MX; jx++)
        dencopy(Jbd[jx][jz], P[jx][jz], NUM_SPECIES);

  *jcurPtr = FALSE;

  }

  else {
  /* jok = FALSE: Generate Jbd from scratch and copy to P */

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  delz = data->dz;
  verdco = data->vdco;
  hordco  = data->hdco;

  /* Compute 2x2 diagonal Jacobian blocks (using q4 values 
     computed on the last f call).  Load into P. */

    for (jz=0; jz < MZ; jz++) {
      zdn = ZMIN + (jz - .5)*delz;
      zup = zdn + delz;
      czdn = verdco*exp(0.2*zdn);
      czup = verdco*exp(0.2*zup);
      diag = -(czdn + czup + 2.0*hordco);
      for (jx=0; jx < MX; jx++) {
        c1 = IJKth(ydata,1,jx,jz);
        c2 = IJKth(ydata,2,jx,jz);
        j = Jbd[jx][jz];
        a = P[jx][jz];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + q4coef;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4coef) + diag;
        dencopy(j, a, NUM_SPECIES);
      }
    }

  *jcurPtr = TRUE;

  }

  /* Scale by -gamma */

    for (jz=0; jz < MZ; jz++)
      for (jx=0; jx < MX; jx++)
        denscale(-gamma, P[jx][jz], NUM_SPECIES);

  /* Add identity matrix and do LU decompositions on blocks in place. */

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      denaddI(P[jx][jz], NUM_SPECIES);
      ier = gefa(P[jx][jz], NUM_SPECIES, pivot[jx][jz]);
      if (ier != 0) return(1);
    }
  }

  return(0);
}


/* ======================================================================= */
/* Preconditioner solve routine */

static int PSolve(realtype tn, N_Vector y, N_Vector fy, N_Vector vtemp,
                  realtype gamma, N_Vector ewt, realtype delta, long int *nfePtr,
                  N_Vector r, int lr, void *P_data, N_Vector z)
{
  realtype **(*P)[MZ];
  integertype *(*pivot)[MZ];
  int jx, jz;
  realtype *zdata, *v;
  UserData data;

  /* Extract the P and pivot arrays from P_data. */

  data = (UserData) P_data;
  P = data->P;
  pivot = data->pivot;
  zdata = NV_DATA_S(z);

  N_VScale(1.0, r, z);

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      v = &(IJKth(zdata, 1, jx, jz));
      gesl(P[jx][jz], NUM_SPECIES, pivot[jx][jz], v);
    }
  }

  return(0);
}
 
