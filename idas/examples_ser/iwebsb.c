/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2004-07-22 23:01:56 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example program for IDA/IDAS: Food web, serial, band solve IDABAND.
 *
 * This example program for IDA/IDAS (serial version) uses IDABAND as the
 * linear solver, and IDACalcIC for initial condition calculation.
 *
 * The mathematical problem solved in this example is a DAE system
 * that arises from a system of partial differential equations after
 * spatial discretization. The PDE system is a food web population
 * model, with predator-prey interaction and diffusion on the unit
 * square in two dimensions. The dependent variable vector is:
 *
 *         1   2         ns
 *   c = (c , c ,  ..., c  ) , ns = 2 * np
 *
 * and the PDE's are as follows:
 *
 *     i             i      i
 *   dc /dt = d(i)*(c    + c  )  +  R (x,y,c)   (i = 1,...,np)
 *                   xx     yy       i
 *
 *              i      i
 *   0 = d(i)*(c    + c  )  +  R (x,y,c)   (i = np+1,...,ns)
 *              xx     yy       i
 *
 *   where the reaction terms R are:
 *
 *                   i             ns         j
 *   R  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *   a(i,i) = -AA   (all i)
 *   a(i,j) = -GG   (i <= np , j >  np)
 *   a(i,j) =  EE   (i >  np, j <= np)
 *   all other a(i,j) = 0
 *   b(i) = BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))   (i <= np)
 *   b(i) =-BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))   (i  > np)
 *   d(i) = DPREY   (i <= np)
 *   d(i) = DPRED   (i > np)
 *
 * Note: The above equations are written in 1-based indices,
 * whereas the code has 0-based indices, being written in C.
 *
 * The various scalar parameters required are set using '#define'
 * statements or directly in routine InitUserData. In this program,
 * np = 1, ns = 2. The boundary conditions are homogeneous Neumann:
 * normal derivative = 0.
 *
 * A polynomial in x and y is used to set the initial values of the
 * first np variables (the prey variables) at each x,y location,
 * while initial values for the remaining (predator) variables are
 * set to a flat value, which is corrected by IDACalcIC.
 *
 * The PDEs are discretized by central differencing on a MX by MY
 * mesh.
 *
 * The DAE system is solved by IDA using the IDABAND linear solver.
 * Output is printed at t = 0, .001, .01, .1, .4, .7, 1.
 * -----------------------------------------------------------------
 * References:
 * [1] Peter N. Brown and Alan C. Hindmarsh,
 *     Reduced Storage Matrix Methods in Stiff ODE systems, Journal
 *     of Applied Mathematics and Computation, Vol. 31 (May 1989),
 *     pp. 40-91.
 *
 * [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Using Krylov Methods in the Solution of Large-Scale
 *     Differential-Algebraic Systems, SIAM J. Sci. Comput., 15
 *     (1994), pp. 1467-1488.
 *
 * [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Consistent Initial Condition Calculation for Differential-
 *     Algebraic Systems, SIAM J. Sci. Comput., 19 (1998),
 *     pp. 1495-1512.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"   /* Definitions of realtype and booleantype       */
#include "idas.h"            /* Main header file                              */
#include "idaband.h"         /* Use IDABAND linear solver                     */
#include "nvector_serial.h"  /* Definitions of type N_Vector, macro NV_DATA_S */
#include "sundialsmath.h"    /* Contains RSqrt routine                        */
#include "smalldense.h"      /* Contains definitions for denalloc routine     */

/* Problem Constants. */

#define NPREY       1        /* Number of prey (= number of predators). */
#define NUM_SPECIES 2*NPREY

#define PI          3.1415926535898   /* pi */ 
#define FOURPI      (4.0*PI)          /* 4 pi */

#define MX          20                /* MX = number of x mesh points */
#define MY          20                /* MY = number of y mesh points */
#define NSMX        (NUM_SPECIES * MX)
#define NEQ         (NUM_SPECIES*MX*MY) /* Number of equations in system */
#define AA          RCONST(1.0)    /* Coefficient in above eqns. for a */
#define EE          RCONST(10000.) /* Coefficient in above eqns. for a */
#define GG          RCONST(0.5e-6) /* Coefficient in above eqns. for a */
#define BB          RCONST(1.0)    /* Coefficient in above eqns. for b */
#define DPREY       RCONST(1.0)    /* Coefficient in above eqns. for d */
#define DPRED       RCONST(0.05)   /* Coefficient in above eqns. for d */
#define ALPHA       RCONST(50.)    /* Coefficient alpha in above eqns. */
#define BETA        RCONST(1000.)  /* Coefficient beta in above eqns. */
#define AX          RCONST(1.0)    /* Total range of x variable */
#define AY          RCONST(1.0)    /* Total range of y variable */
#define RTOL        RCONST(1.e-5)  /*  rtol tolerance */
#define ATOL        RCONST(1.e-5)  /*  atol tolerance */
#define ZERO        RCONST(0.)     /* 0. */
#define ONE         RCONST(1.0)    /* 1. */
#define NOUT        6
#define TMULT       RCONST(10.0)   /* Multiplier for tout values */
#define TADD        RCONST(0.3)    /* Increment for tout values */

/* 
 * User-defined vector and accessor macro: IJ_Vptr.
 * IJ_Vptr is defined in order to express the underlying 3-D structure of the 
 * dependent variable vector from its underlying 1-D storage (an N_Vector).
 * IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
 * species index is = 0, x-index ix = i, and y-index jy = j.                
 */

#define IJ_Vptr(vv,i,j) (&NV_Ith_S(vv, i*NUM_SPECIES + j*NSMX))

/* Type: UserData.  Contains problem constants, etc. */

typedef struct {
  long int Neq, ns, np, mx, my;
  realtype dx, dy, **acoef;
  realtype cox[NUM_SPECIES], coy[NUM_SPECIES], bcoef[NUM_SPECIES];
  N_Vector rates;
} *UserData;

/* Prototypes for private Helper Functions. */

static UserData AllocUserData(long int SystemSize);
static void InitUserData(UserData webdata);
static void FreeUserData(UserData webdata);
static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
                               UserData webdata);
static void PrintOutput(void *mem, N_Vector cc, realtype time,
                        UserData webdata);
static void PrintFinalStats(void *mem);
static void Fweb(realtype tcalc, N_Vector cc, N_Vector crate, UserData webdata);
static void WebRates(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
                     UserData webdata);
static realtype dotprod(long int size, realtype *x1, realtype *x2);


/* Prototypes for functions called by the IDA Solver. */

static int resweb(realtype time, N_Vector cc, N_Vector cp, N_Vector resval, 
                  void *rdata);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

/*
 *----------------------------
 * Main program
 *----------------------------
 */

int main()
{ 
  void *mem;
  UserData webdata;
  N_Vector cc, cp, id;
  int iout, flag, retval, itol, itask;
  long int SystemSize, mu, ml;
  realtype rtol, atol, t0, tout, tret;

  mem = NULL;
  webdata = NULL;
  cc = cp = id = NULL;

  /* Set up and load user data block webdata. */

  SystemSize = NEQ;
  webdata = AllocUserData(SystemSize);
  if(check_flag((void *)webdata, "AllocUserData", 2)) return(1);

  InitUserData(webdata);

  /* Allocate N-vectors and initialize cc, cp, and id. */

  cc  = N_VNew_Serial(NEQ);
  if(check_flag((void *)cc, "N_VNew_Serial", 0)) return(1);

  cp  = N_VNew_Serial(NEQ);
  if(check_flag((void *)cp, "N_VNew_Serial", 0)) return(1);

  id  = N_VNew_Serial(NEQ);
  if(check_flag((void *)id, "N_VNew_Serial", 0)) return(1);
  
  SetInitialProfiles(cc, cp, id, webdata);
  
  /* Set remaining inputs to IDAMalloc. */
  
  t0 = ZERO;
  itol = SS; rtol = RTOL; atol = ATOL;

  /* Call IDACreate and IDAMalloc to initialize IDA. */
  
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);

  retval = IDASetRdata(mem, webdata);
  if(check_flag(&retval, "IDASetRdata", 1)) return(1);

  retval = IDASetId(mem, id);
  if(check_flag(&retval, "IDASetId", 1)) return(1);

  retval = IDAMalloc(mem, resweb, t0, cc, cp, itol, &rtol, &atol);
  if(check_flag(&retval, "IDAMalloc", 1)) return(1);

  /* Call IDABand to specify the IDA linear solver. */

  mu = ml = NSMX;
  retval = IDABand(mem, SystemSize, mu, ml);
  if(check_flag(&retval, "IDABand", 1)) return(1);

  /* Call IDACalcIC (with default options) to correct the initial values. */

  tout = 0.001;
  retval = IDACalcIC(mem, CALC_YA_YDP_INIT, tout);
  if(check_flag(&retval, "IDACalcIC", 1)) return(1);
  
  /* Print heading, basic parameters, and initial values. */

  printf("iwebsb: Predator-prey DAE serial example problem for IDA \n\n");
  printf("Number of species ns: %d", NUM_SPECIES);
  printf("     Mesh dimensions: %d x %d", MX, MY);
  printf("     System size: %ld\n",SystemSize);
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
  printf("Linear solver: IDABAND,  Band parameters mu = %ld, ml = %ld\n",mu,ml);
  printf("CalcIC called to correct initial predator concentrations \n\n");
  
  PrintOutput(mem, cc, ZERO, webdata);
  
  /* Loop over iout, call IDASolve (normal mode), print selected output. */
  
  itask = NORMAL;
  for (iout = 1; iout <= NOUT; iout++) {
    
    flag = IDASolve(mem, tout, &tret, cc, cp, itask);
    if(check_flag(&flag, "IDASolve", 1)) return(flag);
    
    PrintOutput(mem, cc, tret, webdata);
    
    if (iout < 3) tout *= TMULT; else tout += TADD;
    
  } /* End of iout loop. */
  
  /* Print final statistics and free memory. */  
  
  PrintFinalStats(mem);

  N_VDestroy(cc);
  N_VDestroy(cp);
  N_VDestroy(id);
  IDAFree(mem);
  FreeUserData(webdata);

  return(0);
}

/*
 *----------------------------
 * Private Helper Functions
 *----------------------------
 */

/* 
 * AllocUserData: Allocate memory for data structure of type UserData.   
 */

static UserData AllocUserData(long int SystemSize)
{
  UserData webdata;

  webdata = (UserData) malloc(sizeof *webdata);

  webdata->rates = N_VNew_Serial(SystemSize);

  webdata->acoef = denalloc(NUM_SPECIES);
 
  return(webdata);
}

/* Define lines for readability in later routines */ 

#define acoef  (webdata->acoef)
#define bcoef  (webdata->bcoef)
#define cox    (webdata->cox)
#define coy    (webdata->coy)

/*
 * InitUserData: Load problem constants in webdata (of type UserData).   
 */

static void InitUserData(UserData webdata)
{
  int i, j, np;
  realtype *a1,*a2, *a3, *a4, dx2, dy2;

  webdata->mx = MX;
  webdata->my = MY;
  webdata->ns = NUM_SPECIES;
  webdata->np = NPREY;
  webdata->dx = AX/(MX-1);
  webdata->dy = AY/(MY-1);
  webdata->Neq= NEQ;
  
  /* Set up the coefficients a and b, and others found in the equations. */
  np = webdata->np;
  dx2 = (webdata->dx)*(webdata->dx); dy2 = (webdata->dy)*(webdata->dy);
  
  for (i = 0; i < np; i++) {
    a1 = &(acoef[i][np]);
    a2 = &(acoef[i+np][0]);
    a3 = &(acoef[i][0]);
    a4 = &(acoef[i+np][np]);
    /*  Fill in the portion of acoef in the four quadrants, row by row. */
    for (j = 0; j < np; j++) {
      *a1++ =  -GG;
      *a2++ =   EE;
      *a3++ = ZERO;
      *a4++ = ZERO;
    }
    
    /* Reset the diagonal elements of acoef to -AA. */
    acoef[i][i] = -AA; acoef[i+np][i+np] = -AA;
    
    /* Set coefficients for b and diffusion terms. */
    bcoef[i] = BB; bcoef[i+np] = -BB;
    cox[i] = DPREY/dx2; cox[i+np] = DPRED/dx2;
    coy[i] = DPREY/dy2; coy[i+np] = DPRED/dy2;
  }
  
}


/* 
 * FreeUserData: Free webdata memory.
 */

static void FreeUserData(UserData webdata)
{
  denfree(acoef);
  N_VDestroy(webdata->rates);
  
  free(webdata);
}

/* 
 * SetInitialProfiles: Set initial conditions in cc, cp, and id.
 * A polynomial profile is used for the prey cc values, and a constant
 * (1.0e5) is loaded as the initial guess for the predator cc values.
 * The id values are set to 1 for the prey and 0 for the predators.
 * The prey cp values are set according to the given system, and
 * the predator cp values are set to zero.                               
 */

static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
                               UserData webdata)
{
  long int loc, yloc, is, jx, jy, np;
  realtype xx, yy, xyfactor, fac;
  realtype *ccv, *cpv, *idv;
  
  ccv = NV_DATA_S(cc);
  cpv = NV_DATA_S(cp);
  idv = NV_DATA_S(id);
  np = webdata->np;
  
  /* Loop over grid, load cc values and id values. */
  for (jy = 0; jy < MY; jy++) {
    yy = jy * webdata->dy;
    yloc = NSMX * jy;
    for (jx = 0; jx < MX; jx++) {
      xx = jx * webdata->dx;
      xyfactor = 16.*xx*(1.-xx)*yy*(1.-yy);
      xyfactor *= xyfactor;
      loc = yloc + NUM_SPECIES*jx;
      fac = ONE + ALPHA * xx * yy + BETA * sin(FOURPI*xx) * sin(FOURPI*yy);
      
      for (is = 0; is < NUM_SPECIES; is++) {
        if (is < np) {
          ccv[loc+is] = 10. + (realtype)(is+1) * xyfactor;
          idv[loc+is] = ONE;
        }
        else {
          ccv[loc+is] = 1.0e5;
          idv[loc+is] = ZERO;
        }
      }
    }
  }
  
  /* Set c' for the prey by calling the function Fweb. */
  Fweb(ZERO, cc, cp, webdata);
  
  /* Set c' for predators to 0. */
  for (jy = 0; jy < MY; jy++) {
    yloc = NSMX * jy;
    for (jx = 0; jx < MX; jx++) {
      loc = yloc + NUM_SPECIES * jx;
      for (is = np; is < NUM_SPECIES; is++) {
        cpv[loc+is] = ZERO;
      }
    }
  }
}

/* 
 * PrintOutput: Print output values at output time t = tt.
 * Selected run statistics are printed.  Then values of c1 and c2
 * are printed for the bottom left and top right grid points only.  
 * (NOTE: This routine is specific to the case NUM_SPECIES = 2.)         
 */

static void PrintOutput(void *mem, N_Vector cc, realtype tt,
                        UserData webdata)
{
  int jx, jy, kused, flag;
  long int nst, nni, nje, nre, nreB;
  realtype *ct, hused;

  flag = IDAGetLastOrder(mem, &kused);
  check_flag(&flag, "IDAGetLastOrder", 1);
  flag = IDAGetNumSteps(mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&flag, "IDAGetNumNonlinSolvIters", 1);
  flag = IDAGetNumResEvals(mem, &nre);
  check_flag(&flag, "IDAGetNumResEvals", 1);
  flag = IDAGetLastStep(mem, &hused);
  check_flag(&flag, "IDAGetLastStep", 1);
  flag = IDABandGetNumJacEvals(mem, &nje);
  check_flag(&flag, "IDABandGetNumJacEvals", 1);
  flag = IDABandGetNumResEvals(mem, &nreB);
  check_flag(&flag, "IDABandGetNumResEvals", 1);
  
  printf("\nTIME t = %e.     NST = %ld, k = %d, h = %e\n",
         tt, nst, kused, hused);
  printf("NRE = %ld, NRE_B = %ld, NNI = %ld, NJE = %ld\n", 
         nre, nreB, nni, nje);
  
  jx = 0;    jy = 0;    ct = IJ_Vptr(cc,jx,jy);
  printf("At bottom left:  c1, c2 = %e %e \n",   ct[0],ct[1]);
  
  jx = MX-1; jy = MY-1; ct = IJ_Vptr(cc,jx,jy);
  printf("At top right:    c1, c2 = %e %e \n\n", ct[0],ct[1]);
  
}

/* 
 * PrintFinalStats: Print final run data contained in iopt.              
 */

static void PrintFinalStats(void *mem)
{ 
  long int nst, nre, nreB, nni, nje, netf, ncfn;
  int flag;

  flag = IDAGetNumSteps(mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&flag, "IDAGetNumNonlinSolvIters", 1);
  flag = IDAGetNumResEvals(mem, &nre);
  check_flag(&flag, "IDAGetNumResEvals", 1);
  flag = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&flag, "IDAGetNumErrTestFails", 1);
  flag = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&flag, "IDAGetNumNonlinSolvConvFails", 1);
  flag = IDABandGetNumJacEvals(mem, &nje);
  check_flag(&flag, "IDABandGetNumJacEvals", 1);
  flag = IDABandGetNumResEvals(mem, &nreB);
  check_flag(&flag, "IDABandGetNumResEvals", 1);

  printf("\nFinal run statistics: \n\n");
  printf("NST  = %5ld     NRE  = %5ld \n", nst, nre+nreB);
  printf("NNI  = %5ld     NJE  = %5ld \n", nni, nje);
  printf("NETF = %5ld     NCFN = %5ld \n", netf, ncfn);
  
}

/* 
 * resweb: System residual function for predator-prey system.
 * This routine calls Fweb to get all the right-hand sides of the
 * equations, then loads the residual vector accordingly,
 * using cp in the case of prey species.                 
 */

static int resweb(realtype tt, N_Vector cc, N_Vector cp, 
                  N_Vector res,  void *rdata)
{
  long int jx, jy, is, yloc, loc, np;
  realtype *resv, *cpv;
  UserData webdata;
  
  webdata = (UserData)rdata;
  
  cpv = NV_DATA_S(cp);
  resv = NV_DATA_S(res);
  np = webdata->np;
  
  /* Call Fweb to set res to vector of right-hand sides. */
  Fweb(tt, cc, res, webdata);
  
  /* Loop over all grid points, setting residual values appropriately
     for differential or algebraic components.                        */
  
  for (jy = 0; jy < MY; jy++) {
    yloc = NSMX * jy;
    for (jx = 0; jx < MX; jx++) {
      loc = yloc + NUM_SPECIES * jx;
      for (is = 0; is < NUM_SPECIES; is++) {
        if (is < np)
          resv[loc+is] = cpv[loc+is] - resv[loc+is];
        else
          resv[loc+is] = -resv[loc+is];
      } /* End is (species) loop */
    } /* End of jx loop */
  } /* End of jy loop */
  
  return(0);
  
}

/* 
 * Fweb: Rate function for the food-web problem.                        
 * This routine computes the right-hand sides of the system equations,   
 * consisting of the diffusion term and interaction term.                
 * The interaction term is computed by the function WebRates.            
 */

static void Fweb(realtype tcalc, N_Vector cc, N_Vector crate,  
                 UserData webdata)
{ 
  long int jx, jy, is, idyu, idyl, idxu, idxl;
  realtype xx, yy, *cxy, *ratesxy, *cratexy, dcyli, dcyui, dcxli, dcxui;
  
  /* Loop over grid points, evaluate interaction vector (length ns),
     form diffusion difference terms, and load crate.                    */
  
  for (jy = 0; jy < MY; jy++) {
    yy = (webdata->dy) * jy ;
    idyu = (jy!=MY-1) ? NSMX : -NSMX;
    idyl = (jy!= 0  ) ? NSMX : -NSMX;
    
    for (jx = 0; jx < MX; jx++) {
      xx = (webdata->dx) * jx;
      idxu = (jx!= MX-1) ?  NUM_SPECIES : -NUM_SPECIES;
      idxl = (jx!=  0  ) ?  NUM_SPECIES : -NUM_SPECIES;
      cxy = IJ_Vptr(cc,jx,jy);
      ratesxy = IJ_Vptr(webdata->rates,jx,jy);
      cratexy = IJ_Vptr(crate,jx,jy);
      
      /* Get interaction vector at this grid point. */
      WebRates(xx, yy, cxy, ratesxy, webdata);
      
      /* Loop over species, do differencing, load crate segment. */
      for (is = 0; is < NUM_SPECIES; is++) {
        
        /* Differencing in y. */
        dcyli = *(cxy+is) - *(cxy - idyl + is) ;
        dcyui = *(cxy + idyu + is) - *(cxy+is);
        
        /* Differencing in x. */
        dcxli = *(cxy+is) - *(cxy - idxl + is);
        dcxui = *(cxy + idxu +is) - *(cxy+is);
        
        /* Compute the crate values at (xx,yy). */
        cratexy[is] = coy[is] * (dcyui - dcyli) +
          cox[is] * (dcxui - dcxli) + ratesxy[is];
        
      } /* End is loop */
    } /* End of jx loop */
  } /* End of jy loop */
  
}

/* 
 * WebRates: Evaluate reaction rates at a given spatial point.
 * At a given (x,y), evaluate the array of ns reaction terms R.
 */

static void WebRates(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy,
                     UserData webdata)
{
  int is;
  realtype fac;
  
  for (is = 0; is < NUM_SPECIES; is++)
    ratesxy[is] = dotprod(NUM_SPECIES, cxy, acoef[is]);
  
  fac = ONE + ALPHA*xx*yy + BETA*sin(FOURPI*xx)*sin(FOURPI*yy);
  
  for (is = 0; is < NUM_SPECIES; is++)  
    ratesxy[is] = cxy[is]*( bcoef[is]*fac + ratesxy[is] );
  
}

/*
 * dotprod: dot product routine for realtype arrays, for use by WebRates.    
 */

static realtype dotprod(long int size, realtype *x1, realtype *x2)
{
  long int i;
  realtype *xx1, *xx2, temp = ZERO;
  
  xx1 = x1; xx2 = x2;
  for (i = 0; i < size; i++) temp += (*xx1++) * (*xx2++);
  return(temp);
  
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  return(0);
}
