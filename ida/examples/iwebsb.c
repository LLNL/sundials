/*************************************************************************
 * File:       iwebsb.c                                                  *
 * Written by: Allan G. Taylor and Alan C. Hindmarsh @ LLNL              *
 * Version of: 10 January 2001                                           *
 *-----------------------------------------------------------------------*
 *
 * Example program for IDA: Food web, serial, band solve IDABAND.
 *
 * This example program for IDA (serial version) uses IDABAND as the
 * linear solver, and IDACalcIC for initial condition calculation.
 *                                         
 * The mathematical problem solved in this example is a DAE system that 
 * arises from a system of partial differential equations after spatial
 * discretization.  The PDE system is a food web population model, with
 * predator-prey interaction and diffusion on the unit square in two 
 * dimensions. The dependent variable vector is:
 *
 *       1   2         ns
 * c = (c , c ,  ..., c  ) ,   ns = 2 * np
 * 
 * and the PDE's are as follows:
 *
 *     i             i      i
 *   dc /dt = d(i)*(c    + c  )  +  R (x,y,c)    (i=1,...,np)
 *                   xx     yy       i
 *
 *
 *                   i      i      
 *  0       = d(i)*(c    + c  )  +  R  (x,y,c)   (i=np+1,...,ns)
 *                   xx     yy       i
 *
 *   where the reaction terms R are:
 *
 *                   i             ns         j  
 *   R  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being prey and
 * the last np being predators. The coefficients a(i,j), b(i), d(i) are:
 *
 *   a(i,i) = -AA  (all i)
 *   a(i,j) = -GG  (i <= np , j >  np)
 *   a(i,j) =  EE  (i >  np,  j <= np)
 *   all other a(i,j) = 0
 *   b(i) = BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))  (i <= np)
 *   b(i) =-BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))  (i  > np)
 *   d(i) = DPREY  (i <= np)
 *   d(i) = DPRED  (i > np)
 *  
 *  NOTE: The above equations are written in 1-based indices, whereas the
 *  code has 0-based indices, being written in C.
 *
 *  The various scalar parameters required are set using 'define' statements 
 *  or directly in routine InitUserData.  In this program, np = 1, ns = 2.
 *  The boundary conditions are homogeneous Neumann: normal derivative = 0.
 *
 *  A polynomial in x and y is used to set the initial values of the first
 *  np variables (the prey variables) at each x,y location, while initial
 *  values for the remaining (predator) variables are set to a flat value,
 *  which is corrected by IDACalcIC.
 *
 *  The PDEs are discretized by central differencing on a MX by MY mesh.
 * 
 *  The DAE system is solved by IDA using the IDABAND linear solver.
 *  Output is printed at t = 0, .001, .01, .1, .4, .7, 1.
 *
 * References:
 * [1] Peter N. Brown and Alan C. Hindmarsh,
 *     Reduced Storage Matrix Methods in Stiff ODE systems, Journal of
 *     Applied Mathematics and Computation, Vol. 31 (May 1989), pp. 40-91.
 *
 * [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Using Krylov Methods in the Solution of Large-Scale Differential-
 *     Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488.
 * 
 * [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Consistent Initial Condition Calculation for Differential-
 *     Algebraic Systems, SIAM J. Sci. Comput., 19 (1998), pp. 1495-1512.
 * 
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"   /* Definitions of real, integer, boole, TRUE, FALSE.*/
#include "ida.h"        /* Main IDA header file.                            */
#include "idaband.h"    /* Use IDABAND linear solver.                       */
#include "nvector.h"    /* Definitions of type N_Vector, macro N_VDATA.     */
#include "llnlmath.h"   /* Contains RSqrt and UnitRoundoff routines.        */
#include "smalldense.h" /* Contains definitions for denalloc routine.       */

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


/* User-defined vector and accessor macro: IJ_Vptr.
   IJ_Vptr is defined in order to express the underlying 3-D structure of the 
   dependent variable vector from its underlying 1-D storage (an N_Vector).
   IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
   species index is = 0, x-index ix = i, and y-index jy = j.                */

#define IJ_Vptr(vv,i,j)   (&(((vv)->data)[(i)*NUM_SPECIES + (j)*NSMX]))


/* Type: UserData.  Contains problem constants, etc. */

typedef struct {
  integer Neq, ns, np, mx, my;
  real dx, dy, **acoef;
  real cox[NUM_SPECIES], coy[NUM_SPECIES], bcoef[NUM_SPECIES];
  N_Vector rates;
} *UserData;

/* Prototypes for private Helper Functions. */

static UserData AllocUserData(void);
static void InitUserData(UserData webdata);
static void FreeUserData(UserData webdata);
static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
			       UserData webdata);
static void PrintOutput(long int iopt[], real ropt[], N_Vector cc, real time,
                        UserData webdata);
static void PrintFinalStats(long int iopt[]);
static void Fweb(real tcalc, N_Vector cc, N_Vector crate, UserData webdata);
static void WebRates(real xx, real yy, real *cxy, real *ratesxy, 
                     UserData webdata);
static real dotprod(integer size, real *x1, real *x2);


/* Prototypes for functions called by the IDA Solver. */

static int resweb(integer Neq, real time, N_Vector cc, N_Vector cp,
                  N_Vector resval, void *rdata);

/***************************** Main Program ******************************/

main()
{ 
  int iout, flag, retval, i, itol, itask;
  integer SystemSize, mu, ml;
  long int iopt[OPT_SIZE];
  boole optIn;
  real rtol, atol, ropt[OPT_SIZE], t0, tout, tret;
  N_Vector cc, cp, id, res;
  UserData webdata;
  void *mem;
  IDAMem idamem;

  /* Set up and load user data block webdata. */

  SystemSize = NEQ;
  webdata = AllocUserData();
  InitUserData(webdata);

  /* Allocate N-vectors and initialize cc, cp, and id.
     The vector res is used temporarily only.           */ 

  SystemSize = NEQ;
  cc  = N_VNew(SystemSize, NULL);
  cp  = N_VNew(SystemSize, NULL);
  res = N_VNew(SystemSize, NULL);
  id  = N_VNew(SystemSize, NULL);

  SetInitialProfiles(cc, cp, id, webdata);
  N_VFree(res);

  /* Set remaining inputs to IDAMalloc. */

  t0 = ZERO;
  itol = SS; rtol = RTOL; atol = ATOL;
  optIn = FALSE;

  /* Call IDAMalloc to initialize IDA.
  First NULL argument  = constraints vector, not used here.
  Second NULL argument = file pointer for error messages (sent to stdout).
  Third NULL argument = machEnv, not used in serial version.
  A pointer to IDA problem memory is returned and stored in idamem.     */

  mem = IDAMalloc(SystemSize, resweb, webdata, t0, cc, cp, itol,&rtol,&atol,
                  id, NULL, NULL, optIn, iopt, ropt, NULL);

  if (mem == NULL) { printf("IDAMalloc failed."); return(1); }
  idamem = (IDAMem)mem;

  /* Call IDABand to specify the IDA linear solver. */

  mu = ml = NSMX;
  retval = IDABand(idamem, mu, ml, NULL, NULL);

  if (retval != 0) {
    printf("IDABand failed, returning %d \n",retval);
    return(1); }

  /* Call IDACalcIC (with default options) to correct the initial values. */

  tout = 0.001;
  retval = IDACalcIC (idamem, CALC_YA_YDP_INIT, tout, ZERO, 0,0,0,0, ZERO);

  if (retval != SUCCESS) {
    printf("IDACalcIC failed. retval = %d\n",retval);
    return(1); }

  /* Print heading, basic parameters, and initial values. */

  printf("iwebsb: Predator-prey DAE serial example problem for IDA \n\n");
  printf("Number of species ns: %d", NUM_SPECIES);
  printf("     Mesh dimensions: %d x %d", MX, MY);
  printf("     System size: %d\n",SystemSize);
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
  printf("Linear solver: IDABAND,  Band parameters mu = %d, ml = %d\n",mu,ml);
  printf("CalcIC called to correct initial predator concentrations \n\n");

  PrintOutput(iopt, ropt, cc, ZERO, webdata);

  /* Loop over iout, call IDASolve (normal mode), print selected output. */

  itask = NORMAL;
  for (iout = 1; iout <= NOUT; iout++) {

    flag = IDASolve(mem, tout, t0, &tret, cc, cp, itask);

    if(flag != SUCCESS) { 
      printf("IDA failed, flag = %d.\n", flag); 
      return(flag); }

    PrintOutput(iopt, ropt, cc, tret,webdata);

    if (iout < 3) tout *= TMULT; else tout += TADD;

  } /* End of iout loop. */

  /* Print final statistics and free memory. */  

  PrintFinalStats(iopt);
  N_VFree(cc); N_VFree(cp); N_VFree(id);
  IDAFree(idamem);
  FreeUserData(webdata);

  return(0);

}  /* End of iwebsb main program. */


/*********************** Private Helper Functions ************************/


/*************************************************************************/
/* AllocUserData: Allocate memory for data structure of type UserData.   */

static UserData AllocUserData(void)
{
  int jx, jy;
  UserData webdata;

  webdata = (UserData) malloc(sizeof *webdata);

  webdata->rates = N_VNew(NEQ,NULL);

  webdata->acoef = denalloc(NUM_SPECIES);
 
  return(webdata);
}


/* Define lines for readability in later routines */ 

#define acoef  (webdata->acoef)
#define bcoef  (webdata->bcoef)
#define cox    (webdata->cox)
#define coy    (webdata->coy)


/*************************************************************************/
/* InitUserData: Load problem constants in webdata (of type UserData).   */

static void InitUserData(UserData webdata)
{
  int i, j, np;
  real *a1,*a2, *a3, *a4, *b, dx2, dy2;

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

} /* End of InitUserData */


/*************************************************************************/
/* FreeUserData: Free webdata memory.                                    */

static void FreeUserData(UserData webdata)
{
  denfree(acoef);
  N_VFree(webdata->rates);

  free(webdata);
}


/*************************************************************************/
/* SetInitialProfiles: Set initial conditions in cc, cp, and id.
   A polynomial profile is used for the prey cc values, and a constant
   (1.0e5) is loaded as the initial guess for the predator cc values.
   The id values are set to 1 for the prey and 0 for the predators.
   The prey cp values are set according to the given system, and
   the predator cp values are set to zero.                               */

static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
			       UserData webdata)
{
  integer loc, yloc, is, jx, jy, np;
  real xx, yy, xyfactor, fac;
  real *ccv, *cpv, *idv, *crate;

  ccv = N_VDATA(cc);
  cpv = N_VDATA(cp);
  idv = N_VDATA(id);
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
	  ccv[loc+is] = 10. + (real)(is+1) * xyfactor;
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
}  /* End of SetInitialProfiles. */


/*************************************************************************/
/* PrintOutput: Print output values at output time t = tt.
   Selected run statistics are printed.  Then values of c1 and c2
   are printed for the bottom left and top right grid points only.  
   (NOTE: This routine is specific to the case NUM_SPECIES = 2.)         */


static void PrintOutput(long int iopt[], real ropt[], N_Vector cc, real tt,
                        UserData webdata)
{
  int is, jx, jy;
  real  *ct;

  printf("\nTIME t = %e.     NST = %d, k = %d, h = %e\n",
         tt, iopt[NST], iopt[KUSED], ropt[HUSED]);
  printf("NRE = %d, NNI = %d, NJE = %d\n", iopt[NRE], iopt[NNI], 
	 iopt[BAND_NJE]);

  jx = 0;    jy = 0;    ct = IJ_Vptr(cc,jx,jy);
  printf("At bottom left:  c1, c2 = %e %e \n",   ct[0],ct[1]);

  jx = MX-1; jy = MY-1; ct = IJ_Vptr(cc,jx,jy);
  printf("At top right:    c1, c2 = %e %e \n\n", ct[0],ct[1]);

} /* End of PrintOutput. */


/*************************************************************************/
/* PrintFinalStats: Print final run data contained in iopt.              */

static void PrintFinalStats(long int iopt[])
{ 
  printf("\nFinal run statistics: \n\n");
  printf("NST  = %5ld     NRE  = %5ld \n", iopt[NST], iopt[NRE]);
  printf("NNI  = %5ld     NJE  = %5ld \n", iopt[NNI], iopt[BAND_NJE]);
  printf("NETF = %5ld     NCFN = %5ld \n", iopt[NETF], iopt[NCFN]);

} /* End of PrintFinalStats. */


/*********** Functions called by IDA and supporting functions ************/


/*************************************************************************/
/* resweb: System residual function for predator-prey system.            */
/* This routine calls Fweb to get all the right-hand sides of the        */
/* equations, then loads the residual vector accordingly,                */
/* using cp in the case of prey species.                                 */

static int resweb(integer Neq, real tt, N_Vector cc, N_Vector cp, 
                  N_Vector res,  void *rdata)
{
  integer jx, jy, is, yloc, loc, np;
  real *resv, *cpv;
  UserData webdata;
 
  webdata = (UserData)rdata;

  cpv = N_VDATA(cp);
  resv = N_VDATA(res);
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

}  /* End of residual function resweb. */


/*************************************************************************/
/* Fweb: Rate function for the food-web problem.                         */
/* This routine computes the right-hand sides of the system equations,   */
/* consisting of the diffusion term and interaction term.                */
/* The interaction term is computed by the function WebRates.            */

static void Fweb(real tcalc, N_Vector cc, N_Vector crate,  
                UserData webdata)
{ 
  integer jx, jy, is, idyu, idyl, idxu, idxl;
  real xx, yy, *cxy, *ratesxy, *cratexy, dcyli, dcyui, dcxli, dcxui;

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

}  /* End of Fweb. */


/*************************************************************************/
/* WebRates: Evaluate reaction rates at a given spatial point.           */
/* At a given (x,y), evaluate the array of ns reaction terms R.          */

static void WebRates(real xx, real yy, real *cxy, real *ratesxy,
                     UserData webdata)
{
  int is;
  real fac;

  for (is = 0; is < NUM_SPECIES; is++)
       ratesxy[is] = dotprod(NUM_SPECIES, cxy, acoef[is]);

  fac = ONE + ALPHA*xx*yy + BETA*sin(FOURPI*xx)*sin(FOURPI*yy);

  for (is = 0; is < NUM_SPECIES; is++)  
       ratesxy[is] = cxy[is]*( bcoef[is]*fac + ratesxy[is] );

} /* End of WebRates. */


/*************************************************************************/
/* dotprod: dot product routine for real arrays, for use by WebRates.    */

static real dotprod(integer size, real *x1, real *x2)
{
  integer i;
  real *xx1, *xx2, temp = ZERO;

  xx1 = x1; xx2 = x2;
  for (i = 0; i < size; i++) temp += (*xx1++) * (*xx2++);
  return(temp);

} /* End of dotprod. */
