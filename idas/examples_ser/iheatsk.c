/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2004-07-22 23:01:56 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem for IDA/IDAS: 2D heat equation, serial, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver IDASpgmr.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
 * PDE is treated with central differences on a uniform M x M grid.
 * The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = M^2. Here M = 10.
 *
 * The system is solved with IDA/IDAS using the Krylov linear solver
 * IDASPGMR. The preconditioner uses the diagonal elements of the
 * Jacobian only. Routines for preconditioning, required by
 * IDASPGMR, are supplied here. The constraints u >= 0 are posed
 * for all components. Output is taken at t = 0, .01, .02, .04,
 * ..., 10.24. Two cases are run -- with the Gram-Schmidt type
 * being Modified in the first case, and Classical in the second.
 * The second run uses IDAReInit and IDAReInitSpgmr.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sundialstypes.h"
#include "sundialsmath.h"
#include "nvector_serial.h"
#include "idas.h"
#include "idaspgmr.h"
#include "iterative.h"

typedef struct {  
  /* mm is the number of points in the grid. */
  long int mm;
  realtype dx;
  realtype coeff;
  /* pp is the vector of diagonal preconditioner elements. */
  N_Vector pp;
} *UserData;

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res);

int heatres(realtype tres, N_Vector uu, N_Vector up,
            N_Vector resval, void *rdata);

int PrecondHeateq(realtype tt, N_Vector uu,
                  N_Vector up, N_Vector rr, realtype cj,
                  void *pdata, 
                  N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

int PSolveHeateq(realtype tt, N_Vector uu,
                 N_Vector up, N_Vector rr, 
                 N_Vector rvec, N_Vector zvec,
                 realtype cj, realtype delta,
                 void *pdata, N_Vector tempv);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

#define NOUT  11
#define MGRID 10
#define NEQ   (MGRID*MGRID)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define FOUR  RCONST(4.0)

int main()
{
  void *mem;
  UserData data;
  N_Vector uu, up, constraints, id, res;
  int ier, iout, itol, itask, kused;
  realtype rtol, atol, t0, t1, tout, tret, umax, hused;
  long int nst, nni, nje, nre, nreS, nli, npe, nps, netf, ncfn, ncfl;

  mem = NULL;
  data = NULL;
  uu = up = constraints = id = res = NULL;

  /* Allocate N-vectors and the user data structure. */

  uu = N_VNew_Serial(NEQ);
  if(check_flag((void *)uu, "N_VNew_Serial", 0)) return(1);
  up = N_VNew_Serial(NEQ);
  if(check_flag((void *)up, "N_VNew_Serial", 0)) return(1);
  res = N_VNew_Serial(NEQ);
  if(check_flag((void *)res, "N_VNew_Serial", 0)) return(1);
  constraints = N_VNew_Serial(NEQ);
  if(check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);
  id = N_VNew_Serial(NEQ);
  if(check_flag((void *)id, "N_VNew_Serial", 0)) return(1);

  data = (UserData) malloc(sizeof *data);
  data->pp = NULL;
  if(check_flag((void *)data, "malloc", 2)) return(1);

  /* Assign parameters in the user data structure. */

  data->mm  = MGRID;
  data->dx = ONE/(MGRID-ONE);
  data->coeff = ONE/(data->dx * data->dx);
  data->pp = N_VNew_Serial(NEQ);
  if(check_flag((void *)data->pp, "N_VNew_Serial", 0)) return(1);

  /* Initialize uu, up, id. */
  SetInitialProfile(data, uu, up, id, res);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);

  /* Assign various parameters. */
  t0   = ZERO;
  t1   = 0.01;
  itol = SS; /* Specify scalar relative and absolute tolerance. */
  rtol = ZERO;
  atol = 1.e-3; 
  itask = NORMAL;

  /* Call IDACreate and IDAMalloc to initialize solution */
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);
  ier = IDASetRdata(mem, data);
  if(check_flag(&ier, "IDASetRdata", 1)) return(1);
  ier = IDASetId(mem, id);
  if(check_flag(&ier, "IDASetId", 1)) return(1);
  ier = IDASetConstraints(mem, constraints);
  if(check_flag(&ier, "IDASetConstraints", 1)) return(1);
  ier = IDAMalloc(mem, heatres, t0, uu, up, itol, &rtol, &atol);
  if(check_flag(&ier, "IDAMalloc", 1)) return(1);

  /* Call IDASpgmr to specify the linear solver. */

  ier = IDASpgmr(mem, 0.0);
  if(check_flag(&ier, "IDASpgmr", 1)) return(1);
  ier = IDASpgmrSetPrecSetupFn(mem, PrecondHeateq);
  if(check_flag(&ier, "IDASpgmrSetPrecSetupFn", 1)) return(1);
  ier = IDASpgmrSetPrecSolveFn(mem, PSolveHeateq);
  if(check_flag(&ier, "IDASpgmrSetPrecSolveFn", 1)) return(1);
  ier = IDASpgmrSetPrecData(mem, data);
  if(check_flag(&ier, "IDASpgmrSetPrecData", 1)) return(1);

  /* Print output heading. */
  
  printf("iheatsk: Heat equation, serial example problem for IDA \n");
  printf("         Discretized heat equation on 2D unit square. \n");
  printf("         Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("         Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("        Total system size: %d\n\n", NEQ);
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: IDASPGMR, preconditioner using diagonal elements. \n");

  /* Case 1 */

  /* Print case number, output table heading, and initial line of table. */
  printf("\n\nCase 1: gsytpe = MODIFIED_GS\n");
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nje   nre   nreS     h      npe nps\n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");
  umax = N_VMaxNorm(uu);
  
  printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e  %3d %3d\n",
         t0, umax, 0, 0, 0, 0, 0, 0, 0.0, 0, 0);

  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1,iout = 1; iout <= NOUT ; iout++, tout *= TWO) {

    ier = IDASolve(mem, tout, &tret, uu, up, itask);
    if(check_flag(&ier, "IDASolve", 1)) return(1);

    umax = N_VMaxNorm(uu);

    ier = IDAGetLastOrder(mem, &kused);
    check_flag(&ier, "IDAGetLastOrder", 1);
    ier = IDAGetNumSteps(mem, &nst);
    check_flag(&ier, "IDAGetNumSteps", 1);
    ier = IDAGetNumNonlinSolvIters(mem, &nni);
    check_flag(&ier, "IDAGetNumNonlinSolvIters", 1);
    ier = IDAGetNumResEvals(mem, &nre);
    check_flag(&ier, "IDAGetNumResEvals", 1);
    ier = IDAGetLastStep(mem, &hused);
    check_flag(&ier, "IDAGetLastStep", 1);
    ier = IDASpgmrGetNumJtimesEvals(mem, &nje);
    check_flag(&ier, "IDASpgmrGetNumJtimesEvals", 1);
    ier = IDASpgmrGetNumLinIters(mem, &nli);
    check_flag(&ier, "IDASpgmrGetNumLinIters", 1);
    ier = IDASpgmrGetNumResEvals(mem, &nreS);
    check_flag(&ier, "IDASpgmrGetNumResEvals", 1);
    ier = IDASpgmrGetNumPrecEvals(mem, &npe);
    check_flag(&ier, "IDASpgmrGetNumPrecEvals", 1);
    ier = IDASpgmrGetNumPrecSolves(mem, &nps);
    check_flag(&ier, "IDASpgmrGetNumPrecSolves", 1);

    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           tret, umax, kused, nst, nni, nje, nre, nreS, hused, npe, nps);

    } /* End of tout loop. */

  /* Print remaining counters. */

  ier = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&ier, "IDAGetNumErrTestFails", 1);
  ier = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&ier, "IDAGetNumNonlinSolvConvFails", 1);
  ier = IDASpgmrGetNumConvFails(mem, &ncfl);
  check_flag(&ier, "IDASpgmrGetNumConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld,   ncfl = %ld \n", netf, ncfn, ncfl);

  /* Case 2. */

  /* Re-initialize uu, up, id. */
  SetInitialProfile(data, uu, up, id, res);
  
  /* Re-initialize IDA and IDASPGMR */
  ier = IDAReInit(mem, heatres, t0, uu, up, itol, &rtol, &atol);
  if(check_flag(&ier, "IDAReInit", 1)) return(1);
  
  ier = IDASpgmrSetGSType(mem, CLASSICAL_GS);
  if(check_flag(&ier, "IDASpgmrSetGSType",1)) return(1); 
  
  /* Print case number, output table heading, and initial line of table. */
  printf("\n\n\n\nCase 2: gstype = CLASSICAL_GS\n");
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nje   nre   nreS     h      npe nps\n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .\n");
  umax = N_VMaxNorm(uu);
  
  printf(" %5.2f %13.5e  %d  %3d  %3d  %3d  %4d  %4d  %9.2e  %3d %3d\n",
         t0, umax, 0, 0, 0, 0, 0, 0, 0.0, 0, 0);

   
  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1,iout = 1; iout <= NOUT ; iout++, tout *= TWO) {

    ier = IDASolve(mem, tout, &tret, uu, up, itask);
    if(check_flag(&ier, "IDASolve", 1)) return(1);

    umax = N_VMaxNorm(uu);

    ier = IDAGetLastOrder(mem, &kused);
    check_flag(&ier, "IDAGetLastOrder", 1);
    ier = IDAGetNumSteps(mem, &nst);
    check_flag(&ier, "IDAGetNumSteps", 1);
    ier = IDAGetNumNonlinSolvIters(mem, &nni);
    check_flag(&ier, "IDAGetNumNonlinSolvIters", 1);
    ier = IDAGetNumResEvals(mem, &nre);
    check_flag(&ier, "IDAGetNumResEvals", 1);
    ier = IDAGetLastStep(mem, &hused);
    check_flag(&ier, "IDAGetLastStep", 1);
    ier = IDASpgmrGetNumJtimesEvals(mem, &nje);
    check_flag(&ier, "IDASpgmrGetNumJtimesEvals", 1);
    ier = IDASpgmrGetNumLinIters(mem, &nli);
    check_flag(&ier, "IDASpgmrGetNumLinIters", 1);
    ier = IDASpgmrGetNumResEvals(mem, &nreS);
    check_flag(&ier, "IDASpgmrGetNumResEvals", 1);
    ier = IDASpgmrGetNumPrecEvals(mem, &npe);
    check_flag(&ier, "IDASpgmrGetPrecEvals", 1);
    ier = IDASpgmrGetNumPrecSolves(mem, &nps);
    check_flag(&ier, "IDASpgmrGetNumPrecSolves", 1);

    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           tret, umax, kused, nst, nni, nje, nre, nreS, hused, npe, nps);

  } /* End of tout loop. */

  /* Print remaining counters. */

  ier = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&ier, "IDAGetNumErrTestFails", 1);
  ier = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&ier, "IDAGetNumNonlinSolvConvFails", 1);
  ier = IDASpgmrGetNumConvFails(mem, &ncfl);
  check_flag(&ier, "IDASpgmrGetNumConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld,   ncfl = %ld \n", netf, ncfn, ncfl);

  /* Free Memory */

  IDAFree(mem);
  N_VDestroy(uu);
  N_VDestroy(up);
  N_VDestroy(constraints);
  N_VDestroy(id);
  N_VDestroy(res);
  N_VDestroy(data->pp);
  free(data);

  return(0);
}

/*
 * SetInitialProfile: routine to initialize u, up, and id vectors.       
 */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector id, N_Vector res)
{
  long int mm, mm1, i, j, offset, loc;
  realtype xfact, yfact, *udata, *updata, *iddata;

  mm = data->mm;

  udata = NV_DATA_S(uu);
  updata = NV_DATA_S(up);
  iddata = NV_DATA_S(id);

  /* Initialize id to 1's. */
  N_VConst(ONE, id);

  /* Initialize uu on all grid points. */ 
  mm1 = mm - 1;
  for (j = 0; j < mm; j++) {
    yfact = data->dx * j;
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      xfact = data->dx * i;
      loc = offset + i;
      udata[loc] = 16. * xfact * (ONE - xfact) * yfact * (ONE - yfact);
    }
  }
  
  /* Initialize up vector to 0. */
  N_VConst(ZERO, up);

  /* heatres sets res to negative of ODE RHS values at interior points. */
  heatres(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct interior initial up values. */
  N_VScale(-ONE, res, up);

  /* Set up and id at boundary points to zero. */
  for (j = 0; j < mm; j++) {
    offset = mm*j;
    for (i = 0; i < mm; i++) {
      loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) {
        updata[loc] = ZERO; iddata[loc] = ZERO; }
    }
  }
  
  return(SUCCESS);
  
}

/*
 * heatres: heat equation system residual function (user-supplied)      
 * This uses 5-point central differencing on the interior points, and   
 * includes algebraic equations for the boundary values.                
 * So for each interior point, the residual component has the form      
 *    res_i = u'_i - (central difference)_i                             
 * while for each boundary point, it is res_i = u_i.                     
 */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector res, void *rdata)
{
  long int i, j, offset, loc, mm;
  realtype *uv, *upv, *resv, coeff;
  UserData data;
  
  uv = NV_DATA_S(uu); upv = NV_DATA_S(up); resv = NV_DATA_S(res);

  data = (UserData)rdata;
  
  coeff = data->coeff;
  mm    = data->mm;
  
  /* Initialize res to uu, to take care of boundary equations. */
  N_VScale(ONE, uu, res);
  
  /* Loop over interior points; set res = up - (central difference). */
  for (j = 1; j < MGRID-1; j++) {
    offset = mm*j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      resv[loc]= upv[loc] - coeff * 
       (uv[loc-1] + uv[loc+1] + uv[loc-mm] + uv[loc+mm] - 4*uv[loc]);
    }
  }

  return(SUCCESS);

}

/*
 * PrecondHeateq: setup for diagonal preconditioner for iheatsk.   
 *                                                                 
 * The optional user-supplied functions PrecondHeateq and          
 * PSolveHeateq together must define the left preconditoner        
 * matrix P approximating the system Jacobian matrix               
 *                   J = dF/du + cj*dF/du'                         
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   
 * systems P z = r.   This is done in this case by keeping only    
 * the diagonal elements of the J matrix above, storing them as    
 * inverses in a vector pp, when computed in PrecondHeateq, for    
 * subsequent use in PSolveHeateq.                                 
 *                                                                 
 * In this instance, only cj and data (user data structure, with    
 * pp etc.) are used from the PrecondHeateq argument list.         
 */
  
int PrecondHeateq(realtype tt, N_Vector uu,
                  N_Vector up, N_Vector rr, realtype cj,
                  void *pdata, 
                  N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  
  long int i, j, offset, loc, mm;
  realtype *ppv, pelinv;
  UserData data;
  
  data = (UserData) pdata;
  ppv = NV_DATA_S(data->pp);
  mm = data->mm;

  /* Initialize the entire vector to 1., then set the interior points to the
     correct value for preconditioning. */
  
  N_VConst(ONE,data->pp);
  
  /* Compute the inverse of the preconditioner diagonal elements. */
  pelinv = ONE/(cj + FOUR*data->coeff); 
  
  for (j = 1; j < mm-1; j++) {
    offset = mm * j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      ppv[loc] = pelinv;
    }
  }
  
  return(SUCCESS);
  
}

/*
 * PSolveHeateq: solve preconditioner linear system.              
 * This routine multiplies the input vector rvec by the vector pp 
 * containing the inverse diagonal Jacobian elements (previously  
 * computed in PrecondHeateq), returning the result in zvec.      
 */

int PSolveHeateq(realtype tt, N_Vector uu,
                 N_Vector up, N_Vector rr, 
                 N_Vector rvec, N_Vector zvec,
                 realtype cj, realtype delta,
                 void *pdata, N_Vector tempv)
{
  UserData data;
  
  data = (UserData) pdata;
  
  N_VProd(data->pp, rvec, zvec);
  
  return(SUCCESS);

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
