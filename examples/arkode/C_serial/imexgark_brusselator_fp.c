/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Example problem:
 * 
 * The following test simulates a brusselator problem from chemical 
 * kinetics.  This is an ODE system with 3 components, Y = [u,v,w], 
 * satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions 
 * Y0 = [u0,v0,w0]. 
 * 
 * We have 3 different testing scenarios:
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 *    Here, all three components exhibit a rapid transient change 
 *    during the first 0.2 time units, followed by a slow and 
 *    smooth evolution.
 *
 * Test 2:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
 *    Here, w experiences a fast initial transient, jumping 0.5 
 *    within a few steps.  All values proceed smoothly until 
 *    around t=6.5, when both u and v undergo a sharp transition, 
 *    with u increaseing from around 0.5 to 5 and v decreasing 
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat 
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * Test 3:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3,  ep=5.0e-4
 *    Here, all components undergo very rapid initial transients 
 *    during the first 0.3 time units, and all then proceed very 
 *    smoothly for the remainder of the simulation.
 *
 * This file is hard-coded to use test 3.
 * 
 * This program solves the problem with the ARK method, using an
 * accelerated fixed-point iteration for the nonlinear solver.
 *
 * 100 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>               /* prototypes for ARKode fcts., consts.      */
#include <arkode/arkode_imexgarkstep.h>  /* prototypes for IMEXGARKStep fcts., consts */
#include <nvector/nvector_serial.h>      /* serial N_Vector types, fcts., macros      */
#include <sundials/sundials_types.h>     /* def. of type 'realtype'                   */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied Functions Called by the Solver */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);     /* initial time */
  realtype Tf = RCONST(10.0);    /* final time */
  realtype dTout = RCONST(1.0);  /* time between outputs */
  sunindextype NEQ = 3;          /* number of dependent vars. */
  int Nt = ceil(Tf/dTout);       /* number of output times */
  int test = 3;                  /* test problem to run */
  realtype reltol = 1.0e-6;      /* tolerances */
  realtype abstol = 1.0e-10;
  int fp_m = 3;                  /* dimension of acceleration subspace */
  int maxcor = 10;               /* maximum # of nonlinear iterations/step */
  realtype a, b, ep, u0, v0, w0;
  realtype rdata[3];

  /* Butcher tables */
  int stages = 2;
  realtype *Aee, *Aei, *Aie, *Aii;
  realtype *ce, *ci, *be, *bi, *de, *di;
  realtype gamma = RCONST(1.0) - RCONST(sqrt(2))/RCONST(2.0);

  /* general problem variables */
  int flag;                      /* reusable error-checking flag */
  N_Vector y = NULL;             /* empty vector for storing solution */
  void *arkode_mem = NULL;       /* empty ARKode memory structure */
  FILE *UFID;
  realtype t, tout;
  int iout;
  long int nst, nst_a, nfe, nfi, nni, ncfn, netf;

  /* set up the test problem according to the desired test */
  if (test == 1) {
    u0 = RCONST(3.9);
    v0 = RCONST(1.1);
    w0 = RCONST(2.8);
    a  = RCONST(1.2);
    b  = RCONST(2.5);
    ep = RCONST(1.0e-5);
  } else if (test == 3) {
    u0 = RCONST(3.0);
    v0 = RCONST(3.0);
    w0 = RCONST(3.5);
    a  = RCONST(0.5);
    b  = RCONST(3.0);
    ep = RCONST(5.0e-4);
  } else {
    u0 = RCONST(1.2);
    v0 = RCONST(3.1);
    w0 = RCONST(3.0);
    a  = RCONST(1.0);
    b  = RCONST(3.5);
    ep = RCONST(5.0e-6);
  }

  /* Initial problem output */
  printf("\nBrusselator ODE test problem, fixed-point solver:\n");
  printf("    initial conditions:  u0 = %"GSYM",  v0 = %"GSYM",  w0 = %"GSYM"\n",u0,v0,w0);
  printf("    problem parameters:  a = %"GSYM",  b = %"GSYM",  ep = %"GSYM"\n",a,b,ep);
  printf("    reltol = %.1"ESYM",  abstol = %.1"ESYM"\n\n",reltol,abstol);

  /* Initialize data structures */
  rdata[0] = a;    /* set user data  */
  rdata[1] = b;
  rdata[2] = ep;
  y = N_VNew_Serial(NEQ);           /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = u0;               /* Set initial conditions */
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;
  arkode_mem = ARKodeCreate();      /* Create the solver memory */
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;

  /* Call ARKStepCreate to initialize the ARK timestepper module and 
     specify the right-hand side functions in y'=fe(t,y)+fi(t,y), 
     the inital time T0, and the initial dependent variable vector y. */
  flag = IMEXGARKStepCreate(arkode_mem, fe, fi, T0, y);
  if (check_flag(&flag, "IMEXARKStepCreate", 1)) return 1;

  Aee = (realtype *) calloc( stages*stages, sizeof(realtype*) );
  be  = (realtype *) calloc( stages, sizeof(realtype) );
  ce  = (realtype *) calloc( stages, sizeof(realtype) );

  Aee[0] = RCONST(0.0);
  Aee[1] = RCONST(0.0);
  Aee[2] = RCONST(1.0) / (RCONST(2.0) * gamma);
  Aee[3] = RCONST(0.0);

  be[0] = RCONST(1.0) - gamma;
  be[1] = gamma;

  ce[0] = RCONST(0.0);
  ce[1] = RCONST(1.0) / gamma;

  Aei = (realtype *) calloc( stages*stages, sizeof(realtype*) );
  Aei[0] = RCONST(0.0);
  Aei[1] = RCONST(0.0);
  Aei[2] = RCONST(1.0) / (RCONST(2.0) * gamma);
  Aei[3] = RCONST(0.0);

  Aii = (realtype *) calloc( stages*stages, sizeof(realtype*) );
  bi  = (realtype *) calloc( stages, sizeof(realtype) );
  ci  = (realtype *) calloc( stages, sizeof(realtype) );

  Aii[0] = gamma;
  Aii[1] = RCONST(0.0);
  Aii[2] = RCONST(1.0) - gamma;
  Aii[3] = gamma;

  bi[0] = RCONST(1.0) - gamma;
  bi[1] = gamma;

  ci[0] = gamma;
  ci[1] = RCONST(1.0);

  Aie = (realtype *) calloc( stages*stages, sizeof(realtype*) );
  Aie[0] = gamma;
  Aie[1] = RCONST(0.0);
  Aie[2] = RCONST(1.0) - gamma;
  Aie[3] = gamma;

  flag = IMEXGARKStepSetButcherTables(arkode_mem, stages, 2, -1,
                                      ce, ci,
                                      Aee, Aei,
                                      Aie, Aii,
                                      be, bi, 
                                      NULL, NULL);

  if (check_flag(&flag, "IMEXGARKStepSetButcherTables", 1)) return 1;

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem, (void *) rdata);     /* Pass rdata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);    /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  flag = IMEXGARKStepSetFixedPoint(arkode_mem, fp_m);             /* Specify fixed-point solver */
  if (check_flag(&flag, "IMEXGARKStepSetFixedPoint", 1)) return 1;
  flag = IMEXGARKStepSetMaxNonlinIters(arkode_mem, maxcor);       /* Increase default iterations */
  if (check_flag(&flag, "IMEXGARKStepSetMaxNonlinIters", 1)) return 1;

  flag = ARKodeSetInitStep(arkode_mem, RCONST(0.01));
  if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt","w");
  fprintf(UFID,"# t u v w\n");

  /* output initial condition to disk */
  fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM"\n", 
          T0, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));  

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  tout = T0+dTout;
  printf("        t           u           v           w\n");
  printf("   ----------------------------------------------\n");
  for (iout=0; iout<Nt; iout++) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);      /* call integrator */
    if (check_flag(&flag, "ARKode", 1)) break;
    printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"\n",             /* access/print solution */
           t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
    fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM"\n", 
            t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));  
    if (flag >= 0) {                                         /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                 /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);

  /* Print some final statistics */
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);

  flag = IMEXGARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKStepGetNumStepAttempts", 1);
  flag = IMEXGARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1);
  flag = IMEXGARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKStepGetNumErrTestFails", 1);
  flag = IMEXGARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1);
  flag = IMEXGARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total number of fixed-point iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  IMEXGARKStepWriteButcher(arkode_mem, stdout);

  /* Clean up and return with successful completion */
  N_VDestroy(y);        /* Free y vector */
  ARKodeFree(&arkode_mem);     /* Free integrator memory */
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* fi routine to compute the implicit portion of the ODE RHS. */
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
  realtype b  = rdata[1];                     /* access data entries */
  realtype ep = rdata[2];
  realtype w = NV_Ith_S(y,2);                 /* access solution values */

  /* fill in the RHS function */
  NV_Ith_S(ydot,0) = 0.0;
  NV_Ith_S(ydot,1) = 0.0;
  NV_Ith_S(ydot,2) = (b-w)/ep;

  return 0;                                  /* Return with success */
}

/* fe routine to compute the explicit portion of the ODE RHS. */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
  realtype a  = rdata[0];                     /* access data entries */
  realtype u = NV_Ith_S(y,0);                 /* access solution values */
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* fill in the RHS function */
  NV_Ith_S(ydot,0) = a - (w+1.0)*u + v*u*u;
  NV_Ith_S(ydot,1) = w*u - v*u*u;
  NV_Ith_S(ydot,2) = -w*u;

  return 0;                                  /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer  
*/
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/

  /* /\* ERK *\/ */
  /* Aee = (realtype *) calloc( stages*stages, sizeof(realtype*) ); */
  /* Aie = (realtype *) calloc( stages*stages, sizeof(realtype*) ); */
  /* Aei = (realtype *) calloc( stages*stages, sizeof(realtype*) ); */
  /* be  = (realtype *) calloc( stages, sizeof(realtype) ); */
  /* ce  = (realtype *) calloc( stages, sizeof(realtype) ); */
  /* de  = (realtype *) calloc( stages, sizeof(realtype) ); */

  /* Aee[1*stages+0] = Aie[1*stages+0] = RCONST(0.5); */
  /* Aee[2*stages+0] = Aie[2*stages+0] = RCONST(13861.0)/RCONST(62500.0); */
  /* Aee[2*stages+1] = Aie[2*stages+1] = RCONST(6889.0)/RCONST(62500.0); */
  /* Aee[3*stages+0] = Aie[3*stages+0] = RCONST(-116923316275.0)/RCONST(2393684061468.0); */
  /* Aee[3*stages+1] = Aie[3*stages+1] = RCONST(-2731218467317.0)/RCONST(15368042101831.0); */
  /* Aee[3*stages+2] = Aie[3*stages+2] = RCONST(9408046702089.0)/RCONST(11113171139209.0); */
  /* Aee[4*stages+0] = Aie[4*stages+0] = RCONST(-451086348788.0)/RCONST(2902428689909.0); */
  /* Aee[4*stages+1] = Aie[4*stages+1] = RCONST(-2682348792572.0)/RCONST(7519795681897.0); */
  /* Aee[4*stages+2] = Aie[4*stages+2] = RCONST(12662868775082.0)/RCONST(11960479115383.0); */
  /* Aee[4*stages+3] = Aie[4*stages+3] = RCONST(3355817975965.0)/RCONST(11060851509271.0); */
  /* Aee[5*stages+0] = Aie[5*stages+0] = RCONST(647845179188.0)/RCONST(3216320057751.0); */
  /* Aee[5*stages+1] = Aie[5*stages+1] = RCONST(73281519250.0)/RCONST(8382639484533.0); */
  /* Aee[5*stages+2] = Aie[5*stages+2] = RCONST(552539513391.0)/RCONST(3454668386233.0); */
  /* Aee[5*stages+3] = Aie[5*stages+3] = RCONST(3354512671639.0)/RCONST(8306763924573.0); */
  /* Aee[5*stages+4] = Aie[5*stages+4] = RCONST(4040.0)/RCONST(17871.0); */

  /* be[0] = RCONST(82889.0)/RCONST(524892.0); */
  /* be[2] = RCONST(15625.0)/RCONST(83664.0); */
  /* be[3] = RCONST(69875.0)/RCONST(102672.0); */
  /* be[4] = RCONST(-2260.0)/RCONST(8211.0); */
  /* be[5] = RCONST(1.0)/RCONST(4.0); */

  /* de[0] = RCONST(4586570599.0)/RCONST(29645900160.0); */
  /* de[2] = RCONST(178811875.0)/RCONST(945068544.0); */
  /* de[3] = RCONST(814220225.0)/RCONST(1159782912.0); */
  /* de[4] = RCONST(-3700637.0)/RCONST(11593932.0); */
  /* de[5] = RCONST(61727.0)/RCONST(225920.0); */

  /* ce[1] = RCONST(1.0)/RCONST(2.0); */
  /* ce[2] = RCONST(83.0)/RCONST(250.0); */
  /* ce[3] = RCONST(31.0)/RCONST(50.0); */
  /* ce[4] = RCONST(17.0)/RCONST(20.0); */
  /* ce[5] = RCONST(1.0); */

  /* /\* ESDIRK *\/ */
  /* Aii = (realtype *) calloc( stages*stages, sizeof(realtype*) ); */
  /* bi  = (realtype *) calloc( stages, sizeof(realtype) ); */
  /* ci  = (realtype *) calloc( stages, sizeof(realtype) ); */
  /* di  = (realtype *) calloc( stages, sizeof(realtype) ); */

  /* Aii[1*stages+0] = RCONST(1.0)/RCONST(4.0); */
  /* Aii[1*stages+1] = RCONST(1.0)/RCONST(4.0); */
  /* Aii[2*stages+0] = RCONST(8611.0)/RCONST(62500.0); */
  /* Aii[2*stages+1] = RCONST(-1743.0)/RCONST(31250.0); */
  /* Aii[2*stages+2] = RCONST(1.0)/RCONST(4.0); */
  /* Aii[3*stages+0] = RCONST(5012029.0)/RCONST(34652500.0); */
  /* Aii[3*stages+1] = RCONST(-654441.0)/RCONST(2922500.0); */
  /* Aii[3*stages+2] = RCONST(174375.0)/RCONST(388108.0); */
  /* Aii[3*stages+3] = RCONST(1.0)/RCONST(4.0); */
  /* Aii[4*stages+0] = RCONST(15267082809.0)/RCONST(155376265600.0); */
  /* Aii[4*stages+1] = RCONST(-71443401.0)/RCONST(120774400.0); */
  /* Aii[4*stages+2] = RCONST(730878875.0)/RCONST(902184768.0); */
  /* Aii[4*stages+3] = RCONST(2285395.0)/RCONST(8070912.0); */
  /* Aii[4*stages+4] = RCONST(1.0)/RCONST(4.0); */
  /* Aii[5*stages+0] = RCONST(82889.0)/RCONST(524892.0); */
  /* Aii[5*stages+2] = RCONST(15625.0)/RCONST(83664.0); */
  /* Aii[5*stages+3] = RCONST(69875.0)/RCONST(102672.0); */
  /* Aii[5*stages+4] = RCONST(-2260.0)/RCONST(8211.0); */
  /* Aii[5*stages+5] = RCONST(1.0)/RCONST(4.0); */

  /* bi[0] = RCONST(82889.0)/RCONST(524892.0); */
  /* bi[2] = RCONST(15625.0)/RCONST(83664.0); */
  /* bi[3] = RCONST(69875.0)/RCONST(102672.0); */
  /* bi[4] = RCONST(-2260.0)/RCONST(8211.0); */
  /* bi[5] = RCONST(1.0)/RCONST(4.0); */

  /* ci[1] = RCONST(1.0)/RCONST(2.0); */
  /* ci[2] = RCONST(83.0)/RCONST(250.0); */
  /* ci[3] = RCONST(31.0)/RCONST(50.0); */
  /* ci[4] = RCONST(17.0)/RCONST(20.0); */
  /* ci[5] = RCONST(1.0); */

  /* di[0] = RCONST(4586570599.0)/RCONST(29645900160.0); */
  /* di[2] = RCONST(178811875.0)/RCONST(945068544.0); */
  /* di[3] = RCONST(814220225.0)/RCONST(1159782912.0); */
  /* di[4] = RCONST(-3700637.0)/RCONST(11593932.0); */
  /* di[5] = RCONST(61727.0)/RCONST(225920.0); */
