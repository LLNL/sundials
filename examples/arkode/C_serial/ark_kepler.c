/* ----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * We consider the perturbed Kepler problem
 *    dq/dt = [ p1 ]
 *            [ p2 ]
 *    dp/dt = [ -q1 / (q1^2 + q2^2)^(3/2) - delta*q1 / (q1^2 + q2^2)^(5/2) ]
 *          = [ -q2 / (q1^2 + q2^2)^(3/2) - delta*q2 / (q1^2 + q2^2)^(5/2) ]
 * (delta = 0.015) with the initial conditions
 *    q(0) = [ 1 - e ],  p(0) = [        0          ]
 *           [   0   ]          [ sqrt((1+e)/(1-e)) ]
 * where e = 0.6 is the eccentricity.
 *
 * The Hamiltonian for the system,
 *    H(p, q) = 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2)
 *            - 1/200 / (2 * sqrt(q1^2 + q2^2)^3))
 * is conserved as well as the angular momentum,
 *    L(p, q) = q1*p2 - q2*p1.
 *
 * We solve the problem by letting y = [ q, p ]^T then using
 * ARKStep.
 *
 * References:
 *    Ernst Hairer, Christain Lubich, Gerhard Wanner
 *    Geometric Numerical Integration: Structure-Preserving
 *    Algorithms for Ordinary Differential Equations
 *    Springer, 2006,
 *    ISSN 0179-3632
 * ----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <arkode/arkode_mristep.h>      /* prototypes for MRIStep fcts., consts */
#include <arkode/arkode_arkstep.h>      /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>     /* serial N_Vector type, fcts., macros  */
#include <sundials/sundials_math.h>     /* def. math fcns, 'realtype'           */

static int check_retval(void *returnvalue, const char *funcname, int opt);

static void InitialConditions(N_Vector y0, sunrealtype ecc);
static sunrealtype Hamiltonian(N_Vector y);
static sunrealtype AngularMomentum(N_Vector y);

static int dydt(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int dqdt(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int dpdt(realtype t, N_Vector y, N_Vector ydot, void *user_data);

typedef struct {
  sunrealtype ecc;
  sunrealtype delta;
} *UserData;


int main(int argc, char* argv[])
{
  SUNContext sunctx;
  N_Vector y;
  SUNNonlinearSolver NLS;
  UserData udata;
  sunrealtype tout, tret;
  sunrealtype H0, L0;
  ARKodeButcherTable Mp, Mq;
  void* arkode_mem;
  int argi, iout, retval;

  /* Default problem parameters */
  const sunrealtype T0    = SUN_RCONST(0.0);
  sunrealtype Tf          = SUN_RCONST(200.0);
  // sunrealtype Tf          = SUN_RCONST(4000.0);
  const sunrealtype dt    = SUN_RCONST(1e-4);
  const sunrealtype ecc   = SUN_RCONST(0.6);
  const sunrealtype delta = SUN_RCONST(0.015);

  /* Default integrator Options */
  int fixed_step_mode = 1;
  int method          = 1;
  const sunrealtype dTout = SUN_RCONST(10.0);
  // const sunrealtype dTout = SUN_RCONST(100.0);
  const int num_output_times = (int) ceil(Tf/dTout);

  /* Parse CLI args */
  argi = 0;
  if (argc > 1) {
    fixed_step_mode = atoi(argv[++argi]);
  }
  if (argc > 2) {
    method = atoi(argv[++argi]);
  }

  /* Allocate and fill udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->ecc = ecc;
  udata->delta = delta;

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  /* Allocate our state vector */
  y = N_VNew_Serial(4, sunctx);

  /* Fill the initial conditions */
  InitialConditions(y, ecc);

  /* Create ARKStep integrator where we treat dqdt explicitly and dpdt implicitly */
  if (method == 0) {
    arkode_mem = ARKStepCreate(dqdt, dpdt, T0, y, sunctx);

    /* Attach custom Butcher tables for 3rd order scheme [Candy, 1991] */
    Mq = ARKodeButcherTable_Alloc(3, SUNTRUE);
    if (check_retval((void *) Mq, "ARKodeButcherTable_Alloc", 0)) return 1;
    Mq->b[0] = RCONST(2.0)/RCONST(3.0);
    Mq->b[1] = -RCONST(2.0)/RCONST(3.0);
    Mq->b[2] = RCONST(1.0);
    Mq->A[0][0] = Mq->b[0];
    Mq->A[1][0] = Mq->b[0];
    Mq->A[1][1] = Mq->b[1];
    Mq->A[2][0] = Mq->b[0];
    Mq->A[2][1] = Mq->b[1];
    Mq->A[2][2] = Mq->b[2];
    Mq->c[0] = Mq->b[0];
    Mq->c[1] = Mq->b[0] + Mq->b[1];
    Mq->c[2] = Mq->b[0] + Mq->b[1] + Mq->b[2];
    Mq->q = 3;
    Mq->p = 0;

    Mp = ARKodeButcherTable_Alloc(3, SUNTRUE);
    if (check_retval((void *) Mp, "ARKodeButcherTable_Alloc", 0)) return 1;
    Mp->b[0] = RCONST(7.0)/RCONST(24.0);
    Mp->b[1] = RCONST(3.0)/RCONST(4.0);
    Mp->b[2] = -RCONST(1.0)/RCONST(24.0);
    Mp->A[1][0] = Mp->b[0];
    Mp->A[2][0] = Mp->b[0];
    Mp->A[2][1] = Mp->b[1];
    Mp->c[0] = Mp->b[0];
    Mp->c[1] = Mp->b[0] + Mp->b[1];
    Mp->c[2] = Mp->b[0] + Mp->b[1] + Mp->b[2];
    Mp->q = 3;
    Mp->p = 0;

    // Mq = ARKodeButcherTable_Alloc(1, SUNTRUE);
    // if (check_retval((void *) Mq, "ARKodeButcherTable_Alloc", 0)) return 1;
    // Mq->b[0] = RCONST(1.0);
    // Mq->A[0][0] = RCONST(1.0);
    // Mq->c[0] = RCONST(1.0);
    // Mq->q = 1;
    // Mq->p = 0;

    // Mp = ARKodeButcherTable_Alloc(1, SUNTRUE);
    // if (check_retval((void *) Mp, "ARKodeButcherTable_Alloc", 0)) return 1;
    // Mp->b[0] = RCONST(1.0);
    // Mp->A[0][0] = RCONST(0.0);
    // Mp->c[0] = RCONST(0.0);
    // Mp->q = 1;
    // Mp->p = 0;

    retval = ARKStepSetTables(arkode_mem, 1, 0, Mq, Mp);
    if (check_retval(&retval, "ARKStepSetTables", 1)) return 1;

    retval = ARKStepSetSeparableRhs(arkode_mem, SUNTRUE);
    if (check_retval(&retval, "ARKStepSetSeparableRhs", 1)) return 1;
  } else {
    arkode_mem = ARKStepCreate(dydt, NULL, T0, y, sunctx);
  }

  /* Setup ARKStep */
  retval = ARKStepSetUserData(arkode_mem, (void *) udata);
  if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

  retval = ARKStepSetMaxNumSteps(arkode_mem, ((long int) ceil(Tf/dt)) + 1);
  if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return 1;

  if (fixed_step_mode) {
    retval = ARKStepSetFixedStep(arkode_mem, dt);
    if (check_retval(&retval, "ARKStepSetFixedStep", 1)) return 1;
  } else {
    retval = ARKStepSStolerances(arkode_mem, SUN_RCONST(10e-6), SUN_RCONST(10e-12));
    if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
  }

  tret = T0;
  tout = T0+dTout;
  H0 = Hamiltonian(y);
  L0 = AngularMomentum(y);
  fprintf(stdout, "t = %.2f, H(p, q) = %.16f, L(p, q) = %.16f\n", tret, H0, L0);
  for (iout = 0; iout < num_output_times; iout++) {
    retval = ARKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
    fprintf(stdout, "t = %.2f, H(p, q)-H0 = %.16f, L(p, q)-L0 = %.16f\n", tret, Hamiltonian(y)-H0, AngularMomentum(y)-L0);
    N_VPrint(y);
    if (retval >= 0) {  /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {            /* unsuccessful solve: break */
      fprintf(stderr, "Solver failure, stopping integration\n");
      break;
    }
  }

  ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  // ARKStepFree(arkode_mem); // TODO: figure out why this is segfaulting!
  if (Mq) ARKodeButcherTable_Free(Mq);
  if (Mp) ARKodeButcherTable_Free(Mp);
  N_VDestroy(y);
  free(udata);
  SUNContext_Free(&sunctx);

  return 0;
}

void InitialConditions(N_Vector y0vec, sunrealtype ecc)
{
  const sunrealtype zero = SUN_RCONST(0.0);
  const sunrealtype one  = SUN_RCONST(1.0);
  sunrealtype* y0 = N_VGetArrayPointer(y0vec);

  y0[0] = one - ecc;
  y0[1] = zero;
  y0[2] = zero;
  y0[3] = SUNRsqrt((one + ecc)/(one - ecc));
}

sunrealtype Hamiltonian(N_Vector yvec)
{
  sunrealtype H = 0.0;
  sunrealtype* y = N_VGetArrayPointer(yvec);
  const sunrealtype sqrt_q1q1_plus_q2q2 = SUNRsqrt(y[0]*y[0] + y[1]*y[1]);
  const sunrealtype p1p1_plus_p2p2 = y[2]*y[2] + y[3]*y[3];
  H = SUN_RCONST(0.5)*p1p1_plus_p2p2 - SUN_RCONST(1.0)/sqrt_q1q1_plus_q2q2
    - SUN_RCONST(0.005) / SUN_RCONST(2.0) / SUNRpowerR(sqrt_q1q1_plus_q2q2, SUN_RCONST(3.0));

  return H;
}

sunrealtype AngularMomentum(N_Vector yvec)
{
  sunrealtype L = 0.0;
  sunrealtype* y = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  L = q1*p2 - q2*p1;

  return L;
}

int dydt(realtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  int retval = 0;

  retval += dpdt(t, yvec, ydotvec, user_data);
  retval += dqdt(t, yvec, ydotvec, user_data);

  return retval;
}

int dqdt(realtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y = N_VGetArrayPointer(yvec);
  sunrealtype* ydot = N_VGetArrayPointer(ydotvec);
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  ydot[0] = p1;
  ydot[1] = p2;

  // printf("dqdt(t=%g, p=[%g, %g]) = [%g, %g]\n", t, p1, p2, ydot[0], ydot[1]);

  return 0;
}

int dpdt(realtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  UserData udata = (UserData) user_data;
  sunrealtype* y = N_VGetArrayPointer(yvec);
  sunrealtype* ydot = N_VGetArrayPointer(ydotvec);
  const sunrealtype delta = udata->delta;
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype q1q1_plus_q2q2 = q1*q1 + q2*q2;

  ydot[2] = -q1 / SUNRpowerR(q1q1_plus_q2q2, SUN_RCONST(3.0)/SUN_RCONST(2.0))
            -delta*q1 / SUNRpowerR(q1q1_plus_q2q2, SUN_RCONST(5.0)/SUN_RCONST(2.0));
  ydot[3] = -q2 / SUNRpowerR(q1q1_plus_q2q2, SUN_RCONST(3.0)/SUN_RCONST(2.0))
            -delta*q2 / SUNRpowerR(q1q1_plus_q2q2, SUN_RCONST(5.0)/SUN_RCONST(2.0));

  // printf("dpdt(t=%g, q=[%g, %g]) = [%g, %g]\n", t, q1, q2, ydot[2], ydot[3]);

  return 0;
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}
