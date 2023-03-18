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
 *    H(p,q) = 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2)
 *            - 1/200 / (2 * sqrt(q1^2 + q2^2)^3))
 * is conserved as well as the angular momentum,
 *    L(p,q) = q1*p2 - q2*p1.
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
#include <arkode/arkode_sprkstep.h>      /* prototypes for MRIStep fcts., consts */
#include <arkode/arkode_arkstep.h>      /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>     /* serial N_Vector type, fcts., macros  */
#include <sundials/sundials_math.h>     /* def. math fcns, 'sunrealtype'           */
#include "sundials/sundials_nonlinearsolver.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

static int check_retval(void *returnvalue, const char *funcname, int opt);

static void InitialConditions(N_Vector y0, sunrealtype ecc);
static sunrealtype Hamiltonian(N_Vector y);
static sunrealtype AngularMomentum(N_Vector y);

static int dydt(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);
static int dqdt(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);
static int dpdt(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);

static sunrealtype Q(N_Vector yvec, sunrealtype alpha);
static sunrealtype G(N_Vector yvec, sunrealtype alpha);
static int Adapt(N_Vector y, sunrealtype t, sunrealtype h1, sunrealtype h2,
                 sunrealtype h3, sunrealtype e1, sunrealtype e2,
                 sunrealtype e3, int q, int p, sunrealtype *hnew,
                 void *user_data);

typedef struct {
  sunrealtype ecc;
  sunrealtype delta;

  /* for time-step control */
  sunrealtype eps;
  sunrealtype alpha;
  sunrealtype rho_nmhalf;
  sunrealtype rho_nphalf;
  sunrealtype rho_n;
  sunrealtype rho_np1;

  FILE *hhist_fp;
} *UserData;

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  N_Vector y;
  SUNNonlinearSolver NLS;
  UserData udata;
  sunrealtype tout, tret;
  sunrealtype H0, L0;
  void* arkode_mem;
  FILE *conserved_fp, *solution_fp, *times_fp;
  int argi, iout, retval;

  NLS = NULL;
  y = NULL;

    /* Default problem parameters */
  const sunrealtype T0    = SUN_RCONST(0.0);
  sunrealtype Tf          = SUN_RCONST(1000.0);
  const sunrealtype ecc   = SUN_RCONST(0.6);
  // const sunrealtype delta = SUN_RCONST(0.015); // perturbed
  const sunrealtype delta = SUN_RCONST(0.0); // unperturbed

  /* Default integrator Options */
  int step_mode = 0;
  int method = 0;
  int use_compsums = 0;
  sunrealtype dt = SUN_RCONST(1e-2);
  const sunrealtype dTout = dt;
  // const sunrealtype dTout = SUN_RCONST(100.0);
  const int num_output_times = (int) ceil(Tf/dTout);

  /* Parse CLI args */
  argi = 0;
  if (argc > 1) {
    step_mode = atoi(argv[++argi]);
  }
  if (argc > 2) {
    method = atoi(argv[++argi]);
  }
  if (argc > 3) {
    dt = atof(argv[++argi]);
  }
  if (argc > 4) {
    use_compsums = atoi(argv[++argi]);
  }

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  const int num_orders = 12;
  int orders[num_orders] = { 1, 2, 22, 222, 3, 33, 4, 44, 5, 6, 8, 10 };

  for (int test_idx = 0; test_idx < num_orders; test_idx++) {
    int order = orders[test_idx];

    /* Allocate our state vector */
    y = N_VNew_Serial(4, sunctx);

    /* Allocate and fill udata structure */
    udata = (UserData) malloc(sizeof(*udata));
    udata->ecc   = ecc;
    udata->delta = delta;
    udata->alpha = SUN_RCONST(3.0)/SUN_RCONST(2.0);
    udata->eps   = dt;
    udata->rho_n = SUN_RCONST(1.0);

    /* Fill the initial conditions */
    InitialConditions(y, ecc);
    Tf = SUN_RCONST(1000.0);

    /* Create SPRKStep integrator where we treat dqdt explicitly and dpdt implicitly */
    if (method == 0) {
      arkode_mem = SPRKStepCreate(dpdt, dqdt, T0, y, sunctx);

      retval = SPRKStepSetUseCompSums(arkode_mem, use_compsums);
      if (check_retval(&retval, "SPRKStepSetUseCompSums", 1)) return 1;

      switch (order) {
        case 1:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_EULER_1);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 2:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_LEAPFROG_2);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 22:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_PSEUDO_LEAPFROG_2);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 222:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_MCLACHLAN_2);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 3:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_RUTH_3);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 33:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_MCLACHLAN_3);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 4:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_CANDY_ROZMUS_4);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 44:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_MCLACHLAN_4);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 5:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_MCLACHLAN_5);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 6:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_YOSHIDA_6);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 8:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_MCLACHLAN_8);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        case 10:
          retval = SPRKStepSetMethod(arkode_mem, ARKODE_SYMPLECTIC_SOFRONIOU_10);
          if (check_retval(&retval, "SPRKStepSetMethod", 1)) return 1;
          break;
        default:
          fprintf(stderr, "Not a valid method\n");
          return 1;
      }

      if (step_mode == 0) {
        retval = SPRKStepSetFixedStep(arkode_mem, dt);
        if (check_retval(&retval, "SPRKStepSetFixedStep", 1)) return 1;

        retval = SPRKStepSetMaxNumSteps(arkode_mem, ((long int) ceil(Tf/dt)) + 1);
        if (check_retval(&retval, "SPRKStepSetMaxNumSteps", 1)) return 1;
      } else {
        /*  Adaptivity based on [Hairer and Soderlind, 2005] */
        retval = SPRKStepSetAdaptivityFn(arkode_mem, Adapt, udata);
        if (check_retval(&retval, "SPRKStepSetFixedStep", 1)) return 1;

        udata->rho_nmhalf = udata->rho_n - udata->eps*G(y, udata->alpha)/SUN_RCONST(2.0);
        udata->rho_nphalf = udata->rho_nmhalf + udata->eps*G(y, udata->alpha);
        retval = SPRKStepSetInitStep(arkode_mem, udata->eps/udata->rho_nphalf);
        if (check_retval(&retval, "SPRKStepSetInitStep", 1)) return 1;

        retval = SPRKStepSetMaxNumSteps(arkode_mem, (long int) 100*(ceil(Tf/dt) + 1));
        if (check_retval(&retval, "SPRKStepSetMaxNumSteps", 1)) return 1;
      }

      retval = SPRKStepSetUserData(arkode_mem, (void *) udata);
      if (check_retval(&retval, "SPRKStepSetUserData", 1)) return 1;
    } else if (method >= 1) {
      if (method == 1) {
        arkode_mem = ARKStepCreate(dydt, NULL, T0, y, sunctx);

        retval = ARKStepSetOrder(arkode_mem, order);
        if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;
      } else {
        arkode_mem = ARKStepCreate(dqdt, dpdt, T0, y, sunctx);

        retval = ARKStepSetOrder(arkode_mem, order);
        if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;

        NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
        ARKStepSetNonlinearSolver(arkode_mem, NLS);
      }

      retval = ARKStepSetUserData(arkode_mem, (void *) udata);
      if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

      retval = ARKStepSetMaxNumSteps(arkode_mem, 1000000);
      if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return 1;

      if (step_mode == 0) {
        retval = ARKStepSetFixedStep(arkode_mem, dt);
      } else {
        retval = ARKStepSStolerances(arkode_mem, SUN_RCONST(10e-12), SUN_RCONST(10e-14));
        if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
      }
    }

    /* Open output files */
    if (method == 0) {
      const char* fmt1 = "ark_kepler_conserved_sprk-%d-dt-%.6f.txt";
      const char* fmt2 = "ark_kepler_solution_sprk-%d-dt-%.6f.txt";
      const char* fmt3 = "ark_kepler_times_sprk-%d-dt-%.6f.txt";
      // const char* fmt4 = "ark_kepler_hhist_sprk-%d-dt-%.6f.txt";
      // const char* fmt1 = "ark_kepler_conserved_sprkinc-%d.txt";
      // const char* fmt2 = "ark_kepler_solution_sprkinc-%d.txt";
      // const char* fmt3 = "ark_kepler_times_sprkinc-%d.txt";
      char fname[64];
      sprintf(fname, fmt1, order, dt);
      conserved_fp = fopen(fname, "w+");
      sprintf(fname, fmt2, order, dt);
      solution_fp = fopen(fname, "w+");
      sprintf(fname, fmt3, order, dt);
      times_fp = fopen(fname, "w+");
      // sprintf(fname, fmt4, order, dt);
      // udata->hhist_fp = fopen(fname, "w+");
    } else {
      const char* fmt1 = "ark_kepler_conserved_erk-%d-dt-%.6f.txt";
      const char* fmt2 = "ark_kepler_solution_erk-%d-dt-%.6f.txt";
      const char* fmt3 = "ark_kepler_times_erk-%d-dt-%.6f.txt";
      // const char* fmt4 = "ark_kepler_hhist_erk-%d-dt-%.6f.txt";
      char fname[64];
      sprintf(fname, fmt1, order, dt);
      conserved_fp = fopen(fname, "w+");
      sprintf(fname, fmt2, order, dt);
      solution_fp = fopen(fname, "w+");
      sprintf(fname, fmt3, order, dt);
      times_fp = fopen(fname, "w+");
      // sprintf(fname, fmt4, order, dt);
      // udata->hhist_fp = fopen(fname, "w+");
    }

    printf("\n   Begin Kepler Problem\n\n");
    printf("        order: %d\n", order);

    /* Print out starting energy, momentum before integrating */
    tret = T0;
    tout = T0+dTout;
    H0 = Hamiltonian(y);
    L0 = AngularMomentum(y);
    fprintf(stdout, "t = %.4Lf, H(p,q) = %.16Lf, L(p,q) = %.16Lf\n", tret, H0, L0);
    fprintf(times_fp, "%.16Lf\n", tret);
    fprintf(conserved_fp, "%.16Lf, %.16Lf\n", H0, L0);
    N_VPrintFile(y, solution_fp);

    /* Do integration */
    if (method == 0) {
      for (iout = 0; iout < num_output_times; iout++) {
        // SPRKStepSetStopTime(arkode_mem, tout);
        
        retval = SPRKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);

        /* Output current integration status */
        // fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
        //         tret, Hamiltonian(y)-H0, AngularMomentum(y)-L0);
        fprintf(times_fp, "%.16Lf\n", tret);
        fprintf(conserved_fp, "%.16Lf, %.16Lf\n", Hamiltonian(y), AngularMomentum(y));
        N_VPrintFile(y, solution_fp);
        
        /* Check if the solve was successful, if so, update the time and continue */
        if (retval >= 0) {
          tout += dTout;
          tout = (tout > Tf) ? Tf : tout;
        } else {
          fprintf(stderr, "Solver failure, stopping integration\n");
          break;
        }
      }
    } else {
      for (iout = 0; iout < num_output_times; iout++) {
        ARKStepSetStopTime(arkode_mem, tout);
        if (step_mode == 3) {
          while(tret < tout) {
            retval = ARKStepEvolve(arkode_mem, tout, y, &tret, ARK_ONE_STEP);
            if (retval < 0) break;
            
            /* Output current integration status */
            // fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf, Q(p,q)-Q0 = %.16Lf\n",
            //         tret, Hamiltonian(y)-H0, AngularMomentum(y)-L0,
            //         Q(y, udata->alpha)/udata->rho_np1-Q0);
            fprintf(times_fp, "%.16Lf\n", tret);
            fprintf(conserved_fp, "%.16Lf, %.16Lf\n", Hamiltonian(y), AngularMomentum(y));
            N_VPrintFile(y, solution_fp);
          }
        } else {
          retval = ARKStepEvolve(arkode_mem, tout, y, &tret, ARK_NORMAL);
          
          /* Output current integration status */
          // fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf, Q(p,q)-Q0 = %.16Lf\n",
          //         tret, Hamiltonian(y)-H0, AngularMomentum(y)-L0,
          //         Q(y, udata->alpha)/udata->rho_np1-Q0);
          fprintf(times_fp, "%.16Lf\n", tret);
          fprintf(conserved_fp, "%.16Lf, %.16Lf\n", Hamiltonian(y), AngularMomentum(y));
          N_VPrintFile(y, solution_fp);
        }

        /* Check if the solve was successful, if so, update the time and continue */
        if (retval >= 0) { 
          tout += dTout;
          tout = (tout > Tf) ? Tf : tout;
        } else {
          fprintf(stderr, "Solver failure, stopping integration\n");
          break;
        }
      }
    }

    fprintf(stdout, "t = %.4Lf, H(p,q)-H0 = %.16Lf, L(p,q)-L0 = %.16Lf\n",
        tret, Hamiltonian(y)-H0, AngularMomentum(y)-L0);

    // fclose(udata->hhist_fp);
    free(udata);
    fclose(times_fp);
    fclose(conserved_fp);
    fclose(solution_fp);
    if (NLS) {
      SUNNonlinSolFree(NLS);
    }
    N_VDestroy(y);
    if (method == 0) {
      SPRKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
      SPRKStepFree(&arkode_mem); 
    } else {
      ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
      ARKStepFree(&arkode_mem);  
    }
  }

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
  const sunrealtype sqrt_qTq = SUNRsqrt(y[0]*y[0] + y[1]*y[1]);
  const sunrealtype pTp = y[2]*y[2] + y[3]*y[3];

  // Perturbed
  // H = SUN_RCONST(0.5)*pTp - SUN_RCONST(1.0)/sqrt_qTq
    // - SUN_RCONST(0.005) / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0)) / SUN_RCONST(2.0);

  // Unperturbed
  H = SUN_RCONST(0.5)*pTp - SUN_RCONST(1.0)/sqrt_qTq;

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

int dydt(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  int retval = 0;

  retval += dpdt(t, yvec, ydotvec, user_data);
  retval += dqdt(t, yvec, ydotvec, user_data);

  return retval;
}

int dqdt(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  sunrealtype* y = N_VGetArrayPointer(yvec);
  sunrealtype* ydot = N_VGetArrayPointer(ydotvec);
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  ydot[0] = p1;
  ydot[1] = p2;
  // ydot[2] = ydot[3] = SUN_RCONST(0.0);

  return 0;
}

int dpdt(sunrealtype t, N_Vector yvec, N_Vector ydotvec, void* user_data)
{
  UserData udata = (UserData) user_data;
  sunrealtype* y = N_VGetArrayPointer(yvec);
  sunrealtype* ydot = N_VGetArrayPointer(ydotvec);
  const sunrealtype delta = udata->delta;
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype sqrt_qTq = SUNRsqrt(q1*q1 + q2*q2);

  // ydot[0] = ydot[1] = SUN_RCONST(0.0);
  ydot[2] =         - q1 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0))
            - delta * q1 / SUNRpowerR(sqrt_qTq, SUN_RCONST(5.0));
  ydot[3] =         - q2 / SUNRpowerR(sqrt_qTq, SUN_RCONST(3.0))
            - delta * q2 / SUNRpowerR(sqrt_qTq, SUN_RCONST(5.0));

  return 0;
}

sunrealtype G(N_Vector yvec, sunrealtype alpha)
{
  sunrealtype* y = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];
  const sunrealtype p1 = y[2];
  const sunrealtype p2 = y[3];

  const sunrealtype pTq = p1*q1 + p2*q2;
  const sunrealtype qTq = q1*q1 + q2*q2;

  return (-alpha * pTq / qTq);
}

sunrealtype Q(N_Vector yvec, sunrealtype alpha)
{
  sunrealtype* y = N_VGetArrayPointer(yvec);
  const sunrealtype q1 = y[0];
  const sunrealtype q2 = y[1];

  const sunrealtype qTq = q1*q1 + q2*q2;

  return SUNRpowerR(qTq, -alpha/SUN_RCONST(2.0));
}

int Adapt(N_Vector y, sunrealtype t, sunrealtype h1, sunrealtype h2,
          sunrealtype h3, sunrealtype e1, sunrealtype e2,
          sunrealtype e3, int q, int p, sunrealtype *hnew,
          void *user_data)
{
  UserData udata = (UserData) user_data;

  // fprintf(udata->hhist_fp, "%.16Lf\n", h1);

  const sunrealtype G_np1 = G(y, udata->alpha);
  udata->rho_np1 = udata->rho_nphalf + udata->eps*G_np1/SUN_RCONST(2.0);

  udata->rho_nmhalf = udata->rho_nphalf;
  const sunrealtype rho_nphalf_next = udata->rho_nmhalf + udata->eps*G_np1;

  *hnew = udata->eps/rho_nphalf_next;

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
