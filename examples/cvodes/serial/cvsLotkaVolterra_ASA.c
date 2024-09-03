/* -----------------------------------------------------------------------------
  * SUNDIALS Copyright Start
  * Copyright (c) 2002-2024, Lawrence Livermore National Security
  * and Southern Methodist University.
  * All rights reserved.
  *
  * See the top-level LICENSE and NOTICE files for details.
  *
  * SPDX-License-Identifier: BSD-3-Clause
  * SUNDIALS Copyright End
  * -----------------------------------------------------------------------------
  * This example solves the Lotka-Volterra ODE with four parameters,
  *
  *     u = [dx/dt] = [ p_0*x - p_1*x*y  ]
  *         [dy/dt]   [ -p_2*y + p_3*x*y ].
  *
  * The initial condition is u(t_0) = 1.0 and we use the parameters
  * p  = [1.5, 1.0, 3.0, 1.0]. The integration interval is t \in [0, 10.].
  * The implicit BDF method from CVODES is used to solve the forward problem.
  * Afterwards, the continuous adjoint sensitivity analysis capabilites of CVODES
  * are used to obtain the gradient of the cost function,
  *
  *    g(u(t_f), p) = (sum(u)^2) / 2
  *
  * with respect to the initial condition and the parameters.
  * -----------------------------------------------------------------------------
  */

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_core.h>
#include "cvodes/cvodes_ls.h"
#include "sundials/sundials_context.h"
#include "sundials/sundials_iterative.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunlinsol/sunlinsol_spgmr.h"

/* Problem Constants */
#define NEQ   2                   /* number of equations  */
#define T0    SUN_RCONST(0.0)     /* initial time         */
#define TF    SUN_RCONST(10.0)    /* final time           */
#define RTOL  SUN_RCONST(1.0e-10) /* relative tolerance   */
#define ATOL  SUN_RCONST(1.0e-14) /* absolute tolerance  */
#define STEPS 5                   /* checkpoint interval */

static int check_retval(void* retval_ptr, const char* funcname, int opt);

static const sunrealtype params[4] = {1.5, 1.0, 3.0, 1.0};

/* Function to compute the ODE right-hand side */
static int lotka_volterra(sunrealtype t, N_Vector uvec, N_Vector udotvec,
                          void* user_data)
{
  sunrealtype* p    = (sunrealtype*)user_data;
  sunrealtype* u    = N_VGetArrayPointer(uvec);
  sunrealtype* udot = N_VGetArrayPointer(udotvec);

  udot[0] = p[0] * u[0] - p[1] * u[0] * u[1];
  udot[1] = -p[2] * u[1] + p[3] * u[0] * u[1];

  return 0;
}

int vjp(N_Vector vvec, N_Vector Jvvec, sunrealtype t, N_Vector uvec,
        void* user_data)
{
  sunrealtype* p  = (sunrealtype*)user_data;
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* v  = N_VGetArrayPointer(vvec);
  sunrealtype* Jv = N_VGetArrayPointer(Jvvec);

  Jv[0] = (p[0] - p[1] * u[1]) * v[0] + p[3] * u[1] * v[1];
  Jv[1] = -p[1] * u[0] * v[0] + (-p[2] + p[3] * u[0]) * v[1];

  return 0;
}

/* Gradient of the cost function */
void dgdu(N_Vector uvec, N_Vector dgvec)
{
  sunrealtype* u  = N_VGetArrayPointer(uvec);
  sunrealtype* dg = N_VGetArrayPointer(dgvec);

  dg[0] = u[0] + u[1];
  dg[1] = u[0] + u[1];
}

// /* Function to compute the adjoint ODE right-hand side:
//     -lambda*(df/du) - (dg/du)
//  */
// static int adjoint_rhs(sunrealtype t, N_Vector uvec, N_Vector lvec,
//                        N_Vector ldotvec, void* user_data)
// {
//   sunrealtype* p         = (sunrealtype*)user_data;
//   sunrealtype* u         = N_VGetArrayPointer(uvec);
//   sunrealtype* lambda    = N_VGetArrayPointer(lvec);
//   sunrealtype* lambdaDot = N_VGetArrayPointer(ldotvec);

//   N_Vector dg = N_VClone(uvec);
//   vjp(lvec, ldotvec, t, uvec, user_data);
//   dgdu(uvec, dg);
//   N_VLinearSum(-1.0, ldotvec, -1.0, dg, ldotvec);
//   N_VDestroy(dg);

//   return 0;
// }

/* Function to compute the adjoint ODE right-hand side:
    -mu*(df/du)
 */
static int adjoint_rhs(sunrealtype t, N_Vector uvec, N_Vector lvec,
                       N_Vector ldotvec, void* user_data)
{
  sunrealtype* p         = (sunrealtype*)user_data;
  sunrealtype* u         = N_VGetArrayPointer(uvec);
  sunrealtype* lambda    = N_VGetArrayPointer(lvec);
  sunrealtype* lambdaDot = N_VGetArrayPointer(ldotvec);

  vjp(lvec, ldotvec, t, uvec, user_data);
  N_VScale(-1.0, ldotvec, ldotvec);

  return 0;
}

int main(int argc, char* argv[])
{
  SUNContext sunctx;
  sunrealtype reltol, abstol, t, tout;
  N_Vector u, uB;
  void *cvode_mem, *cvode_memB;
  int which, retval;

  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Allocate memory for the solution vector */
  u = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)u, "N_VNew_Serial", 0)) { return 1; }

  /* Initialize the solution vector */
  N_VConst(1.0, u);

  /* Set the tolerances */
  reltol = RTOL;
  abstol = ATOL;

  /* Create the CVODES object */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) { return 1; }

  /* Initialize the CVODES solver */
  retval = CVodeInit(cvode_mem, lotka_volterra, T0, u);
  if (check_retval(&retval, "CVodeInit", 1)) { return 1; }

  /* Set the user data */
  retval = CVodeSetUserData(cvode_mem, (void*)params);
  if (check_retval(&retval, "CVodeSetUserData", 1)) { return 1; }

  /* Set the tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) { return 1; }

  // SUNLinearSolver LS = SUNLinSol_Dense(y, NULL, sunctx);
  SUNLinearSolver LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, 3, sunctx);

  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return 1; }

  retval = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) { return 1; }

  /* Initialize ASA */
  retval = CVodeAdjInit(cvode_mem, STEPS, CV_HERMITE);
  if (check_retval(&retval, "CVodeAdjInit", 1)) { return 1; }

  /* Integrate the ODE */
  tout = TF;
  int ncheck;
  retval = CVodeF(cvode_mem, tout, u, &t, CV_NORMAL, &ncheck);
  if (check_retval(&retval, "CVode", 1)) { return 1; }

  /* Print the final solution */
  printf("Forward Solution at t = %g\n", t);
  N_VPrint(u);

  /* Allocate memory for the adjoint solution vector */
  uB = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)uB, "N_VNew_Serial", 0)) { return 1; }

  /* Initialize the adjoint solution vector */
  dgdu(u, uB);

  printf("Adjoint terminal condition:\n");
  N_VPrint(uB);

  /* Create the CVODES object for the backward problem */
  retval = CVodeCreateB(cvode_mem, CV_BDF, &which);

  /* Initialize the CVODES solver for the backward problem */
  retval = CVodeInitB(cvode_mem, which, adjoint_rhs, TF, uB);
  if (check_retval(&retval, "CVodeInitB", 1)) { return 1; }

  /* Set the user data for the backward problem */
  retval = CVodeSetUserDataB(cvode_mem, which, (void*)params);
  if (check_retval(&retval, "CVodeSetUserDataB", 1)) { return 1; }

  /* Set the tolerances for the backward problem */
  retval = CVodeSStolerancesB(cvode_mem, which, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerancesB", 1)) { return 1; }

  SUNLinearSolver LSB = SUNLinSol_SPGMR(uB, SUN_PREC_NONE, 3, sunctx);

  retval = CVodeSetLinearSolverB(cvode_mem, which, LSB, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) { return 1; }

  /* Integrate the adjoint ODE */
  retval = CVodeB(cvode_mem, T0, CV_NORMAL);
  if (check_retval(&retval, "CVodeB", 1)) { return 1; }

  /* Get the final adjoint solution */
  retval = CVodeGetB(cvode_mem, which, &t, uB);
  if (check_retval(&retval, "CVodeGetB", 1)) { return 1; }

  /* Print the final adjoint solution */
  printf("Adjoint Solution at t = %g:\n", t);
  N_VPrint(uB);

  /* Free memory */
  N_VDestroy(u);
  N_VDestroy(uB);
  SUNLinSolFree(LS);
  SUNLinSolFree(LSB);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  return 0;
}

/* Check function return value */
int check_retval(void* retval_ptr, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && retval_ptr == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)retval_ptr;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  return (0);
}
