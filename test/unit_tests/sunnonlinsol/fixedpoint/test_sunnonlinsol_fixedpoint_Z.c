/* -----------------------------------------------------------------------------
 * Programmer(s): Sylvia Amihere @ SMU, Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves the complex-valued nonlinear system
 *
 * 4x        - sin(y) - zi     - 1  = 0; 
 * -x^2      + 5y     - cos(z) - 2i = 0;
 * - exp(-x) -y       +6z      - 3  = 0;
 *
 * using the accelerated fixed pointer solver in KINSOL. The nonlinear fixed
 * point function is
 *
 * g1(x,y,z) = (1/4) (sin(y)  + zi     + 1)
 * g2(x,y,z) = (1/5) (x^2     + cos(z) + 2i)
 * g3(x,y,z) = (1/6) (exp(-x) + y      + 3)
 *
 * This system has the analytic solution: x = 0.28443101049565 + 0.27031686078054i
 *                                        y = 0.16117132843381 + 0.42622240595676i
 *                                        z = 0.64771494226506 + 0.03754877135588i.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_types.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

/* precision specific formatting macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

/* precision specific math function macros */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsin(x) (sin((x)))
#define SUNRcos(x) (cos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsin(x) (sinf((x)))
#define SUNRcos(x) (cosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsin(x) (sinl((x)))
#define SUNRcos(x) (cosl((x)))
#endif

/* problem constants */
#define NEQ 3 /* number of equations */

#define ZERO         SUN_RCONST(0.0)             /* real 0.0  */
#define PTONE        SUN_RCONST(0.1)             /* real 0.1  */
#define HALF         SUN_RCONST(0.5)             /* real 0.5  */
#define PTNINE       SUN_RCONST(0.9)             /* real 0.9  */
#define ONE          SUN_RCONST(1.0)             /* real 1.0  */
#define ONEPTZEROSIX SUN_RCONST(1.06)            /* real 1.06 */
#define THREE        SUN_RCONST(3.0)             /* real 3.0  */
#define SIX          SUN_RCONST(6.0)             /* real 6.0  */
#define NINE         SUN_RCONST(9.0)             /* real 9.0  */
#define TEN          SUN_RCONST(10.0)            /* real 10.0 */
#define TWENTY       SUN_RCONST(20.0)            /* real 20.0 */
#define SIXTY        SUN_RCONST(60.0)            /* real 60.0 */
#define PI           SUN_RCONST(3.1415926535898) /* real pi   */

/* analytic solution */
#define XTRUE HALF
#define YTRUE ONE
#define ZTRUE -PI / SIX

/* Check the system solution */
static int check_ans(N_Vector ycur, sunrealtype tol);

/* Check function return values */
static int check_retval(void* flagvalue, const char* funcname, int opt);

/* Nonlinear fixed point function */
static int FPFunction(N_Vector y, N_Vector f, void* mem);

/* Convergence test function */
static int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                    sunrealtype tol, N_Vector ewt, void* mem);

/*
 * Proxy for integrator memory struct
 */

/* Integrator memory structure */
typedef struct IntegratorMemRec
{
  N_Vector y0;
  N_Vector ycor;
  N_Vector ycur;
  N_Vector w;
}* IntegratorMem;

/* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  IntegratorMem Imem     = NULL;
  int retval             = 0;
  SUNNonlinearSolver NLS = NULL;
  sunrealtype tol        = 100 * SUNRsqrt(SUN_UNIT_ROUNDOFF);
  int mxiter             = 20;
  int maa                = 0;               /* no acceleration */
  sunrealtype damping    = SUN_RCONST(1.0); /* no damping      */
  long int niters        = 0;
  sunscalartype* data    = NULL;
  SUNContext sunctx      = NULL;

  /* Check if a acceleration/damping values were provided */
  if (argc > 1) { maa = atoi(argv[1]); }
  if (argc > 2) { damping = (sunrealtype)atof(argv[2]); }

  /* Print problem description */
  printf("Solve the nonlinear system:\n");
  printf("    4x - sin(y) - zi - 1  = 0\n");
  printf("    -x^2 + 5y - cos(z) - 2i = 0\n");
  printf("    - exp(-x) -y +6z - 3 = 0\n");
  printf("Analytic solution:\n");
  printf("    x = %f + %fI\n", creal(XTRUE), cimag(XTRUE));
  printf("    y = %f + %fI\n", creal(YTRUE), cimag(YTRUE));
  printf("    z = %f + %fI\n", creal(ZTRUE), cimag(ZTRUE));
  printf("Solution method: Anderson accelerated fixed point iteration.\n");
  printf("    tolerance = %" GSYM "\n", tol);
  printf("    max iters = %d\n", mxiter);
  printf("    accel vec = %d\n", maa);
  printf("    damping   = %" GSYM "\n", damping);

  /* create SUNDIALS context */
  retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) { return (1); }

  /* create proxy for integrator memory */
  Imem = (IntegratorMem)malloc(sizeof(struct IntegratorMemRec));
  if (check_retval((void*)Imem, "Creating Integrator Memory", 0))
  {
    return (1);
  }

  /* create vectors */
  Imem->y0 = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void*)Imem->y0, "N_VNew_Serial", 0)) { return (1); }

  Imem->ycor = N_VClone(Imem->y0);
  if (check_retval((void*)Imem->ycor, "N_VClone", 0)) { return (1); }

  Imem->ycur = N_VClone(Imem->y0);
  if (check_retval((void*)Imem->ycur, "N_VClone", 0)) { return (1); }

  Imem->w = N_VClone(Imem->y0);
  if (check_retval((void*)Imem->w, "N_VClone", 0)) { return (1); }

  /* set initial guess */
  data = N_VGetArrayPointer(Imem->y0);
  if (check_retval((void*)data, "N_VGetArrayPointer", 0)) { return (1); }

  data[0] = PTONE;
  data[1] = PTONE;
  data[2] = -PTONE;

  /* set initial correction */
  N_VConst(ZERO, Imem->ycor);

  /* set weights */
  N_VConst(ONE, Imem->w);

  /* create nonlinear solver */
  NLS = SUNNonlinSol_FixedPoint(Imem->y0, maa, sunctx);
  if (check_retval((void*)NLS, "SUNNonlinSol_FixedPoint", 0)) { return (1); }

  /* set the nonlinear residual function */
  retval = SUNNonlinSolSetSysFn(NLS, FPFunction);
  if (check_retval(&retval, "SUNNonlinSolSetSysFn", 1)) { return (1); }

  /* set the convergence test function */
  retval = SUNNonlinSolSetConvTestFn(NLS, ConvTest, NULL);
  if (check_retval(&retval, "SUNNonlinSolSetConvTestFn", 1)) { return (1); }

  /* set the maximum number of nonlinear iterations */
  retval = SUNNonlinSolSetMaxIters(NLS, mxiter);
  if (check_retval(&retval, "SUNNonlinSolSetMaxIters", 1)) { return (1); }

  /* set the damping parameter */
  retval = SUNNonlinSolSetDamping_FixedPoint(NLS, damping);
  if (check_retval(&retval, "SUNNonlinSolSetDamping", 1)) { return (1); }

  /* solve the nonlinear system */
  retval = SUNNonlinSolSolve(NLS, Imem->y0, Imem->ycor, Imem->w, tol, SUNTRUE,
                             Imem);
  if (check_retval(&retval, "SUNNonlinSolSolve", 1)) { return (1); }

  /* update the initial guess with the final correction */
  N_VLinearSum(ONE, Imem->y0, ONE, Imem->ycor, Imem->ycur);

  /* get the number of linear iterations */
  retval = SUNNonlinSolGetNumIters(NLS, &niters);
  if (check_retval(&retval, "SUNNonlinSolGetNumIters", 1)) { return (1); }

  printf("Number of nonlinear iterations: %ld\n", niters);

  /* check solution */
  retval = check_ans(Imem->ycur, tol);

  /* Free vector, matrix, linear solver, and nonlinear solver */
  N_VDestroy(Imem->y0);
  N_VDestroy(Imem->ycor);
  N_VDestroy(Imem->ycur);
  N_VDestroy(Imem->w);
  SUNNonlinSolFree(NLS);
  free(Imem);
  SUNContext_Free(&sunctx);

  return (retval);
}

/* Proxy for integrator convergence test function */
int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del, sunrealtype tol,
             N_Vector ewt, void* mem)
{
  sunrealtype delnrm;

  /* compute the norm of the correction */
  delnrm = N_VMaxNorm(del);

  if (delnrm <= tol) { return (SUN_SUCCESS); /* success       */ }
  else { return (SUN_NLS_CONTINUE); /* not converged */ }
}

/* -----------------------------------------------------------------------------
 * Complex-valued Nonlinear system
 *
 * 4x       - sin(y) - zi     - 1  = 0; 
 * -x^2      + 5y     - cos(z) - 2i = 0;
 * -exp(-x)  -y       + 6z     - 3  = 0;
 *
 * Nonlinear fixed point function
 *
 * g1(x,y,z) = (1/4) (sin(y)  + zi     + 1)
 * g2(x,y,z) = (1/5) (x^2     + cos(z) + 2i)
 * g3(x,y,z) = (1/6) (exp(-x) + y      + 3)
 * ---------------------------------------------------------------------------*/
int FPFunction(N_Vector ycor, N_Vector gvec, void* mem)
{
  IntegratorMem Imem;
  sunscalartype* ydata = NULL;
  sunscalartype* gdata = NULL;
  sunrealtype x, y, z;

  if (mem == NULL)
  {
    printf("ERROR: Integrator memory is NULL");
    return (-1);
  }
  Imem = (IntegratorMem)mem;

  /* update state based on current correction */
  N_VLinearSum(ONE, Imem->y0, ONE, ycor, Imem->ycur);

  /* Get vector data arrays */
  ydata = N_VGetArrayPointer(Imem->ycur);
  if (check_retval((void*)ydata, "N_VGetArrayPointer", 0)) { return (-1); }

  gdata = N_VGetArrayPointer(gvec);
  if (check_retval((void*)gdata, "N_VGetArrayPointer", 0)) { return (-1); }

  /* get vector components */
  x = ydata[0];
  y = ydata[1];
  z = ydata[2];

  /* compute fixed point function */
  gdata[0] = (1.0 / 4.0) * (SIN(y) + SUN_CCONST(0.0, 1.0) * z + 1.0);
  gdata[1] = (1.0 / 5.0) * (x * x + COS(z) + SUN_CCONST(0.0, 2.0));
  gdata[2] = (1.0 / 6.0) * (EXP(-x) + y + 3.0);

  N_VLinearSum(ONE, gvec, -ONE, Imem->y0, gvec);

  return (0);
}

/* -----------------------------------------------------------------------------
 * Check the solution of the nonlinear system and return PASS or FAIL
 * ---------------------------------------------------------------------------*/
static int check_ans(N_Vector u, sunrealtype tol)
{
  sunscalartype* data = NULL;
  // sunrealtype ex, ey, ez;
  sunrealtype exR, eyR, ezR;
  sunrealtype exI, eyI, ezI;

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void*)data, "N_VGetArrayPointer", 0)) { return (1); }

  /* print the solution */
  printf("Computed solution:\n");
  printf("    x = %f + %fI\n", creal(data[0]), cimag(data[0]));
  printf("    y = %f + %fI\n", creal(data[1]), cimag(data[1]));
  printf("    z = %f + %fI\n", creal(data[2]), cimag(data[2]));

  /* solution error */
  exR = ABS(creal(data[0]) - creal(XTRUE));
  eyR = ABS(creal(data[1]) - creal(YTRUE));
  ezR = ABS(creal(data[2]) - creal(ZTRUE));

  exI = ABS(cimag(data[0]) - cimag(XTRUE));
  eyI = ABS(cimag(data[1]) - cimag(YTRUE));
  ezI = ABS(cimag(data[2]) - cimag(ZTRUE));

  // /* print the solution error */
  printf("Solution error:\n");
  printf("    ex = %f + %fI\n", exR, exI);
  printf("    ey = %f + %fI\n", eyR, eyI);
  printf("    ez = %f + %fI\n", ezR, ezI);

  tol *= TEN;
  if (exR > tol && exI > tol || eyR > tol && eyI > tol || ezR > tol && ezI > tol)
  {
    printf("FAIL\n");
    return (1);
  }

  printf("PASS\n");
  return (0);
}

/* -----------------------------------------------------------------------------
 * Check function return value
 *   opt == 0 check if returned NULL pointer
 *   opt == 1 check if returned a non-zero value
 * ---------------------------------------------------------------------------*/
static int check_retval(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if the function returned a NULL pointer -- no memory allocated */
  if (opt == 0)
  {
    if (flagvalue == NULL)
    {
      fprintf(stderr, "\nERROR: %s() failed -- returned NULL\n\n", funcname);
      return (1);
    }
    else { return (0); }
  }

  /* Check if the function returned an non-zero value -- internal failure */
  if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag != 0)
    {
      fprintf(stderr, "\nERROR: %s() failed -- returned %d\n\n", funcname,
              *errflag);
      return (1);
    }
    else { return (0); }
  }

  /* if we make it here then opt was not 0 or 1 */
  fprintf(stderr, "\nERROR: check_retval failed -- Invalid opt value\n\n");
  return (1);
}
