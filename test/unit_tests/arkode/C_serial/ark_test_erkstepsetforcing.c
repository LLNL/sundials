/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Unit test for erkStep_SetInnerForcing
 *
 * erkStep_SetInnerForcing is used to provide a polynomial forcing term in
 * ERKStep when it is used as the inner integrator under MRIStep. To check that
 * the forcing is computed and applied correctly we integrate an ODE in time
 * with ERKStep where the RHS consists of only the polynomial forcing term. The
 * solution should be exact when the method order is greater than or equal to
 * the polynomial order.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "arkode/arkode_erkstep.h"
#include "arkode/arkode_erkstep_impl.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_types.h"

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
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

/* Private function to check function return values */
static int check_flag(void* flagvalue, const char* funcname, int opt);

/* Private function to check computed solution */
static int compute_ans(sunrealtype t, sunrealtype tshift, sunrealtype tscale,
                       N_Vector* forcing, int order, N_Vector ans);
static int compute_error(N_Vector y, N_Vector ans, N_Vector tmp, sunrealtype rtol,
                         sunrealtype atol);

/* Main Program */
int main(int argc, char* argv[])
{
  SUNContext sunctx = NULL;

  /* default input values */
  sunindextype NEQ   = 1;               /* number of dependent vars.    */
  int order          = 3;               /* order of polynomial forcing  */
  sunrealtype T0     = SUN_RCONST(0.0); /* initial time                 */
  sunrealtype Tf     = SUN_RCONST(1.0); /* final time                   */
  sunrealtype tshift = T0;              /* time shift for normalization */
  sunrealtype tscale = Tf;              /* time scale for normalization */

  /* tolerances */
  sunrealtype reltol = SUNRsqrt(SUN_UNIT_ROUNDOFF);
  sunrealtype abstol = SUNRsqrt(SUN_UNIT_ROUNDOFF) / 100;

  /* general problem variables */
  int flag;                  /* reusable error-checking flag             */
  N_Vector y        = NULL;  /* vector for storing the computed solution */
  N_Vector ans      = NULL;  /* vector for storing the true solution     */
  N_Vector tmp      = NULL;  /* temporary workspace vector               */
  N_Vector* forcing = NULL;  /* array of forcing vectors                 */
  void* arkode_mem  = NULL;  /* ARKode memory structure                  */
  int i, j;                  /* loop counters                            */
  sunrealtype tret;          /* integrator return time                   */
  sunrealtype* data;         /* array for accessing vector data          */
  long int nst, nst_a;       /* number of integrator steps               */
  long int mxsteps = 100000; /* max steps before output                  */

  /* check inputs */
  if (argc > 1)
  {
    NEQ = (sunindextype)atol(argv[1]);
    if (NEQ <= 0)
    {
      printf("ERROR: The problem size must be a positive integer\n");
      return (1);
    }
  }

  if (argc > 2)
  {
    order = atoi(argv[2]);
    if (order < 0)
    {
      printf("ERROR: The polynomial order must be a non-negative integer\n");
      return (1);
    }
  }

  if (argc > 3)
  {
    if (argc < 5)
    {
      printf("ERROR: Both the initial and final time are required\n");
      return (1);
    }
    T0 = SUNStrToReal(argv[3]);
    Tf = SUNStrToReal(argv[4]);
    if (SUNRabs(T0) >= SUNRabs(Tf))
    {
      printf("ERROR: |T0| must be less than |Tf|\n");
      return (1);
    }
  }

  if (argc > 5)
  {
    if (argc < 7)
    {
      printf("ERROR: Both tshift and tscale are required\n");
      return (1);
    }
    tshift = SUNStrToReal(argv[5]);
    tscale = SUNStrToReal(argv[6]);
    if (SUNRabs(tscale) < TINY)
    {
      printf("ERROR: |tscale| must be greater than %" GSYM "\n", TINY);
      return (1);
    }
  }

  /* Output test setup */
  printf("\nerkStep_SetInnerForcing unit test:\n");
  printf("   NEQ    = %li\n", (long int)NEQ);
  printf("   order  = %i\n", order);
  printf("   t0     = %.1" GSYM "\n", T0);
  printf("   tf     = %.1" GSYM "\n", Tf);
  printf("   tshift = %.1" GSYM "\n", tshift);
  printf("   tscale = %.1" GSYM "\n", tscale);
  printf("   reltol = %.1" ESYM "\n", reltol);
  printf("   abstol = %.1" ESYM "\n\n", abstol);

  /* Create the SUNDIALS context object for this simulation. */
  SUNContext_Create(NULL, &sunctx);

  /* Create solution vector and initialize to zero */
  y = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void*)y, "N_VNew_Serial", 0)) return 1;

  ans = N_VClone(y);
  if (check_flag((void*)ans, "N_VClone", 0)) return 1;

  tmp = N_VClone(y);
  if (check_flag((void*)tmp, "N_VClone", 0)) return 1;

  /* allocate vector array for polynomial forcing */
  forcing = N_VCloneVectorArray(order + 1, y);
  if (check_flag((void*)forcing, "N_VCloneVectorArray", 0)) return 1;
  for (i = 0; i < order + 1; i++)
  {
    if (check_flag((void*)forcing[i], "N_VCloneVectorArray", 0)) return 1;
  }

  /* fill forcing vectors with random data */
  for (i = 0; i < order + 1; i++)
  {
    data = N_VGetArrayPointer(forcing[i]);
    for (j = 0; j < NEQ; j++)
    {
      data[j] = (sunrealtype)rand() / (sunrealtype)RAND_MAX;
    }
  }

  /* compute the true solution */
  flag = compute_ans(Tf, tshift, tscale, forcing, order, ans);
  if (check_flag(&flag, "compute_ans", 1)) return 1;

  printf("True solution:\n");
  N_VPrint_Serial(ans);

  /* ---------------------------------------------------------------------------
   * explicit test
   * -------------------------------------------------------------------------*/

  /* initialize solution vector */
  flag = compute_ans(T0, tshift, tscale, forcing, order, y);
  if (check_flag(&flag, "compute_ans", 1)) return 1;

  /* Create ARKode mem structure */
  arkode_mem = ERKStepCreate(f, T0, y, sunctx);
  if (check_flag((void*)arkode_mem, "ERKStepCreate", 0)) return 1;

  /* Specify tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Set stop time */
  flag = ARKodeSetStopTime(arkode_mem, Tf);
  if (check_flag(&flag, "ARKodeSetStopTime", 1)) return 1;

  /* Set max steps before output */
  flag = ARKodeSetMaxNumSteps(arkode_mem, mxsteps);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Set forcing */
  flag = erkStep_SetInnerForcing(arkode_mem, tshift, tscale, forcing, order + 1);
  if (check_flag(&flag, "erkStep_SetInnerForcing", 1)) return 1;

  /* Integrate the problem */
  flag = ARKodeEvolve(arkode_mem, Tf, y, &tret, ARK_NORMAL);

  /* check for errors */
  if (flag < 0)
  {
    fprintf(stderr, "ARKodeEvolve failure, flag = %d\n", flag);
    return 1;
  }

  /* get some integrator stats */
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);

  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);

  printf("Stats:\n");
  printf("Steps = %li (attempted = %li)\n\n", nst, nst_a);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);
  arkode_mem = NULL;

  /* print solution */
  printf("Explicit solution:\n");
  N_VPrint_Serial(y);

  /* check the solution error */
  flag = compute_error(y, ans, tmp, reltol, abstol);
  if (flag != 0) return (1);

  /* ---------------------------------------------------------------------------
   * Clean up and return
   * -------------------------------------------------------------------------*/

  N_VDestroy(y); /* Free vectors */
  N_VDestroy(ans);
  N_VDestroy(tmp);

  N_VDestroyVectorArray(forcing, order + 1);

  SUNContext_Free(&sunctx);
  return flag;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* ODE RHS function */
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);
  return 0;
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  /* Check if flag < 0 */
  else if (opt == 1)
  {
    errflag = (int*)flagvalue;
    if (*errflag < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1;
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

/* computed the true solution at a given time */
static int compute_ans(sunrealtype t, sunrealtype tshift, sunrealtype tscale,
                       N_Vector* forcing, int order, N_Vector ans)
{
  int i;
  sunrealtype tau;
  N_Vector* vecs;
  sunrealtype* vals;

  vals = NULL;
  vals = (sunrealtype*)calloc(order + 1, sizeof(sunrealtype));
  if (vals == NULL) return (1);

  vecs = NULL;
  vecs = (N_Vector*)calloc(order + 1, sizeof(N_Vector));
  if (vecs == NULL) return (1);

  /* compute normalized time */
  tau = (t - tshift) / tscale;

  /* compute true solution */
  for (i = 0; i < order + 1; i++)
  {
    vals[i] = ((SUNRpowerI(tau, i + 1)) / (i + 1)) * tscale;
    vecs[i] = forcing[i];
  }
  N_VLinearCombination(order + 1, vals, vecs, ans);

  free(vals);
  free(vecs);

  return (0);
}

/* compure the weighted max norm of the difference of two vectors */
static int compute_error(N_Vector y, N_Vector ans, N_Vector tmp, sunrealtype rtol,
                         sunrealtype atol)
{
  int status; /* success (0) or failure (1) flag */
  sunrealtype error;

  /* compute the error in y */
  N_VLinearSum(ONE, y, -ONE, ans, y);

  /* compute error weights */
  N_VAbs(ans, tmp);
  N_VScale(rtol, tmp, tmp);
  N_VAddConst(tmp, atol, tmp);
  N_VInv(tmp, tmp);

  /* compute weighted max norm */
  N_VProd(tmp, y, y);
  error = N_VMaxNorm(y);

  /* is the solution within the tolerances? */
  status = (error < ONE) ? 0 : 1;

  if (status)
    fprintf(stdout, "ERROR: Test failed with wmax error = %" GSYM "\n\n", error);
  else
    fprintf(stdout, "SUCCESS: Test passed with wmax error = %" GSYM "\n\n",
            error);

  return (status);
}

/*---- end of file ----*/
