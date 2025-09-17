/* -----------------------------------------------------------------------------
 *
 * Brusselator ODE
 *   u' = A + u^2 v - B u - u
 *   v' = B u - u^2 v
 *
 * Constants:
 *   A = 1
 *   B = 3
 *
 * Initial conditions:
 *   u(0) = 1
 *   v(0) = 1
 * ---------------------------------------------------------------------------*/

#include "nvector/nvector_serial.h"
#include "cvode/cvode.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

typedef struct
{
  sunrealtype A;
  sunrealtype B;
} UserData;

static int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* y_data  = N_VGetArrayPointer(y);
  sunrealtype* yd_data = N_VGetArrayPointer(ydot);

  sunrealtype u = y_data[0];
  sunrealtype v = y_data[1];

  UserData* data = (UserData*)user_data;

  yd_data[0] = data->A + u * u * v - data->B * u - u;
  yd_data[1] = data->B * u - u * u * v;

  return 0;
}

int main(int argc, char** argv)
{
  /* Create SUNDIALS context */
  SUNContext ctx;
  int ierr = SUNContext_Create(SUN_COMM_NULL, &ctx);
  if (ierr) return ierr;

  N_Vector y = N_VNew_Serial(2, ctx);
  if (!y) return 1;

  sunrealtype* y_data = N_VGetArrayPointer(y);
  y_data[0] = SUN_RCONST(1.0);
  y_data[1] = SUN_RCONST(1.0);

  void* cvode_mem = CVodeCreate(CV_ADAMS, ctx);
  if (!cvode_mem) return 1;

  ierr = CVodeInit(cvode_mem, ode_rhs, SUN_RCONST(0.0), y);
  if (ierr) return 1;

  ierr = CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (ierr) return 1;

  UserData user_data;
  user_data.A = 1;
  user_data.B = 3;

  ierr = CVodeSetUserData(cvode_mem, (void*)&user_data);
  if (ierr) return 1;

  SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 0, ctx);
  if (!NLS) return 1;

  ierr = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if (ierr) return 1;

  ierr = CVodeSetMaxNonlinIters(cvode_mem, 5);
  if (ierr) return 1;

  sunrealtype t_return = SUN_RCONST(0.0);
  sunrealtype t_final  = SUN_RCONST(30.0);

  FILE *f_out = fopen("output.txt", "w");

  fprintf(f_out, ""SUN_FORMAT_E ""SUN_FORMAT_E ""SUN_FORMAT_E"\n", t_return,
          y_data[0], y_data[1]);

  while (t_return <= t_final)
  {
    ierr = CVode(cvode_mem, t_final, y, &t_return, CV_ONE_STEP);
    if (ierr) return 1;

    fprintf(f_out, ""SUN_FORMAT_E ""SUN_FORMAT_E ""SUN_FORMAT_E"\n", t_return,
            y_data[0], y_data[1]);
  }

  ierr = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (ierr) return 1;

  N_VDestroy(y);
  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);

  ierr = SUNContext_Free(&ctx);

  return ierr;
}
