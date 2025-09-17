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

#include <nvector/nvector_serial.h>
#include <arkode/arkode_arkstep.h>
#include "arkode/arkode.h"
#include "sundials/sundials_nvector.h"
#include "sundials/sundials_types.h"

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

  void* arkode_mem = ARKStepCreate(ode_rhs, NULL, SUN_RCONST(0.0), y, ctx);
  if (!arkode_mem) return 1;

  ierr = ARKodeSStolerances(arkode_mem, SUN_RCONST(1.0e-4), SUN_RCONST(1.0e-8));
  if (ierr) return 1;

  UserData user_data;
  user_data.A = 1;
  user_data.B = 3;

  ierr = ARKodeSetUserData(arkode_mem, (void*)&user_data);
  if (ierr) return 1;

  SUNAdaptController ctrl = SUNAdaptController_ExpGus(ctx);
  if (!ctrl) return 1;

  ierr = ARKodeSetAdaptController(arkode_mem, ctrl);
  if (ierr) return 1;

  sunrealtype t_return = SUN_RCONST(0.0);
  sunrealtype t_final  = SUN_RCONST(30.0);

  FILE *f_out = fopen("output.txt", "w");

  fprintf(f_out, ""SUN_FORMAT_E ""SUN_FORMAT_E ""SUN_FORMAT_E"\n", t_return,
          y_data[0], y_data[1]);

  while (t_return <= t_final)
  {
    ierr = ARKodeEvolve(arkode_mem, t_final, y, &t_return, ARK_ONE_STEP);
    if (ierr) return 1;

    fprintf(f_out, ""SUN_FORMAT_E ""SUN_FORMAT_E ""SUN_FORMAT_E"\n", t_return,
            y_data[0], y_data[1]);
  }

  ierr = ARKodePrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
  if (ierr) return 1;

  N_VDestroy(y);
  ARKodeFree(&arkode_mem);
  SUNAdaptController_Destroy(ctrl);

  ierr = SUNContext_Free(&ctx);

  return ierr;
}
