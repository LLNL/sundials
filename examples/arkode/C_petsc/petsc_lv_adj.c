static char help[] =
  "Performs adjoint sensitivity analysis for the Lotka Volterra equation.\n";

#include <petsctao.h>
#include <petscts.h>

typedef struct _n_User* User;

static const PetscReal p[4] = {1.5, 1.0, 3.0, 1.0};

struct _n_User
{
  PetscReal mu;
  PetscReal next_output;
  PetscBool imex;
  /* Sensitivity analysis support */
  PetscInt steps;
  PetscReal ftime;
  Mat A;                    /* IJacobian matrix */
  Mat B;                    /* RHSJacobian matrix */
  Mat Jacp;                 /* IJacobianP matrix */
  Mat Jacprhs;              /* RHSJacobianP matrix */
  Vec U, lambda[2], mup[2]; /* adjoint variables */
};

/* ----------------------- Explicit form of the ODE  -------------------- */

static PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec U, Vec F, void* ctx)
{
  User user = (User)ctx;
  PetscScalar* f;
  const PetscScalar* u;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  PetscCall(VecGetArray(F, &f));
  f[0] = p[0] * u[0] - p[1] * u[0] * u[1];
  f[1] = -p[2] * u[1] + p[3] * u[0] * u[1];
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscCall(VecRestoreArray(F, &f));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RHSJacobian(TS ts, PetscReal t, Vec U, Mat A, Mat B,
                                  void* ctx)
{
  User user         = (User)ctx;
  PetscReal mu      = user->mu;
  PetscInt rowcol[] = {0, 1};
  PetscScalar J[2][2];
  const PetscScalar* u;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  J[0][0] = p[0] - p[1] * u[1];
  J[0][1] = -p[1] * u[0];
  J[1][0] = p[3] * u[1];
  J[1][1] = p[3] * u[0] - p[2];
  PetscCall(MatSetValues(A, 2, rowcol, 2, rowcol, &J[0][0], INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  if (A != B)
  {
    PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
  }
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RHSJacobianP(TS ts, PetscReal t, Vec U, Mat A, void* ctx)
{
  User user = (User)ctx;

  PetscFunctionBeginUser;
  PetscInt row[] = {0, 1}, col[] = {0};
  PetscScalar J[2][1];
  const PetscScalar* u;
  PetscCall(VecGetArrayRead(U, &u));
  J[0][0] = 0;
  J[1][0] = (1. - u[0] * u[0]) * u[1] - u[0];
  PetscCall(MatSetValues(A, 2, row, 1, col, &J[0][0], INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscReal g(Vec u, const PetscReal* p, PetscReal t)
{
  /* (sum(u) .^ 2) ./ 2 */
  const PetscScalar* uarr;
  PetscInt vec_len;
  VecGetArrayRead(u, &uarr);
  VecGetLocalSize(u, &vec_len);
  PetscReal sum = 0.0;
  for (PetscInt i = 0; i < vec_len; i++) { sum += uarr[i]; }
  VecRestoreArrayRead(u, &uarr);
  return (sum * sum) / 2.0;
}

static void dgdu(Vec uvec, Vec dgvec, const PetscReal* p, PetscReal t)
{
  const PetscScalar* u;
  PetscScalar* dg;
  VecGetArrayRead(uvec, &u);
  VecGetArrayWrite(dgvec, &dg);
  dg[0] = u[0] + u[1];
  dg[1] = u[0] + u[1];
  VecRestoreArrayRead(uvec, &u);
  VecRestoreArrayWrite(dgvec, &dg);
}

/* Monitor timesteps and use interpolation to output at integer multiples of 0.1 */
static PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec U, void* ctx)
{
  const PetscScalar* u;
  PetscReal tfinal, dt;
  User user = (User)ctx;
  Vec interpolatedU;

  PetscFunctionBeginUser;
  PetscCall(TSGetTimeStep(ts, &dt));
  PetscCall(TSGetMaxTime(ts, &tfinal));

  while (user->next_output <= t && user->next_output <= tfinal)
  {
    PetscCall(VecDuplicate(U, &interpolatedU));
    PetscCall(TSInterpolate(ts, user->next_output, interpolatedU));
    PetscCall(VecGetArrayRead(interpolatedU, &u));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "[%g] %" PetscInt_FMT " TS %g (dt = %g) X %g %g\n",
                          (double)user->next_output, step, (double)t,
                          (double)dt, (double)PetscRealPart(u[0]),
                          (double)PetscRealPart(u[1])));
    PetscCall(VecRestoreArrayRead(interpolatedU, &u));
    PetscCall(VecDestroy(&interpolatedU));
    user->next_output += 0.1;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char** argv)
{
  TS ts;
  PetscBool monitor = PETSC_FALSE, implicitform = PETSC_FALSE;
  PetscScalar *x_ptr, *y_ptr, derp;
  PetscMPIInt size;
  struct _n_User user;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE,
             "This is a uniprocessor example only!");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set runtime options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user.next_output = 0.0;
  user.mu          = 1.0e3;
  user.steps       = 0;
  user.ftime       = 1.0;
  user.imex        = PETSC_FALSE;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-monitor", &monitor, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-mu", &user.mu, NULL));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create necessary matrix and vectors, solve same ODE on every process
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &user.A));
  PetscCall(MatSetSizes(user.A, PETSC_DECIDE, PETSC_DECIDE, 2, 2));
  PetscCall(MatSetFromOptions(user.A));
  PetscCall(MatSetUp(user.A));
  PetscCall(MatCreateVecs(user.A, &user.U, NULL));
  PetscCall(MatDuplicate(user.A, MAT_DO_NOT_COPY_VALUES, &user.B));
  PetscCall(MatCreate(PETSC_COMM_WORLD, &user.Jacp));
  PetscCall(MatSetSizes(user.Jacp, PETSC_DECIDE, PETSC_DECIDE, 2, 1));
  PetscCall(MatSetFromOptions(user.Jacp));
  PetscCall(MatSetUp(user.Jacp));
  PetscCall(MatDuplicate(user.Jacp, MAT_DO_NOT_COPY_VALUES, &user.Jacprhs));
  PetscCall(MatZeroEntries(user.Jacprhs));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(
    TSSetEquationType(ts, TS_EQ_ODE_EXPLICIT)); /* less Jacobian evaluations when adjoint BEuler is used, otherwise no effect */
  PetscCall(TSSetRHSFunction(ts, NULL, RHSFunction, &user));
  PetscCall(TSSetRHSJacobian(ts, user.A, user.A, RHSJacobian, &user));
  PetscCall(TSSetRHSJacobianP(ts, user.Jacp, RHSJacobianP, &user));
  PetscCall(TSSetType(ts, TSRK));
  PetscCall(TSSetMaxTime(ts, user.ftime));
  PetscCall(TSSetTimeStep(ts, 0.001));
  PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
  if (monitor) { PetscCall(TSMonitorSet(ts, Monitor, &user, NULL)); }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGetArray(user.U, &x_ptr));
  x_ptr[0] = 1.0;
  x_ptr[1] = 1.0;
  PetscCall(TSSetTimeStep(ts, 0.0001));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Save trajectory of solution so that TSAdjointSolve() may be used
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TSSetSaveTrajectory(ts));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TSSetFromOptions(ts));

  PetscCall(TSSolve(ts, user.U));
  PetscCall(TSGetSolveTime(ts, &user.ftime));
  PetscCall(TSGetStepNumber(ts, &user.steps));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Adjoint model starts here
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* Set initial conditions for the adjoint integration */

  PetscCall(MatCreateVecs(user.A, &user.lambda[0], NULL));
  dgdu(user.U, user.lambda[0], p, user.ftime);

  PetscCall(MatCreateVecs(user.Jacp, &user.mup[0], NULL));
  PetscCall(VecGetArray(user.mup[0], &x_ptr));
  x_ptr[0] = 0.0;
  PetscCall(VecRestoreArray(user.mup[0], &x_ptr));

  PetscCall(TSSetCostGradients(ts, 1, user.lambda, user.mup));

  /* Integrate adjoint system */
  PetscCall(TSAdjointSolve(ts));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                        "\n sensitivity wrt initial conditions: "
                        "d[y(tf)]/d[y0]  d[y(tf)]/d[z0]\n"));
  PetscCall(VecView(user.lambda[0], PETSC_VIEWER_STDOUT_WORLD));

  // PetscCall(VecGetArray(user.mup[0], &x_ptr));
  // PetscCall(VecGetArray(user.lambda[0], &y_ptr));
  // derp = y_ptr[1] * (-10.0 / (81.0 * user.mu * user.mu) +
  //                    2.0 * 292.0 / (2187.0 * user.mu * user.mu * user.mu)) +
  //        x_ptr[0];
  // PetscCall(VecRestoreArray(user.mup[0], &x_ptr));
  // PetscCall(VecRestoreArray(user.lambda[0], &y_ptr));
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD,
  //                       "\n sensitivity wrt parameters: d[y(tf)]/d[mu]\n%g\n",
  //                       (double)PetscRealPart(derp)));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatDestroy(&user.A));
  PetscCall(MatDestroy(&user.B));
  PetscCall(MatDestroy(&user.Jacp));
  PetscCall(MatDestroy(&user.Jacprhs));
  PetscCall(VecDestroy(&user.U));
  PetscCall(VecDestroy(&user.lambda[0]));
  PetscCall(VecDestroy(&user.lambda[1]));
  PetscCall(VecDestroy(&user.mup[0]));
  PetscCall(VecDestroy(&user.mup[1]));
  PetscCall(TSDestroy(&ts));

  PetscCall(PetscFinalize());
  return 0;
}
