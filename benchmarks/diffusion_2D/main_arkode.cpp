/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * ARKODE main for 2D diffusion benchmark problem
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"
#include "arkode/arkode_arkstep.h"

struct UserOptions
{
  // Integrator settings
  realtype rtol        = RCONST(1.0e-5);   // relative tolerance
  realtype atol        = RCONST(1.0e-10);  // absolute tolerance
  realtype hfixed      = ZERO;             // fixed step size
  int      order       = 3;                // ARKode method order
  int      controller  = 0;                // step size adaptivity method
  int      maxsteps    = 0;                // max steps between outputs
  int      onestep     = 0;                // one step mode, number of steps
  bool     linear      = true;             // linearly implicit RHS

  // Linear solver and preconditioner settings
  std::string ls              = "cg";   // linear solver to use
  bool        preconditioning = true;   // preconditioner on/off
  bool        lsinfo          = false;  // output residual history
  int         liniters        = 20;     // number of linear iterations
  int         msbp            = 0;      // preconditioner setup frequency
  realtype    epslin          = ZERO;   // linear solver tolerance factor

  // Helper functions
  int parse_args(vector<string> &args, bool outproc);
  void help();
  void print();
};


// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // Reusable error-checking flag
  int flag;

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;

  // Create SUNDIALS context
  MPI_Comm    comm = MPI_COMM_WORLD;
  SUNContext  ctx  = NULL;
  SUNProfiler prof = NULL;

  flag = SUNContext_Create((void*) &comm, &ctx);
  if (check_flag(&flag, "SUNContextCreate", 1)) return 1;

  flag = SUNContext_GetProfiler(ctx, &prof);
  if (check_flag(&flag, "SUNContext_GetProfiler", 1)) return 1;

  // Add scope so objects are destroyed before MPI_Finalize
  {
    SUNDIALS_CXX_MARK_FUNCTION(prof);

    // MPI process ID
    int myid;
    flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

    bool outproc = (myid == 0);

    // --------------------------
    // Parse command line inputs
    // --------------------------

    UserData    udata(prof);
    UserOptions uopts;
    UserOutput  uout;

    vector<string> args(argv + 1, argv + argc);

    flag = udata.parse_args(args, outproc);
    if (check_flag(&flag, "UserData::parse_args", 1)) return 1;

    flag = uopts.parse_args(args, outproc);
    if (check_flag(&flag, "UserOptions::parse_args", 1)) return 1;

    flag = uout.parse_args(args, outproc);
    if (check_flag(&flag, "UserOutput::parse_args", 1)) return 1;

    // Check for unparsed inputs
    if (args.size() > 0)
    {
      if (find(args.begin(), args.end(), "--help") == args.end())
      {
        cerr << "ERROR: Unknown inputs: ";
        for (auto i = args.begin(); i != args.end(); ++i)
          cerr << *i << ' ';
        cerr << endl;
      }
      return 1;
    }

    // -----------------------------
    // Setup parallel decomposition
    // -----------------------------

    flag = udata.setup();
    if (check_flag(&flag, "UserData::setup", 1)) return 1;

    if (outproc)
    {
      udata.print();
      uopts.print();
      uout.print();
    }

    // ---------------
    // Create vectors
    // ---------------

    // Create vector for solution
#if defined(USE_HIP)
    N_Vector u = N_VMake_MPIPlusX(udata.comm_c,
                                  N_VNew_Hip(udata.nodes_loc, ctx), ctx);
    if (check_flag((void *) u, "N_VMake_MPIPlusX", 0)) return 1;
#elif defined(USE_CUDA)
    N_Vector u = N_VMake_MPIPlusX(udata.comm_c,
                                  N_VNew_Cuda(udata.nodes_loc, ctx), ctx);
    if (check_flag((void *) u, "N_VMake_MPIPlusX", 0)) return 1;
#else
    N_Vector u = N_VNew_Parallel(udata.comm_c, udata.nodes_loc, udata.nodes,
                                 ctx);
    if (check_flag((void *) u, "N_VNew_Parallel", 0)) return 1;
#endif

    // Set initial condition
    flag = Solution(ZERO, u, &udata);
    if (check_flag(&flag, "Solution", 1)) return 1;

    // Create vector for error
    if (udata.forcing)
    {
      uout.error = N_VClone(u);
      if (check_flag((void *) (uout.error), "N_VClone", 0)) return 1;
    }

    // ---------------------
    // Create linear solver
    // ---------------------

    // Create linear solver
    SUNLinearSolver LS = NULL;
    SUNMatrix A = nullptr;
#if defined(USE_SUPERLU_DIST)
    // SuperLU-DIST objects
    SuperMatrix A_super;
    gridinfo_t grid;
    dLUstruct_t A_lu;
    dScalePermstruct_t A_scaleperm;
    dSOLVEstruct_t A_solve;
    SuperLUStat_t A_stat;
    superlu_dist_options_t A_opts;
    sunrealtype* A_data = nullptr;
    sunindextype* A_col_idxs = nullptr;
    sunindextype* A_row_ptrs = nullptr;
#endif

    int prectype = (uopts.preconditioning) ? SUN_PREC_RIGHT : SUN_PREC_NONE;

    if (uopts.ls == "cg")
    {
      LS = SUNLinSol_PCG(u, prectype, uopts.liniters, ctx);
      if (check_flag((void *) LS, "SUNLinSol_PCG", 0)) return 1;

    }
    else if (uopts.ls == "gmres")
    {
      LS = SUNLinSol_SPGMR(u, prectype, uopts.liniters, ctx);
      if (check_flag((void *) LS, "SUNLinSol_SPGMR", 0)) return 1;
    }
    else
    {
#if defined(USE_SUPERLU_DIST)
      // Initialize SuperLU-DIST grid
      superlu_gridinit(udata.comm_c, udata.npx, udata.npy, &grid);

      // Create arrays for CSR matrix: data, column indices, and row pointers
      sunindextype nnz_loc = 5 * udata.nodes_loc;

      A_data = (realtype*)malloc(nnz_loc * sizeof(sunrealtype));
      if (check_flag((void*)A_data, "malloc Adata", 0)) return 1;

      A_col_idxs = (sunindextype*)malloc(nnz_loc * sizeof(sunindextype));
      if (check_flag((void*)A_col_idxs, "malloc Acolind", 0)) return 1;

      A_row_ptrs = (sunindextype*)malloc((udata.nodes_loc + 1) * sizeof(sunindextype));
      if (check_flag((void*)A_row_ptrs, "malloc Arowptr", 0)) return 1;

      // Create and initialize SuperLU_DIST structures
      dCreate_CompRowLoc_Matrix_dist(&A_super, udata.nodes, udata.nodes,
                                     nnz_loc, udata.nodes_loc, 0, A_data,
                                     A_col_idxs, A_row_ptrs, SLU_NR_loc, SLU_D,
                                     SLU_GE);
      dScalePermstructInit(udata.nodes, udata.nodes, &A_scaleperm);
      dLUstructInit(udata.nodes, &A_lu);
      PStatInit(&A_stat);
      set_default_options_dist(&A_opts);
      A_opts.PrintStat = NO;

      // SUNDIALS structures
      A = SUNMatrix_SLUNRloc(&A_super, &grid, ctx);
      if (check_flag((void*)A, "SUNMatrix_SLUNRloc", 0)) return 1;

      LS = SUNLinSol_SuperLUDIST(u, A, &grid, &A_lu, &A_scaleperm, &A_solve,
                                 &A_stat, &A_opts, ctx);
      if (check_flag((void*)LS, "SUNLinSol_SuperLUDIST", 0)) return 1;

      uopts.preconditioning = false;
#else
      std::cerr << "ERROR: Benchmark was not built with SuperLU_DIST enabled\n";
      return 1;
#endif
    }

    // Allocate preconditioner workspace
    if (uopts.preconditioning)
    {
      udata.diag = N_VClone(u);
      if (check_flag((void *) (udata.diag), "N_VClone", 0)) return 1;
    }

    // --------------
    // Setup ARKStep
    // --------------

    // Create integrator
    void *arkode_mem = ARKStepCreate(NULL, diffusion, ZERO, u, ctx);
    if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

    // Specify tolerances
    flag = ARKStepSStolerances(arkode_mem, uopts.rtol, uopts.atol);
    if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

    // Attach user data
    flag = ARKStepSetUserData(arkode_mem, (void *) &udata);
    if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

    // Attach linear solver
    flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
    if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

#if defined(USE_SUPERLU_DIST)
    if (uopts.ls == "sludist")
    {
      ARKStepSetJacFn(arkode_mem, diffusion_jac);
      if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;
    }
#endif

    if (uopts.preconditioning)
    {
      // Attach preconditioner
      flag = ARKStepSetPreconditioner(arkode_mem, PSetup, PSolve);
      if (check_flag(&flag, "ARKStepSetPreconditioner", 1)) return 1;

      // Set linear solver setup frequency (update preconditioner)
      flag = ARKStepSetLSetupFrequency(arkode_mem, uopts.msbp);
      if (check_flag(&flag, "ARKStepSetLSetupFrequency", 1)) return 1;
    }

    // Set linear solver tolerance factor
    flag = ARKStepSetEpsLin(arkode_mem, uopts.epslin);
    if (check_flag(&flag, "ARKStepSetEpsLin", 1)) return 1;

    // Select method order
    flag = ARKStepSetOrder(arkode_mem, uopts.order);
    if (check_flag(&flag, "ARKStepSetOrder", 1)) return 1;

    // Set fixed step size or adaptivity method
    if (uopts.hfixed > ZERO)
    {
      flag = ARKStepSetFixedStep(arkode_mem, uopts.hfixed);
      if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;
    }
    else
    {
      flag = ARKStepSetAdaptivityMethod(arkode_mem, uopts.controller, SUNTRUE,
                                        SUNFALSE, NULL);
      if (check_flag(&flag, "ARKStepSetAdaptivityMethod", 1)) return 1;
    }

    // Specify linearly implicit non-time-dependent RHS
    if (uopts.linear)
    {
      flag = ARKStepSetLinear(arkode_mem, 0);
      if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;
    }

    // Set max steps between outputs
    flag = ARKStepSetMaxNumSteps(arkode_mem, uopts.maxsteps);
    if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;

    // Set stopping time
    flag = ARKStepSetStopTime(arkode_mem, udata.tf);
    if (check_flag(&flag, "ARKStepSetStopTime", 1)) return 1;

    // -----------------------
    // Loop over output times
    // -----------------------

    // Optionally run in one step mode for a fixed number of time steps (helpful
    // for debugging)
    int stepmode = ARK_NORMAL;

    if (uopts.onestep)
    {
      uout.nout = uopts.onestep;
      stepmode  = ARK_ONE_STEP;
    }

    realtype t     = ZERO;
    realtype dTout = udata.tf / uout.nout;
    realtype tout  = dTout;

    // Inital output
    flag = uout.open(&udata);
    if (check_flag(&flag, "UserOutput::open", 1)) return 1;

    flag = uout.write(t, u, &udata);
    if (check_flag(&flag, "UserOutput::write", 1)) return 1;

    for (int iout = 0; iout < uout.nout; iout++)
    {
      SUNDIALS_MARK_BEGIN(prof, "Evolve");

      // Evolve in time
      flag = ARKStepEvolve(arkode_mem, tout, u, &t, stepmode);
      if (check_flag(&flag, "ARKStepEvolve", 1)) break;

      SUNDIALS_MARK_END(prof, "Evolve");

      // Output solution and error
      flag = uout.write(t, u, &udata);
      if (check_flag(&flag, "UserOutput::write", 1)) return 1;

      // Update output time
      tout += dTout;
      tout = (tout > udata.tf) ? udata.tf : tout;
    }

    // Close output
    flag = uout.close(&udata);
    if (check_flag(&flag, "UserOutput::close", 1)) return 1;

    // --------------
    // Final outputs
    // --------------

    // Print final integrator stats
    if (outproc)
    {
      cout << "Final integrator statistics:" << endl;
      flag = ARKStepPrintAllStats(arkode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);
      if (check_flag(&flag, "ARKStepPrintAllStats", 1)) return 1;
    }

    // ---------
    // Clean up
    // ---------

    // Free MPI Cartesian communicator
    MPI_Comm_free(&(udata.comm_c));

    ARKStepFree(&arkode_mem);
    SUNLinSolFree(LS);

    // Free the SuperLU_DIST structures (also frees user allocated arrays
    // A_data, A_col_idxs, and A_row_ptrs)
#if defined(USE_SUPERLU_DIST)
    if (uopts.ls == "sludist")
    {
      PStatFree(&A_stat);
      dScalePermstructFree(&A_scaleperm);
      dLUstructFree(&A_lu);
      Destroy_CompRowLoc_Matrix_dist(&A_super);
      superlu_gridexit(&grid);
    }
#endif

    // Free vectors
#if defined(USE_HIP) || defined(USE_CUDA)
    N_VDestroy(N_VGetLocalVector_MPIPlusX(u));
#endif
    N_VDestroy(u);
  }
  // Close scope so objects are destroyed before MPI_Finalize

  SUNContext_Free(&ctx);

  // Finalize MPI
  flag = MPI_Finalize();
  return 0;
}


// -----------------------------------------------------------------------------
// UserOptions Helper functions
// -----------------------------------------------------------------------------


int UserOptions::parse_args(vector<string> &args, bool outproc)
{
  vector<string>::iterator it;

  it = find(args.begin(), args.end(), "--help");
  if (it != args.end())
  {
    if (outproc) help();
    return 0;
  }

  it = find(args.begin(), args.end(), "--rtol");
  if (it != args.end())
  {
    rtol = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--atol");
  if (it != args.end())
  {
    atol = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--fixedstep");
  if (it != args.end())
  {
    hfixed = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--order");
  if (it != args.end())
  {
    order = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--controller");
  if (it != args.end())
  {
    controller = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--maxsteps");
  if (it != args.end())
  {
    maxsteps = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--onestep");
  if (it != args.end())
  {
    onestep = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--nonlinear");
  if (it != args.end())
  {
    linear = false;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--ls");
  if (it != args.end())
  {
    ls = *(it + 1);
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--lsinfo");
  if (it != args.end())
  {
    lsinfo = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--liniters");
  if (it != args.end())
  {
    liniters = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--msbp");
  if (it != args.end())
  {
    msbp = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--epslin");
  if (it != args.end())
  {
    epslin = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  return 0;
}


// Print command line options
void UserOptions::help()
{
  cout << endl;
  cout << "Integrator command line options:" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --atol <atol>           : absoltue tolerance" << endl;
  cout << "  --nonlinear             : disable linearly implicit flag" << endl;
  cout << "  --order <ord>           : method order" << endl;
  cout << "  --fixedstep <step>      : used fixed step size" << endl;
  cout << "  --controller <ctr>      : time step adaptivity controller" << endl;
  cout << "  --ls <cg|gmres|sludist> : linear solver" << endl;
  cout << "  --lsinfo                : output residual history" << endl;
  cout << "  --liniters <iters>      : max number of iterations" << endl;
  cout << "  --epslin <factor>       : linear tolerance factor" << endl;
  cout << "  --noprec                : disable preconditioner" << endl;
  cout << "  --msbp <steps>          : max steps between prec setups" << endl;
}


// Print user options
void UserOptions::print()
{
  cout << endl;
  cout << " Integrator options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << " rtol        = " << rtol        << endl;
  cout << " atol        = " << atol        << endl;
  cout << " hfixed      = " << hfixed      << endl;
  cout << " order       = " << order       << endl;
  cout << " controller  = " << controller  << endl;
  cout << " max steps   = " << maxsteps    << endl;
  cout << " linear RHS  = " << linear      << endl;
  cout << " --------------------------------- " << endl;

  cout << endl;
  if (ls == "sludist")
  {
    cout << " Linear solver options:" << endl;
    cout << " --------------------------------- " << endl;
#if defined(HAVE_HIP)
    cout << " LS       = SuperLU_DIST (HIP enabled)" << endl;
#elif defined(HAVE_CUDA)
    cout << " LS       = SuperLU_DIST (CUDA enabled)" << endl;
#else
    cout << " LS       = SuperLU_DIST" << endl;
#endif
    cout << " LS info  = " << lsinfo   << endl;
    cout << " msbp     = " << msbp     << endl;
    cout << " --------------------------------- " << endl;
  }
  else
  {
    cout << " Linear solver options:" << endl;
    cout << " --------------------------------- " << endl;
    cout << " LS       = " << ls              << endl;
    cout << " precond  = " << preconditioning << endl;
    cout << " LS info  = " << lsinfo          << endl;
    cout << " LS iters = " << liniters        << endl;
    cout << " msbp     = " << msbp            << endl;
    cout << " epslin   = " << epslin          << endl;
    cout << " --------------------------------- " << endl;
  }
}

//---- end of file ----
