/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
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
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"
#include <chrono>
#include <iostream>
#include <iterator>
#include <fstream>
#include "string.h"

struct UserOptions
{
  // Integrator settings
  realtype rtol             = RCONST(1.0e-5);  				// relative tolerance
  realtype atol             = RCONST(1.0e-10);  			// absolute tolerance
  realtype hfixed           = ZERO;             			// fixed step size
  ARKODE_DIRKTableID method = ARKODE_ARK324L2SA_DIRK_4_2_3; // ARKode method
  int controller            = 0;                			// step size adaptivity method
  int maxsteps              = 0;                			// max steps between outputs
  int onestep               = 0;                			// one step mode, number of steps
  bool linear               = true;             			// linearly implicit RHS
  bool diagnostics          = false;            			// output diagnostics
  bool stats	            = false;						// whether to output integrator statistics
  bool savesols             = false;						// whether to save solutions to disc
  bool userefsols           = false;                        // whether to use generated solutions to compute error

  // Nonlinear solver settings
  int maxncf	      = 10;      // max nonlinear solver convergence failures per step
  realtype nlscoef    = 0.1;     // nonlinear solver convergence safety factor
  int maxcor          = 3;       // max nonlinear solver iterations per step
  bool use_fixedpoint = false;   // use fixedpoint iteration over default newton iteration
  int fixedpoint_m    = 0;       // number of anderson acceleration vectors

  // Linear solver and preconditioner settings
  bool     pcg            = true;   // use PCG (true) or GMRES (false)
  bool     prec           = true;   // preconditioner on/off
  bool     lsinfo         = false;  // output residual history
  int      liniters       = 20;     // number of linear iterations
  int      msbp     	  = 0;      // preconditioner setup frequency
  realtype epslin   	  = ZERO;   // linear solver tolerance factor
  int itype				  = 0; 	    // dense output interpolation type
  int idegree 			  = 0;	    // dense output interpolant degree
  bool deduce 			  = false;  // if implicit stage derivatives are deduced without evaluation
  realtype dgmax 	 	  = 0.2; 	// step size ratio limit signalling re-setup of linear solver
  int msbj 				  = 51; 	// preconditioner jacobian setup frequency
  

  // Helper functions
  int parse_args(vector<string> &args, bool outproc);
  void help();
  void print();
};

// Print integration statistics and timings
static int OutputStats(void *arkode_mem, UserData *udata);

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
  MPI_Comm comm = MPI_COMM_WORLD;

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
    // Output problem setup/options
    FILE *diagfp = NULL;

    if (outproc)
    {
      //udata.print();
      //uopts.print();
      //uout.print();

      // Open diagnostics output file
      if (uopts.diagnostics || uopts.lsinfo)
      {
        diagfp = fopen("diagnostics.txt", "w");
        if (check_flag((void *) diagfp, "fopen", 0)) return 1;
      }
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

	// ------------------------
	// Create Nonlinear Solver
	// ------------------------
	// Set fixed-point solver if selected, otherwise set newton solver
	SUNNonlinearSolver NLS = NULL;

	if (uopts.use_fixedpoint) {
		NLS = SUNNonlinSol_FixedPoint(u, uopts.fixedpoint_m, ctx);
      	if (check_flag((void *) NLS, "SUNLinSol_FixedPoint", 0)) return 1;
	} else {
		NLS = SUNNonlinSol_Newton(u, ctx);
      	if (check_flag((void *) NLS, "SUNLinSol_Newton", 0)) return 1;
	}

    // ---------------------
    // Create linear solver
    // ---------------------

    // Create linear solver
    SUNLinearSolver LS = NULL;

    int prectype = (uopts.prec) ? PREC_RIGHT : PREC_NONE;

    if (uopts.pcg)
    {
      LS = SUNLinSol_PCG(u, prectype, uopts.liniters, ctx);
      if (check_flag((void *) LS, "SUNLinSol_PCG", 0)) return 1;

      if (uopts.lsinfo && outproc)
      {
        flag = SUNLinSolSetPrintLevel_PCG(LS, 1);
        if (check_flag(&flag, "SUNLinSolSetPrintLevel_PCG", 1)) return(1);

        flag = SUNLinSolSetInfoFile_PCG(LS, diagfp);
        if (check_flag(&flag, "SUNLinSolSetInfoFile_PCG", 1)) return(1);
      }
    }
    else
    {
      LS = SUNLinSol_SPGMR(u, prectype, uopts.liniters, ctx);
      if (check_flag((void *) LS, "SUNLinSol_SPGMR", 0)) return 1;

      if (uopts.lsinfo && outproc)
      {
        flag = SUNLinSolSetPrintLevel_SPGMR(LS, 1);
        if (check_flag(&flag, "SUNLinSolSetPrintLevel_SPGMR", 1)) return(1);

        flag = SUNLinSolSetInfoFile_SPGMR(LS, diagfp);
        if (check_flag(&flag, "SUNLinSolSetInfoFile_SPGMR", 1)) return(1);
      }
    }

    // Allocate preconditioner workspace
    if (uopts.prec)
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

	// Set max steps to large value 
	flag = ARKStepSetMaxNumSteps(arkode_mem, 100000);
	if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;
	
    // Specify tolerances
    flag = ARKStepSStolerances(arkode_mem, uopts.rtol, uopts.atol);
    if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

    // Attach user data
    flag = ARKStepSetUserData(arkode_mem, (void *) &udata);
    if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

	// Attach nonlinear solver 
	flag = ARKStepSetNonlinearSolver(arkode_mem, NLS);
    if (check_flag(&flag, "ARKStepSetNonlinearSolver", 1)) return 1;
 
    // Attach linear solver
    flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
    if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

    if (uopts.prec)
    {
      // Attach preconditioner
      flag = ARKStepSetPreconditioner(arkode_mem, PSetup, PSolve);
      if (check_flag(&flag, "ARKStepSetPreconditioner", 1)) return 1;

      // Set linear solver setup frequency (update preconditioner)
      flag = ARKStepSetLSetupFrequency(arkode_mem, uopts.msbp);
      if (check_flag(&flag, "ARKStepSetLSetupFrequency", 1)) return 1;

      // Set jacobian evaluation frequency (preconditioner)
      flag = ARKStepSetJacEvalFrequency(arkode_mem, uopts.msbj);
      if (check_flag(&flag, "ARKStepSetJacEvalFrequency", 1)) return 1;
    }
    
	// Set linear solver tolerance factor
    flag = ARKStepSetEpsLin(arkode_mem, uopts.epslin);
    if (check_flag(&flag, "ARKStepSetEpsLin", 1)) return 1;

    // Set dense output interpolation type
    flag = ARKStepSetInterpolantType(arkode_mem, uopts.itype);
    if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) return 1;	

    // Set dense output interpolation degree
    flag = ARKStepSetInterpolantDegree(arkode_mem, uopts.idegree);
    if (check_flag(&flag, "ARKStepSetInterpolantDegree", 1)) return 1;	

    // Set nonlinear solver convergence safety factor
    flag = ARKStepSetNonlinConvCoef(arkode_mem, uopts.nlscoef);
    if (check_flag(&flag, "ARKStepSetNonlinConvCoef", 1)) return 1;	

    // Set max nonlinear solver convergence failures per step
    flag = ARKStepSetMaxConvFails(arkode_mem, uopts.maxncf);
    if (check_flag(&flag, "ARKStepSetMaxConvFails", 1)) return 1;	

    // Set whether to implicit stage derivatives
    flag = ARKStepSetDeduceImplicitRhs(arkode_mem, uopts.deduce);
    if (check_flag(&flag, "ARKStepSetDeduceImplicitRhs", 1)) return 1;	

    // Set step size ratio limit signalling re-setup of linear solver
    flag = ARKStepSetDeltaGammaMax(arkode_mem, uopts.dgmax);
    if (check_flag(&flag, "ARKStepSetDeltaGammaMax", 1)) return 1;	

    // Set step size ratio limit signalling re-setup of linear solver
    flag = ARKStepSetMaxNonlinIters(arkode_mem, uopts.maxcor);
    if (check_flag(&flag, "ARKStepSetMaxNonlinIters", 1)) return 1;	

    // Select method
    flag = ARKStepSetTableNum(arkode_mem, uopts.method, ARKODE_ERK_NONE);
    if (check_flag(&flag, "ARKStepSetTableNum", 1)) return 1;

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
    //flag = ARKStepSetMaxNumSteps(arkode_mem, uopts.maxsteps);
    //if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;

    // Set stopping time
    flag = ARKStepSetStopTime(arkode_mem, udata.tf);
    if (check_flag(&flag, "ARKStepSetStopTime", 1)) return 1;

    // Set diagnostics output file
    if (diagfp)
    {
      flag = ARKStepSetDiagnostics(arkode_mem, diagfp);
      if (check_flag(&flag, "ARKStepSetDiagnostics", 1)) return 1;
    }

	// Load reference solutions if desired
	vector<N_Vector> solutions_vec;
	if (uopts.userefsols) {
	  for (int iin = 0; iin < uout.nout; iin++) {
		char filename[50];
		sprintf(filename, "diffusion_arkode_u_%d_%d.out", iin, myid);
		ifstream solfile(filename);
		vector<realtype> soldata((istream_iterator<realtype>(solfile)),istream_iterator<realtype>());
		realtype* solrawdata = soldata.data();
		N_Vector solution_vec = N_VNew_Parallel(udata.comm_c, NV_LOCLENGTH_P(u), NV_GLOBLENGTH_P(u), ctx);
		memcpy(NV_DATA_P(solution_vec),solrawdata,NV_LOCLENGTH_P(solution_vec)*sizeof(realtype));
		solutions_vec.push_back(solution_vec);
		solfile.close(); 
	  } 
	}

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
	realtype max_error = ZERO;
	realtype curr_error = ZERO;

    // Inital output
    //flag = uout.open(&udata);
    //if (check_flag(&flag, "UserOutput::open", 1)) return 1;

    //flag = uout.write(t, u, &udata);
    //if (check_flag(&flag, "UserOutput::write", 1)) return 1;

	auto start = std::chrono::steady_clock::now();
    for (int iout = 0; iout < uout.nout; iout++)
    {
      SUNDIALS_MARK_BEGIN(prof, "Evolve");

      // Evolve in time
      flag = ARKStepEvolve(arkode_mem, tout, u, &t, stepmode);
      if (check_flag(&flag, "ARKStepEvolve", 1)) break;

      SUNDIALS_MARK_END(prof, "Evolve");
	
	  // Compute solution error
	  if (uopts.userefsols) {
		N_VLinearSum(ONE, u, -ONE, solutions_vec.at(iout), uout.error);
	  } else {
	  	SolutionError(t, u, uout.error, &udata);
	  }
	  max_error = std::max(max_error, N_VMaxNorm(uout.error)); 

      // Output solution and error
      //flag = uout.write(t, u, &udata);
      //if (check_flag(&flag, "UserOutput::write", 1)) return 1;

	  if (uopts.savesols) {
		char filename[50];
		sprintf(filename, "diffusion_arkode_u_%d_%d.out", iout, myid);
		FILE* solfile = fopen(filename, "w");
		N_VPrintFile_Parallel(u, solfile);
		fclose(solfile);
	  }

      // Update output time
      tout += dTout;
      tout = (tout > udata.tf) ? udata.tf : tout;
    }
	auto end = std::chrono::steady_clock::now();
	auto runtime = std::chrono::duration<double>(end - start);
	double runtime_value = (double) runtime.count();

    // Close output
    //flag = uout.close(&udata);
    //if (check_flag(&flag, "UserOutput::close", 1)) return 1;

    // --------------
    // Final outputs
    // --------------

    // Print final integrator stats
    if (outproc && uopts.stats)
    {
      cout << "Final integrator statistics:" << endl;
      flag = OutputStats(arkode_mem, &udata);
      if (check_flag(&flag, "OutputStats", 1)) return 1;
    }
    
	if (outproc) { 
		printf("%.16f,%.16f\n",runtime_value,max_error);
	}

    // ---------
    // Clean up
    // ---------

    // Free MPI Cartesian communicator
    MPI_Comm_free(&(udata.comm_c));

    // Close diagnostics output file
    if (diagfp) fclose(diagfp);

    // Free integrator and linear solver
    ARKStepFree(&arkode_mem);
    SUNLinSolFree(LS);

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
// Output functions
// -----------------------------------------------------------------------------


// Print integrator statistics
static int OutputStats(void *arkode_mem, UserData* udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, nst_a, netf, nfe, nfi, nni, ncfn, nli, nlcf, nsetups, nfi_ls, nJv;
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  if (check_flag(&flag, "ARKStepGetNumSteps", 1)) return -1;
  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_flag(&flag, "ARKStepGetNumStepAttempts", 1)) return -1;
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  if (check_flag(&flag, "ARKStepGetNumErrTestFails", 1)) return -1;
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_flag(&flag, "ARKStepGetNumRhsEvals", 1)) return -1;
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return -1;
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  if (check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1)) return -1;
  flag = ARKStepGetNumLinIters(arkode_mem, &nli);
  if (check_flag(&flag, "ARKStepGetNumLinIters", 1)) return -1;
  flag = ARKStepGetNumLinConvFails(arkode_mem, &nlcf);
  if (check_flag(&flag, "ARKStepGetNumLinConvFails", 1)) return -1;
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  if (check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1)) return -1;
  flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfi_ls);
  if (check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1)) return -1;
  flag = ARKStepGetNumJtimesEvals(arkode_mem, &nJv);
  if (check_flag(&flag, "ARKStepGetNumJtimesEvals", 1)) return -1;

  cout << fixed;
  cout << setprecision(6);

  cout << "  Steps            = " << nst     << endl;
  cout << "  Step attempts    = " << nst_a   << endl;
  cout << "  Error test fails = " << netf    << endl;
  cout << "  RHS evals        = " << nfi     << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << "  NLS fails        = " << ncfn    << endl;
  cout << "  LS iters         = " << nli     << endl;
  cout << "  LS fails         = " << nlcf    << endl;
  cout << "  LS setups        = " << nsetups << endl;
  cout << "  LS RHS evals     = " << nfi_ls  << endl;
  cout << "  Jv products      = " << nJv     << endl;
  cout << endl;

  // Compute average nls iters per step attempt and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nst_a;
  realtype avgli  = (realtype) nli / (realtype) nni;
  cout << "  Avg NLS iters per step attempt = " << avgnli << endl;
  cout << "  Avg LS iters per NLS iter      = " << avgli  << endl;
  cout << endl;

  // Get preconditioner stats
  long int npe, nps;
  flag = ARKStepGetNumPrecEvals(arkode_mem, &npe);
  if (check_flag(&flag, "ARKStepGetNumPrecEvals", 1)) return -1;
  flag = ARKStepGetNumPrecSolves(arkode_mem, &nps);
  if (check_flag(&flag, "ARKStepGetNumPrecSolves", 1)) return -1;

  cout << "  Preconditioner setups = " << npe << endl;
  cout << "  Preconditioner solves = " << nps << endl;
  cout << endl;

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

  it = find(args.begin(), args.end(), "--method");
  if (it != args.end())
  {
	string option = *(it + 1);

	if (option.compare("ARKODE_SDIRK_2_1_2") == 0) method = ARKODE_SDIRK_2_1_2;
	if (option.compare("ARKODE_BILLINGTON_3_3_2") == 0) method = ARKODE_BILLINGTON_3_3_2;
	if (option.compare("ARKODE_TRBDF2_3_3_2") == 0) method = ARKODE_TRBDF2_3_3_2;
	if (option.compare("ARKODE_KVAERNO_4_2_3") == 0) method = ARKODE_KVAERNO_4_2_3;	
	if (option.compare("ARKODE_ARK324L2SA_DIRK_4_2_3") == 0) method = ARKODE_ARK324L2SA_DIRK_4_2_3;	
	if (option.compare("ARKODE_CASH_5_2_4") == 0) method = ARKODE_CASH_5_2_4;	
	if (option.compare("ARKODE_CASH_5_3_4") == 0) method = ARKODE_CASH_5_3_4;	
	if (option.compare("ARKODE_SDIRK_5_3_4") == 0) method = ARKODE_SDIRK_5_3_4;	
	if (option.compare("ARKODE_KVAERNO_5_3_4") == 0) method = ARKODE_KVAERNO_5_3_4;	
	if (option.compare("ARKODE_ARK436L2SA_DIRK_6_3_4") == 0) method = ARKODE_ARK436L2SA_DIRK_6_3_4;	
	if (option.compare("ARKODE_KVAERNO_7_4_5") == 0) method = ARKODE_KVAERNO_7_4_5;	
	if (option.compare("ARKODE_ARK548L2SA_DIRK_8_4_5") == 0) method = ARKODE_ARK548L2SA_DIRK_8_4_5;	
	if (option.compare("ARKODE_ARK437L2SA_DIRK_7_3_4") == 0) method = ARKODE_ARK437L2SA_DIRK_7_3_4;	
	if (option.compare("ARKODE_ARK548L2SAb_DIRK_8_4_5") == 0) method = ARKODE_ARK548L2SAb_DIRK_8_4_5;
	if (option.compare("ARKODE_ESDIRK324L2SA_4_2_3") == 0) method = ARKODE_ESDIRK324L2SA_4_2_3;	
	if (option.compare("ARKODE_ESDIRK325L2SA_5_2_3") == 0) method = ARKODE_ESDIRK325L2SA_5_2_3;	
	if (option.compare("ARKODE_ESDIRK32I5L2SA_5_2_3") == 0) method = ARKODE_ESDIRK32I5L2SA_5_2_3;	
	if (option.compare("ARKODE_ESDIRK436L2SA_6_3_4") == 0) method = ARKODE_ESDIRK436L2SA_6_3_4;	
	if (option.compare("ARKODE_ESDIRK43I6L2SA_6_3_4") == 0) method = ARKODE_ESDIRK43I6L2SA_6_3_4;	
	if (option.compare("ARKODE_QESDIRK436L2SA_6_3_4") == 0) method = ARKODE_QESDIRK436L2SA_6_3_4;	
	if (option.compare("ARKODE_ESDIRK437L2SA_7_3_4") == 0) method = ARKODE_ESDIRK437L2SA_7_3_4;	
	if (option.compare("ARKODE_ESDIRK547L2SA_7_4_5") == 0) method = ARKODE_ESDIRK547L2SA_7_4_5;	
	if (option.compare("ARKODE_ESDIRK547L2SA2_7_4_5") == 0) method = ARKODE_ESDIRK547L2SA2_7_4_5;	
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

  it = find(args.begin(), args.end(), "--diagnostics");
  if (it != args.end())
  {
    diagnostics = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--stats");
  if (it != args.end())
  {
    stats = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--savesols");
  if (it != args.end())
  {
    savesols = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--userefsols");
  if (it != args.end())
  {
    userefsols = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--gmres");
  if (it != args.end())
  {
    pcg = false;
    args.erase(it);
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

  it = find(args.begin(), args.end(), "--itype");
  if (it != args.end())
  {
    itype = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--idegree");
  if (it != args.end())
  {
    idegree = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--nlscoef");
  if (it != args.end())
  {
    nlscoef = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--maxncf");
  if (it != args.end())
  {
    maxncf = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--deduce");
  if (it != args.end())
  {
    deduce = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--dgmax");
  if (it != args.end())
  {
    dgmax = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--msbj");
  if (it != args.end())
  {
    msbj = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--maxcor");
  if (it != args.end())
  {
    maxcor = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--fixedpoint");
  if (it != args.end())
  {
	use_fixedpoint = true;
    fixedpoint_m = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  return 0;
}


// Print command line options
void UserOptions::help()
{
  cout << endl;
  cout << "Integrator command line options:" << endl;
  cout << "  --rtol <rtol>       : relative tolerance" << endl;
  cout << "  --atol <atol>       : absoltue tolerance" << endl;
  cout << "  --nonlinear         : disable linearly implicit flag" << endl;
  cout << "  --method <method>   : method used" << endl;
  cout << "  --fixedstep <step>  : used fixed step size" << endl;
  cout << "  --controller <ctr>  : time step adaptivity controller" << endl;
  cout << "  --diagnostics       : output diagnostics" << endl;
  cout << "  --stats			 : output integrator statistics" << endl;
  cout << "  --savesols			 : save solutions to disc" << endl;
  cout << "  --userefsols		 : use generated solutions for error computation" << endl;
  cout << "  --gmres             : use GMRES linear solver" << endl;
  cout << "  --lsinfo            : output residual history" << endl;
  cout << "  --liniters <iters>  : max number of iterations" << endl;
  cout << "  --epslin <factor>   : linear tolerance factor" << endl;
  cout << "  --noprec            : disable preconditioner" << endl;
  cout << "  --itype <type>      : dense output inteprolant type" << endl;
  cout << "  --idegree <deg>     : dense output interpoland degree" << endl;
  cout << "  --nlscoef <coef>    : nonlinear solver safety factor " << endl;
  cout << "  --maxncf <factor>   : max nonlinear solver congergence failures per step" << endl;
  cout << "  --deduce            : whether to deduce implicit righthand side" << endl;
  cout << "  --dgmax <factor>    : step size ratio limit to re-setup linear sovler" << endl;
  cout << " --------------------------------- " << endl;
  cout << "  --msbp <steps>      : max steps between prec setups" << endl;
  cout << "  --msbj <steps>      : prec jacobian setup frequency" << endl; 
  cout << "  --maxcor <iters>    : max nonlinear iterations per step" << endl;
  cout << "  --fixedpoint <m>    : use fixed point nonlinear solver with m anderson vectors" << endl;
  
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
  cout << " method      = " << method      << endl;
  cout << " controller  = " << controller  << endl;
  cout << " max steps   = " << maxsteps    << endl;
  cout << " linear RHS  = " << linear      << endl;
  cout << " diagnostics = " << diagnostics << endl;
  cout << " stats       = " << stats       << endl;
  cout << " savesols    = " << savesols    << endl;
  cout << " userefsols  = " << userefsols  << endl;
  cout << " --------------------------------- " << endl;

  cout << endl;
  cout << " Linear solver options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << " PCG      = " << pcg      << endl;
  cout << " precond  = " << prec     << endl;
  cout << " LS info  = " << lsinfo   << endl;
  cout << " LS iters = " << liniters << endl;
  cout << " msbp     = " << msbp     << endl;
  cout << " epslin   = " << epslin   << endl;
  cout << " itype    = " << itype    << endl;
  cout << " idegree  = " << idegree  << endl;
  cout << " deduce   = " << deduce   << endl;
  cout << " dgmax    = " << dgmax    << endl;
  cout << " msbj     = " << msbj     << endl; 
  cout << " --------------------------------- " << endl;

  cout << endl;
  cout << " Nonlinear solver options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << " nlscoef        = " << nlscoef  << endl;
  cout << " maxncf         = " << maxncf   << endl;
  cout << " maxcor   = " << maxcor   << endl;
  cout << " fixedpoint_m   = " << fixedpoint_m   << endl;
  cout << " --------------------------------- " << endl;

}

//---- end of file ----
