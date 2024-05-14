/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Rujeko Chinomona @ SMU
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
 *
 *    [u]' = [ G  e ] [(-1+u^2-r)/(2u)] + [      r'(t)/(2u)        ]
 *    [v]    [ e -1 ] [(-2+v^2-s)/(2v)]   [ s'(t)/(2*sqrt(2+s(t))) ]
 *         = [ fs(t,u,v) ]
 *           [ ff(t,u,v) ]
 *
 * where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.
 *
 * This problem has analytical solution given by
 *    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
 *
 * This program allows a number of parameters:
 *   e:    fast/slow coupling strength [default = 0.5]
 *   G:    stiffness at slow time scale [default = -1e2]
 *   w:    time-scale separation factor [default = 100]
 *
 * The stiffness of the slow time scale is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to a
 * multirate method that is implicit at the slow time scale.
 *
 * Additional input options may be used to select between various
 * solver options:
 * - slow step size:  hs [default = 0.01]
 * - fast step size:  hf [default = 0.0001]
 * - set initial adaptive step size as hs/hf above:  set_h0 [default 0]
 * - relative solution tolerance:  rtol [default = 1e-4]
 * - absolute solution tolerance:  atol [default = 1e-11]
 * - use p (0) vs q (1) for slow adaptivity:  slow_pq [default = 0]
 * - use p (0) vs q (1) for fast adaptivity:  fast_pq [default = 0]
 * - "slow" MRI method:  mri_type [default = 4]
 *      0:  explicit Knoth-Wolke MIS [ARKODE_MIS_KW3]
 *      1:  explicit ERK22a MRI-GARK [ARKODE_MRI_GARK_ERK22a]
 *      2:  explicit ERK22b MRI-GARK [ARKODE_MRI_GARK_ERK22b]
 *      3:  explicit ERK33a MRI-GARK [ARKODE_MRI_GARK_ERK33a]
 *      4:  explicit ERK45a MRI-GARK [ARKODE_MRI_GARK_ERK45a]
 *      5:  implicit IRK21a MRI-GARK [ARKODE_MRI_GARK_IRK21a]
 *      6:  implicit ESDIRK34a MRI-GARK [ARKODE_MRI_GARK_ESDIRK34a]
 *      8.  implicit ESDIRK46a MRI-GARK [ARKODE_MRI_GARK_ESDIRK46a]
 *      9.  ImEx 3a MRI-GARK [ARKODE_IMEX_MRI_GARK3a]
 *      10. ImEx 3b MRI-GARK [ARKODE_IMEX_MRI_GARK3b]
 *      11. ImEx 4 MRI-GARK [ARKODE_IMEX_MRI_GARK4]
 * - "fast" ARKStep method: fast_type [default = 1]
 *      0:  none [all physics operators treated at fast scale
 *      1:  default 2nd-order ERK
 *      2:  default 3rd-order ERK
 *      3:  default 4th-order ERK
 *      4:  default 2nd-order DIRK
 *      5:  default 3rd-order DIRK
 *      6:  default 4th-order DIRK
 * - "slow" MRI temporal adaptivity controller: scontrol [default = 6]
 *      0:  no controller [fixed time steps]
 *      1:  MRI-CC controller
 *      2:  MRI-LL controller
 *      3:  MRI-PI controller
 *      4:  MRI-PID controller
 *      5:  I controller (as part of MRI-HTOL)
 *      6:  I controller (alone)
 *      7:  PI controller (as part of MRI-HTOL)
 *      8:  PI controller (alone)
 *      9:  PID controller (as part of MRI-HTOL)
 *      10: PID controller (alone)
 *      11: ExpGus controller (as part of MRI-HTOL)
 *      12: ExpGus controller (alone)
 *      13: ImpGus controller (as part of MRI-HTOL)
 *      14: ImpGus controller (alone)
 *      15: ImExGus controller (as part of MRI-HTOL)
 *      16: ImExGus controller (alone)
 * - "fast" ARKStep temporal adaptivity controller: fcontrol [default = 1]
 *   Note that this will only be used for 5 <= scontrol <= 16.
 *      0:  no controller [fixed time steps]
 *      1:  I controller
 *      2:  PI controller
 *      3:  PID controller
 *      4:  ExpGus controller
 *      5:  ImpGus controller
 *      6:  ImExGus controller
 *
 * Outputs and solution error values are printed at equal intervals
 * of 0.1 and run statistics are printed at the end.
 * ----------------------------------------------------------------*/

// Header files
#include <arkode/arkode_arkstep.h> // prototypes for ARKStep fcts., consts
#include <arkode/arkode_mristep.h> // prototypes for MRIStep fcts., consts
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nvector/nvector_serial.h> // serial N_Vector type, fcts., macros
#include <sunadaptcontroller/sunadaptcontroller_imexgus.h>
#include <sunadaptcontroller/sunadaptcontroller_mricc.h>
#include <sunadaptcontroller/sunadaptcontroller_mrihtol.h>
#include <sunadaptcontroller/sunadaptcontroller_mrill.h>
#include <sunadaptcontroller/sunadaptcontroller_mripi.h>
#include <sunadaptcontroller/sunadaptcontroller_mripid.h>
#include <sunadaptcontroller/sunadaptcontroller_soderlind.h>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_dense.h> // dense linear solver
#include <sunmatrix/sunmatrix_dense.h> // dense matrix type, fcts., macros
#include <test_utilities.hpp>          // common utility functions

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define ESYM "Le"
#define FSYM "Lf"
#else
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)
#define TWO  SUN_RCONST(2.0)

// Problem options
struct Options
{
  // Problem parameters
  sunrealtype e = SUN_RCONST(0.5);
  sunrealtype G = SUN_RCONST(-100.0);
  sunrealtype w = SUN_RCONST(100.0);

  // Step sizes and tolerances
  int set_h0       = 0;
  sunrealtype hs   = SUN_RCONST(1.0e-2);
  sunrealtype hf   = SUN_RCONST(1.0e-4);
  sunrealtype rtol = SUN_RCONST(1.0e-4);
  sunrealtype atol = SUN_RCONST(1.0e-11);

  // Method selection
  int mri_type  = 3;
  int fast_type = 2;
  int scontrol  = 6;
  int fcontrol  = 1;
  int slow_pq   = 0;
  int fast_pq   = 0;
};

// User-supplied functions called by the solver
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);
static int Jf(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Utility functions
static void InputHelp();
static int ReadInputs(std::vector<std::string>& args, Options& opts,
                      SUNContext ctx);
static void PrintSlowAdaptivity(Options opts);
static void PrintFastAdaptivity(Options opts);

// Private function to check function return values
static sunrealtype r(sunrealtype t, Options* opts);
static sunrealtype s(sunrealtype t, Options* opts);
static sunrealtype rdot(sunrealtype t, Options* opts);
static sunrealtype sdot(sunrealtype t, Options* opts);
static sunrealtype utrue(sunrealtype t, Options* opts);
static sunrealtype vtrue(sunrealtype t, Options* opts);
static int Ytrue(sunrealtype t, N_Vector y, Options* opts);

// Main Program
int main(int argc, char* argv[])
{
  // SUNDIALS context object for this simulation
  sundials::Context sunctx;

  // Read input options
  Options opts;
  std::vector<std::string> args(argv + 1, argv + argc);
  int flag = ReadInputs(args, opts, sunctx);
  if (check_flag(flag, "ReadInputs")) return 1;

  // General problem parameters
  sunrealtype T0    = SUN_RCONST(0.0);       // initial time
  sunrealtype Tf    = SUN_RCONST(5.0);       // final time
  sunrealtype dTout = SUN_RCONST(0.5);       // time between outputs
  sunindextype NEQ  = 2;                     // number of dependent vars.
  int Nt            = (int)ceil(Tf / dTout); // number of output times

  // Initial problem output
  //    While traversing these, set various function pointers, table constants, and method orders.
  ARKRhsFn f_fe, f_fi, f_se, f_si;
  ARKLsJacFn J_s, J_f;
  int retval, P, p;
  ARKODE_MRITableID mri_table;
  sunbooleantype slowimplicit, fastimplicit;
  slowimplicit = fastimplicit = SUNFALSE;
  f_fe = f_fi = f_se = f_si = NULL;
  J_s = J_f = NULL;
  std::cout << "\nAdaptive multirate nonlinear Kvaerno-Prothero-Robinson test "
               "problem:\n";
  std::cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  std::cout << "    G = " << opts.G << std::endl;
  std::cout << "    w = " << opts.w << std::endl;
  std::cout << "    e = " << opts.e << std::endl;
  std::cout << "\n  Slow integrator:\n";
  switch (opts.mri_type)
  {
  case (0):
    f_se      = (opts.fast_type == 0) ? fn : fs;
    P         = 3;
    mri_table = ARKODE_MIS_KW3;
    std::cout << "    explicit Knoth-Wolke MIS\n";
    PrintSlowAdaptivity(opts);
    break;
  case (1):
    f_se      = (opts.fast_type == 0) ? fn : fs;
    P         = 2;
    mri_table = ARKODE_MRI_GARK_ERK22a;
    std::cout << "    explicit ERK22a MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (2):
    f_se      = (opts.fast_type == 0) ? fn : fs;
    P         = 2;
    mri_table = ARKODE_MRI_GARK_ERK22b;
    std::cout << "    explicit ERK22b MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (3):
    f_se      = (opts.fast_type == 0) ? fn : fs;
    P         = 3;
    mri_table = ARKODE_MRI_GARK_ERK33a;
    std::cout << "    explicit ERK33a MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (4):
    f_se      = (opts.fast_type == 0) ? fn : fs;
    P         = 4;
    mri_table = ARKODE_MRI_GARK_ERK45a;
    std::cout << "    explicit ERK45a MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (5):
    f_si         = (opts.fast_type == 0) ? fn : fs;
    J_s          = (opts.fast_type == 0) ? Jn : Js;
    P            = 2;
    slowimplicit = SUNTRUE;
    mri_table    = ARKODE_MRI_GARK_IRK21a;
    std::cout << "    implicit IRK21a MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (6):
    f_si         = (opts.fast_type == 0) ? fn : fs;
    J_s          = (opts.fast_type == 0) ? Jn : Js;
    P            = 3;
    slowimplicit = SUNTRUE;
    mri_table    = ARKODE_MRI_GARK_ESDIRK34a;
    std::cout << "    implicit ESDIRK34a MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (8):
    f_si         = (opts.fast_type == 0) ? fn : fs;
    J_s          = (opts.fast_type == 0) ? Jn : Js;
    P            = 4;
    slowimplicit = SUNTRUE;
    mri_table    = ARKODE_MRI_GARK_ESDIRK46a;
    std::cout << "    implicit ESDIRK46a MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (9):
    f_se         = (opts.fast_type == 0) ? fs : fse;
    f_si         = (opts.fast_type == 0) ? ff : fsi;
    J_s          = (opts.fast_type == 0) ? Jf : Jsi;
    P            = 3;
    slowimplicit = SUNTRUE;
    mri_table    = ARKODE_IMEX_MRI_GARK3a;
    std::cout << "    ImEx 3a MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (10):
    f_se         = (opts.fast_type == 0) ? fs : fse;
    f_si         = (opts.fast_type == 0) ? ff : fsi;
    J_s          = (opts.fast_type == 0) ? Jf : Jsi;
    P            = 3;
    slowimplicit = SUNTRUE;
    mri_table    = ARKODE_IMEX_MRI_GARK3b;
    std::cout << "    ImEx 3b MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  case (11):
    f_se         = (opts.fast_type == 0) ? fs : fse;
    f_si         = (opts.fast_type == 0) ? ff : fsi;
    J_s          = (opts.fast_type == 0) ? Jf : Jsi;
    P            = 4;
    slowimplicit = SUNTRUE;
    mri_table    = ARKODE_IMEX_MRI_GARK4;
    std::cout << "    ImEx 4 MRI-GARK\n";
    PrintSlowAdaptivity(opts);
    break;
  }
  std::cout << "\n  Fast integrator:\n";
  switch (opts.fast_type)
  {
  case (0):
    f_fe = f0;
    p    = 2;
    std::cout
      << "    zero-valued operator, integrated with default 2nd-order ERK\n";
    PrintFastAdaptivity(opts);
    break;
  case (1):
    f_fe = ff;
    p    = 2;
    std::cout << "    Default 2nd-order ERK\n";
    PrintFastAdaptivity(opts);
    break;
  case (2):
    f_fe = ff;
    p    = 3;
    std::cout << "    Default 3rd-order ERK\n";
    PrintFastAdaptivity(opts);
    break;
  case (3):
    f_fe = ff;
    p    = 4;
    std::cout << "    Default 4th-order ERK\n";
    PrintFastAdaptivity(opts);
    break;
  case (4):
    f_fi         = ff;
    J_f          = Jf;
    p            = 2;
    fastimplicit = SUNTRUE;
    std::cout << "    Default 2nd-order DIRK\n";
    PrintFastAdaptivity(opts);
    break;
  case (5):
    f_fi         = ff;
    J_f          = Jf;
    p            = 3;
    fastimplicit = SUNTRUE;
    std::cout << "    Default 3rd-order DIRK\n";
    PrintFastAdaptivity(opts);
    break;
  case (6):
    f_fi         = ff;
    J_f          = Jf;
    p            = 4;
    fastimplicit = SUNTRUE;
    std::cout << "    Default 4th-order DIRK\n";
    PrintFastAdaptivity(opts);
    break;
  }

  // Create and initialize serial vector for the solution
  N_Vector y = NULL;
  y          = N_VNew_Serial(NEQ, sunctx);
  if (check_ptr((void*)y, "N_VNew_Serial")) return 1;
  retval = Ytrue(T0, y, &opts);
  if (check_flag(retval, "Ytrue")) return 1;

  // Create and configure fast controller object
  SUNAdaptController fcontrol = NULL;
  switch (opts.fcontrol)
  {
  case (1):
    fcontrol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_I")) return 1;
    break;
  case (2):
    fcontrol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_PI")) return 1;
    break;
  case (3):
    fcontrol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_PID")) return 1;
    break;
  case (4):
    fcontrol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ExpGus")) return 1;
    break;
  case (5):
    fcontrol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ImpGus")) return 1;
    break;
  case (6):
    fcontrol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)fcontrol, "SUNAdaptController_ImExGus")) return 1;
    break;
  }

  // Create ARKStep (fast) integrator, storing desired adaptivity order in p
  void* inner_arkode_mem = NULL; // ARKode memory structure
  inner_arkode_mem       = ARKStepCreate(f_fe, f_fi, T0, y, sunctx);
  if (check_ptr((void*)inner_arkode_mem, "ARKStepCreate")) return 1;
  retval = ARKodeSetOrder(inner_arkode_mem, p);
  if (check_flag(retval, "ARKodeSetOrder")) return 1;
  SUNMatrix Af        = NULL; // matrix for fast solver
  SUNLinearSolver LSf = NULL; // fast linear solver object
  if (fastimplicit)
  {
    Af = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_ptr((void*)Af, "SUNDenseMatrix")) return 1;
    LSf = SUNLinSol_Dense(y, Af, sunctx);
    if (check_ptr((void*)LSf, "SUNLinSol_Dense")) return 1;
    retval = ARKodeSetLinearSolver(inner_arkode_mem, LSf, Af);
    if (check_flag(retval, "ARKodeSetLinearSolver")) return 1;
    retval = ARKodeSetJacFn(inner_arkode_mem, J_f);
    if (check_flag(retval, "ARKodeSetJacFn")) return 1;
  }
  retval = ARKodeSStolerances(inner_arkode_mem, opts.rtol, opts.atol);
  if (check_flag(retval, "ARKodeSStolerances")) return 1;
  if (opts.fcontrol != 0)
  {
    retval = ARKodeSetAdaptController(inner_arkode_mem, fcontrol);
    if (check_flag(retval, "ARKodeSetAdaptController")) return 1;
    if (opts.set_h0 != 0)
    {
      retval = ARKodeSetInitStep(inner_arkode_mem, opts.hf);
      if (check_flag(retval, "ARKodeSetInitStep")) return 1;
    }
    if (opts.fast_pq == 1)
    {
      retval = ARKodeSetAdaptivityAdjustment(inner_arkode_mem, 0);
      if (check_flag(retval, "ARKodeSetAdaptivityAdjustment")) return 1;
    }
  }
  else
  {
    retval = ARKodeSetFixedStep(inner_arkode_mem, opts.hf);
    if (check_flag(retval, "ARKodeSetFixedStep")) return 1;
  }
  retval = ARKodeSetMaxNumSteps(inner_arkode_mem, 10000);
  if (check_flag(retval, "ARKodeSetMaxNumSteps")) return 1;
  retval = ARKodeSetUserData(inner_arkode_mem, (void*)&opts);
  if (check_flag(retval, "ARKodeSetUserData")) return 1;

  // Create inner stepper
  MRIStepInnerStepper inner_stepper = NULL; // inner stepper
  retval = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &inner_stepper);
  if (check_flag(retval, "ARKStepCreateMRIStepInnerStepper")) return 1;

  // Create slow controller object, and select orders of accuracy as relevant
  SUNAdaptController scontrol       = NULL;
  SUNAdaptController scontrol_inner = NULL;
  switch (opts.scontrol)
  {
  case (1):
    scontrol = SUNAdaptController_MRICC(sunctx, p);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRICC")) return 1;
    break;
  case (2):
    scontrol = SUNAdaptController_MRILL(sunctx, p);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRILL")) return 1;
    break;
  case (3):
    scontrol = SUNAdaptController_MRIPI(sunctx, p);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIPI")) return 1;
    break;
  case (4):
    scontrol = SUNAdaptController_MRIPID(sunctx, p);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIPID")) return 1;
    break;
  case (5):
    scontrol_inner = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol_inner, "SUNAdaptController_I (slow)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(sunctx, scontrol_inner, fcontrol);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    break;
  case (6):
    scontrol = SUNAdaptController_I(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptControllerI (slow)")) return 1;
    break;
  case (7):
    scontrol_inner = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol_inner, "SUNAdaptController_PI (slow)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(sunctx, scontrol_inner, fcontrol);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    break;
  case (8):
    scontrol = SUNAdaptController_PI(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_PI (slow)")) return 1;
    break;
  case (9):
    scontrol_inner = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol_inner, "SUNAdaptController_PID (slow)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(sunctx, scontrol_inner, fcontrol);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    break;
  case (10):
    scontrol = SUNAdaptController_PID(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_PID (slow)")) return 1;
    break;
  case (11):
    scontrol_inner = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol_inner, "SUNAdaptController_ExpGus (slow)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(sunctx, scontrol_inner, fcontrol);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    break;
  case (12):
    scontrol = SUNAdaptController_ExpGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ExpGus (slow)"))
      return 1;
    break;
  case (13):
    scontrol_inner = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol_inner, "SUNAdaptController_ImpGus (slow)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(sunctx, scontrol_inner, fcontrol);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    break;
  case (14):
    scontrol = SUNAdaptController_ImpGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ImpGus (slow)"))
      return 1;
    break;
  case (15):
    scontrol_inner = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol_inner, "SUNAdaptController_ImExGus (slow)"))
      return 1;
    scontrol = SUNAdaptController_MRIHTol(sunctx, scontrol_inner, fcontrol);
    if (check_ptr((void*)scontrol, "SUNAdaptController_MRIHTol")) return 1;
    break;
  case (16):
    scontrol = SUNAdaptController_ImExGus(sunctx);
    if (check_ptr((void*)scontrol, "SUNAdaptController_ImExGus (slow)"))
      return 1;
    break;
  }

  // Create MRI (slow) integrator, storing desired adaptivity order in P
  void* arkode_mem = NULL; // ARKode memory structure
  arkode_mem       = MRIStepCreate(f_se, f_si, T0, y, inner_stepper, sunctx);
  if (check_ptr((void*)arkode_mem, "MRIStepCreate")) return 1;
  MRIStepCoupling C = NULL; // slow coupling coefficients
  C                 = MRIStepCoupling_LoadTable(mri_table);
  if (check_ptr((void*)C, "MRIStepCoupling_LoadTable")) return 1;
  retval = MRIStepSetCoupling(arkode_mem, C);
  if (check_flag(retval, "MRIStepSetCoupling")) return 1;
  SUNMatrix As        = NULL; // matrix for slow solver
  SUNLinearSolver LSs = NULL; // slow linear solver object
  if (slowimplicit)
  {
    As = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_ptr((void*)As, "SUNDenseMatrix")) return 1;
    LSs = SUNLinSol_Dense(y, As, sunctx);
    if (check_ptr((void*)LSs, "SUNLinSol_Dense")) return 1;
    retval = ARKodeSetLinearSolver(arkode_mem, LSs, As);
    if (check_flag(retval, "ARKodeSetLinearSolver")) return 1;
    retval = ARKodeSetJacFn(arkode_mem, J_s);
    if (check_flag(retval, "ARKodeSetJacFn")) return 1;
  }
  retval = ARKodeSStolerances(arkode_mem, opts.rtol, opts.atol);
  if (check_flag(retval, "ARKodeSStolerances")) return 1;
  retval = ARKodeSetMaxNumSteps(arkode_mem, 10000);
  if (check_flag(retval, "ARKodeSetMaxNumSteps")) return 1;
  retval = ARKodeSetUserData(arkode_mem, (void*)&opts);
  if (check_flag(retval, "ARKodeSetUserData")) return 1;
  if (opts.scontrol != 0)
  {
    retval = MRIStepSetAdaptController(arkode_mem, scontrol);
    if (check_flag(retval, "MRIStepSetAdaptController")) return 1;
    if (opts.set_h0 != 0)
    {
      retval = ARKodeSetInitStep(arkode_mem, opts.hs);
      if (check_flag(retval, "ARKodeSetInitStep")) return 1;
    }
    if (opts.slow_pq == 1)
    {
      retval = ARKodeSetAdaptivityAdjustment(arkode_mem, 0);
      if (check_flag(retval, "ARKodeSetAdaptivityAdjustment")) return 1;
    }
  }
  else
  {
    retval = ARKodeSetFixedStep(arkode_mem, opts.hs);
    if (check_flag(retval, "ARKodeSetFixedStep")) return 1;
  }

  //
  // Integrate ODE
  //

  // Open output stream for results, output comment line
  FILE* UFID = NULL;
  UFID       = fopen("ark_kpr_mri_solution.txt", "w");
  fprintf(UFID, "# t u v uerr verr\n");

  // output initial condition to disk
  fprintf(UFID,
          " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n",
          T0, NV_Ith_S(y, 0), NV_Ith_S(y, 1),
          SUNRabs(NV_Ith_S(y, 0) - utrue(T0, &opts)),
          SUNRabs(NV_Ith_S(y, 1) - vtrue(T0, &opts)));

  // Main time-stepping loop: calls ARKodeEvolve to perform the
  // integration, then prints results. Stops when the final time
  // has been reached
  sunrealtype t, tout;
  sunrealtype uerr, verr, uerrtot, verrtot, errtot;
  t       = T0;
  tout    = T0 + dTout;
  uerr    = ZERO;
  verr    = ZERO;
  uerrtot = ZERO;
  verrtot = ZERO;
  errtot  = ZERO;
  printf("        t           u           v       uerr      verr\n");
  printf("   ------------------------------------------------------\n");
  printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %.2" ESYM "  %.2" ESYM
         "\n",
         t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), uerr, verr);

  for (int iout = 0; iout < Nt; iout++)
  {
    // call integrator
    retval = ARKodeEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(retval, "ARKodeEvolve")) break;

    // access/print solution and error
    uerr = SUNRabs(NV_Ith_S(y, 0) - utrue(t, &opts));
    verr = SUNRabs(NV_Ith_S(y, 1) - vtrue(t, &opts));
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %.2" ESYM
           "  %.2" ESYM "\n",
           t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), uerr, verr);
    fprintf(UFID,
            " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM "\n",
            t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), uerr, verr);
    uerrtot += uerr * uerr;
    verrtot += verr * verr;
    errtot += uerr * uerr + verr * verr;

    // successful solve: update time
    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;
  }
  uerrtot = SUNRsqrt(uerrtot / Nt);
  verrtot = SUNRsqrt(verrtot / Nt);
  errtot  = SUNRsqrt(errtot / Nt / 2);
  printf("   ------------------------------------------------------\n");
  fclose(UFID);

  //
  // Finalize
  //

  // Get some slow integrator statistics
  long int nsts, natts, netfs, nfse, nfsi;
  retval = ARKodeGetNumSteps(arkode_mem, &nsts);
  check_flag(retval, "ARKodeGetNumSteps");
  retval = ARKodeGetNumStepAttempts(arkode_mem, &natts);
  check_flag(retval, "ARKodeGetNumStepAttempts");
  retval = ARKodeGetNumErrTestFails(arkode_mem, &netfs);
  check_flag(retval, "ARKodeGetNumErrTestFails");
  retval = MRIStepGetNumRhsEvals(arkode_mem, &nfse, &nfsi);
  check_flag(retval, "MRIStepGetNumRhsEvals");

  // Get some fast integrator statistics
  long int nstf, nattf, netff, nffe, nffi;
  retval = ARKodeGetNumSteps(inner_arkode_mem, &nstf);
  check_flag(retval, "ARKodeGetNumSteps");
  retval = ARKodeGetNumStepAttempts(inner_arkode_mem, &nattf);
  check_flag(retval, "ARKodeGetNumStepAttempts");
  retval = ARKodeGetNumErrTestFails(inner_arkode_mem, &netff);
  check_flag(retval, "ARKodeGetNumErrTestFails");
  retval = ARKStepGetNumRhsEvals(inner_arkode_mem, &nffe, &nffi);
  check_flag(retval, "ARKStepGetNumRhsEvals");

  // Print some final statistics
  std::cout << "\nFinal Solver Statistics:\n";
  std::cout << "   Slow steps = " << nsts << "  (attempts = " << natts
            << ",  fails = " << netfs << ")\n";
  std::cout << "   Fast steps = " << nstf << "  (attempts = " << nattf
            << ",  fails = " << netff << ")\n";
  std::cout << "   u error = " << uerrtot << ", v error = " << verrtot
            << ", total error = " << errtot << std::endl;
  std::cout << "   Total RHS evals:  Fse = " << nfse << ", Fsi = " << nfsi
            << ", Ffe = " << nffe << ", Ffi = " << nffi << std::endl;

  // Get/print slow integrator decoupled implicit solver statistics
  if (slowimplicit)
  {
    long int nnis, nncs, njes;
    retval = ARKodeGetNonlinSolvStats(arkode_mem, &nnis, &nncs);
    check_flag(retval, "ARKodeGetNonlinSolvStats");
    retval = ARKodeGetNumJacEvals(arkode_mem, &njes);
    check_flag(retval, "ARKodeGetNumJacEvals");
    std::cout << "   Slow Newton iters = " << nnis << std::endl;
    std::cout << "   Slow Newton conv fails = " << nncs << std::endl;
    std::cout << "   Slow Jacobian evals = " << njes << std::endl;
  }

  // Get/print fast integrator implicit solver statistics
  if (fastimplicit)
  {
    long int nnif, nncf, njef;
    retval = ARKodeGetNonlinSolvStats(inner_arkode_mem, &nnif, &nncf);
    check_flag(retval, "ARKodeGetNonlinSolvStats");
    retval = ARKodeGetNumJacEvals(inner_arkode_mem, &njef);
    check_flag(retval, "ARKodeGetNumJacEvals");
    std::cout << "   Fast Newton iters = " << nnif << std::endl;
    std::cout << "   Fast Newton conv fails = " << nncf << std::endl;
    std::cout << "   Fast Jacobian evals = " << njef << std::endl;
  }

  // Clean up and return
  N_VDestroy(y);           // Free y vector
  MRIStepCoupling_Free(C); // free coupling coefficients
  if (fastimplicit)
  {
    SUNMatDestroy(Af);  // free fast matrix
    SUNLinSolFree(LSf); // free fast linear solver
  }
  if (slowimplicit)
  {
    SUNMatDestroy(As);  // free slow matrix
    SUNLinSolFree(LSs); // free slow linear solver
  }
  // if (opts.scontrol != 0) {
  //   SUNAdaptControllerDestroy(scontrol);             // free slow controller(s)
  //   if (opts.scontrol > 4) {
  //     SUNAdaptControllerDestroy(scontrol_inner);
  //   }
  // }
  // if (opts.fcontrol != 0) {
  //   SUNAdaptControllerDestroy(fcontrol);             // free slow controller
  // }
  ARKodeFree(&inner_arkode_mem);            // Free fast integrator memory
  MRIStepInnerStepper_Free(&inner_stepper); // Free inner stepper structure
  ARKodeFree(&arkode_mem);                  // Free slow integrator memory

  return 0;
}

// ------------------------------
// Functions called by the solver
// -----------------------------

// ff routine to compute the fast portion of the ODE RHS.
static int ff(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
  //   [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))]
  tmp1              = (-ONE + u * u - r(t, opts)) / (TWO * u);
  tmp2              = (-TWO + v * v - s(t, opts)) / (TWO * v);
  NV_Ith_S(ydot, 0) = ZERO;
  NV_Ith_S(ydot, 1) = opts->e * tmp1 - tmp2 +
                      sdot(t, opts) / (TWO * vtrue(t, opts));

  // Return with success
  return 0;
}

// fs routine to compute the slow portion of the ODE RHS.
static int fs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
  //   [0 0] [(-2+v^2-s(t))/(2*v)]    [      0      ]
  tmp1 = (-ONE + u * u - r(t, opts)) / (TWO * u);
  tmp2 = (-TWO + v * v - s(t, opts)) / (TWO * v);
  NV_Ith_S(ydot, 0) = opts->G * tmp1 + opts->e * tmp2 + rdot(t, opts) / (TWO * u);
  NV_Ith_S(ydot, 1) = ZERO;

  // Return with success
  return 0;
}

// fse routine to compute the slow portion of the ODE RHS.
static int fse(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);

  // fill in the slow explicit RHS function:
  //   [rdot(t)/(2*u)]
  //   [      0      ]
  NV_Ith_S(ydot, 0) = rdot(t, opts) / (TWO * u);
  NV_Ith_S(ydot, 1) = ZERO;

  // Return with success
  return 0;
}

// fsi routine to compute the slow portion of the ODE RHS.(currently same as fse)
static int fsi(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  // fill in the slow implicit RHS function:
  //   [G e]*[(-1+u^2-r(t))/(2*u))]
  //   [0 0] [(-2+v^2-s(t))/(2*v)]
  tmp1              = (-ONE + u * u - r(t, opts)) / (TWO * u);
  tmp2              = (-TWO + v * v - s(t, opts)) / (TWO * v);
  NV_Ith_S(ydot, 0) = opts->G * tmp1 + opts->e * tmp2;
  NV_Ith_S(ydot, 1) = ZERO;

  // Return with success
  return 0;
}

static int fn(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);
  sunrealtype tmp1, tmp2;

  // fill in the RHS function:
  //   [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
  //   [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))]
  tmp1 = (-ONE + u * u - r(t, opts)) / (TWO * u);
  tmp2 = (-TWO + v * v - s(t, opts)) / (TWO * v);
  NV_Ith_S(ydot, 0) = opts->G * tmp1 + opts->e * tmp2 + rdot(t, opts) / (TWO * u);
  NV_Ith_S(ydot, 1) = opts->e * tmp1 - tmp2 +
                      sdot(t, opts) / (TWO * vtrue(t, opts));

  // Return with success
  return 0;
}

static int f0(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VConst(ZERO, ydot);
  return (0);
}

static int Jf(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  // fill in the Jacobian:
  //   [         0                            0         ]
  //   [e/2+e*(1+r(t))/(2*u^2)   -1/2 - (2+s(t))/(2*v^2)]
  SM_ELEMENT_D(J, 0, 0) = ZERO;
  SM_ELEMENT_D(J, 0, 1) = ZERO;
  SM_ELEMENT_D(J, 1, 0) = opts->e / TWO +
                          opts->e * (ONE + r(t, opts)) / (TWO * u * u);
  SM_ELEMENT_D(J, 1, 1) = -ONE / TWO - (TWO + s(t, opts)) / (TWO * v * v);

  // Return with success
  return 0;
}

static int Js(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  // fill in the Jacobian:
  //   [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)   e/2+e*(2+s(t))/(2*v^2)]
  //   [                 0                             0           ]
  SM_ELEMENT_D(J, 0, 0) =
    opts->G / TWO + (opts->G * (ONE + r(t, opts)) - rdot(t, opts)) / (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = opts->e / TWO +
                          opts->e * (TWO + s(t, opts)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;

  // Return with success
  return 0;
}

static int Jsi(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  // fill in the Jacobian:
  //   [G/2 + (G*(1+r(t)))/(2*u^2)   e/2 + e*(2+s(t))/(2*v^2)]
  //   [                 0                       0           ]
  SM_ELEMENT_D(J, 0, 0) = opts->G / TWO +
                          (opts->G * (ONE + r(t, opts))) / (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = opts->e / TWO +
                          opts->e * (TWO + s(t, opts)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = ZERO;
  SM_ELEMENT_D(J, 1, 1) = ZERO;

  // Return with success
  return 0;
}

static int Jn(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  Options* opts       = (Options*)user_data;
  const sunrealtype u = NV_Ith_S(y, 0);
  const sunrealtype v = NV_Ith_S(y, 1);

  // fill in the Jacobian:
  //   [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)     e/2 + e*(2+s(t))/(2*v^2)]
  //   [e/2+e*(1+r(t))/(2*u^2)                -1/2 - (2+s(t))/(2*v^2)  ]
  SM_ELEMENT_D(J, 0, 0) =
    opts->G / TWO + (opts->G * (ONE + r(t, opts)) - rdot(t, opts)) / (2 * u * u);
  SM_ELEMENT_D(J, 0, 1) = opts->e / TWO +
                          opts->e * (TWO + s(t, opts)) / (TWO * v * v);
  SM_ELEMENT_D(J, 1, 0) = opts->e / TWO +
                          opts->e * (ONE + r(t, opts)) / (TWO * u * u);
  SM_ELEMENT_D(J, 1, 1) = -ONE / TWO - (TWO + s(t, opts)) / (TWO * v * v);

  // Return with success
  return 0;
}

// ------------------------------
// Private helper functions
// -----------------------------

static sunrealtype r(sunrealtype t, Options* opts)
{
  return (SUN_RCONST(0.5) * cos(t));
}

static sunrealtype s(sunrealtype t, Options* opts)
{
  return (cos(opts->w * t));
}

static sunrealtype rdot(sunrealtype t, Options* opts)
{
  return (-SUN_RCONST(0.5) * sin(t));
}

static sunrealtype sdot(sunrealtype t, Options* opts)
{
  return (-opts->w * sin(opts->w * t));
}

static sunrealtype utrue(sunrealtype t, Options* opts)
{
  return (SUNRsqrt(ONE + r(t, opts)));
}

static sunrealtype vtrue(sunrealtype t, Options* opts)
{
  return (SUNRsqrt(TWO + s(t, opts)));
}

static int Ytrue(sunrealtype t, N_Vector y, Options* opts)
{
  NV_Ith_S(y, 0) = utrue(t, opts);
  NV_Ith_S(y, 1) = vtrue(t, opts);
  return (0);
}

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Print command line options
void InputHelp()
{
  std::cout << std::endl;
  std::cout << "Command line options:" << std::endl;
  std::cout << "  --help         : print options and exit\n";
  std::cout << "  --e            : fast/slow coupling strength\n";
  std::cout << "  --G            : stiffness at slow time scale\n";
  std::cout << "  --w            : time-scale separation factor\n";
  std::cout << "  --hs           : slow (fixed/initial) step size\n";
  std::cout << "  --hf           : fast (fixed/initial) step size\n";
  std::cout
    << "  --set_h0       : use hs/hf above to set the initial step size\n";
  std::cout << "  --rtol         : relative solution tolerance\n";
  std::cout << "  --atol         : absolute solution tolerance\n";
  std::cout << "  --mri_type     : MRI method, int in [0,11] (see source for "
               "explanation)\n";
  std::cout << "  --fast_type    : fast RK method, int in [0,5] (see source)\n";
  std::cout << "  --scontrol     : slow time step controller, int in [0,16] "
               "(see source)\n";
  std::cout << "  --fcontrol     : fast time step controller, int in [0,6] "
               "(see source)\n";
  std::cout << "  --slow_pq      : use p (0) vs q (1) for slow adaptivity\n";
  std::cout << "  --fast_pq      : use p (0) vs q (1) for fast adaptivity\n";
}

// Read input options
int ReadInputs(std::vector<std::string>& args, Options& opts, SUNContext ctx)
{
  if (find(args.begin(), args.end(), "--help") != args.end())
  {
    InputHelp();
    return 1;
  }

  // Problem options
  find_arg(args, "--e", opts.e);
  find_arg(args, "--G", opts.G);
  find_arg(args, "--w", opts.w);
  find_arg(args, "--hs", opts.hs);
  find_arg(args, "--hf", opts.hf);
  find_arg(args, "--set_h0", opts.set_h0);
  find_arg(args, "--rtol", opts.rtol);
  find_arg(args, "--atol", opts.atol);
  find_arg(args, "--mri_type", opts.mri_type);
  find_arg(args, "--fast_type", opts.fast_type);
  find_arg(args, "--scontrol", opts.scontrol);
  find_arg(args, "--fcontrol", opts.fcontrol);
  find_arg(args, "--slow_pq", opts.slow_pq);
  find_arg(args, "--fast_pq", opts.fast_pq);

  // Check inputs for validity
  //   0 < rtol < 1
  if ((opts.rtol < ZERO) || (opts.rtol > ONE))
  {
    std::cerr << "ERROR: rtol must be in (0,1), (" << opts.rtol << " input)\n";
    return -1;
  }
  //   0 < atol < 1
  if ((opts.atol < ZERO) || (opts.atol > ONE))
  {
    std::cerr << "ERROR: atol must be in (0,1), (" << opts.atol << " input)\n";
    return -1;
  }
  //   slow_pq in {0,1}
  if ((opts.slow_pq < 0) || (opts.slow_pq > 1))
  {
    std::cerr << "ERROR: slow_pq must be in {0,1}, (" << opts.slow_pq
              << " input)\n";
    return -1;
  }
  //   fast_pq in {0,1}
  if ((opts.fast_pq < 0) || (opts.fast_pq > 1))
  {
    std::cerr << "ERROR: fast_pq must be in {0,1}, (" << opts.fast_pq
              << " input)\n";
    return -1;
  }
  //   mri_type in [0,11]
  if ((opts.mri_type < 0) || (opts.mri_type > 11))
  {
    std::cerr << "ERROR: mri_type must be in [0,11], (" << opts.mri_type
              << " input)\n";
    return -1;
  }
  //   scontrol = 0 if mri_type in {0,9,10,11}
  if ((opts.scontrol != 0) && ((opts.mri_type == 0) || (opts.mri_type == 9) ||
                               (opts.mri_type == 10) || (opts.mri_type == 11)))
  {
    std::cerr << "ERROR: scontrol must be 0 for mri_type " << opts.mri_type
              << ", (" << opts.scontrol << " input)\n";
    return -1;
  }
  //   fast_type in [0,6]
  if ((opts.fast_type < 0) || (opts.fast_type > 6))
  {
    std::cerr << "ERROR: fast_type must be in [0,6], (" << opts.fast_type
              << " input)\n";
    return -1;
  }
  //   scontrol in [0,16]
  if ((opts.scontrol < 0) || (opts.scontrol > 16))
  {
    std::cerr << "ERROR: scontrol must be in [0,16], (" << opts.scontrol
              << " input)\n";
    return -1;
  }
  //   fcontrol in [0,6]
  if ((opts.fcontrol < 0) || (opts.fcontrol > 6))
  {
    std::cerr << "ERROR: fcontrol must be in [0,6], (" << opts.fcontrol
              << " input)\n";
    return -1;
  }
  //   hs > 0 if scontrol == 0
  if ((opts.hs <= 0) && (opts.scontrol == 0))
  {
    std::cerr << "ERROR: positive hs required with scontrol = 0, (" << opts.hs
              << " input)\n";
    return -1;
  }
  //   hf > 0 if fcontrol == 0
  if ((opts.hf <= 0) && (opts.fcontrol == 0))
  {
    std::cerr << "ERROR: positive hf required with fcontrol = 0, (" << opts.hf
              << " input)\n";
    return -1;
  }
  //   G < 0.0
  if (opts.G >= ZERO)
  {
    std::cerr << "ERROR: G must be a negative real number, (" << opts.G
              << " input)\n";
    return -1;
  }
  //   hs <= 1/|G| if scontrol == 0 and mri_type in [0,4]
  if ((opts.hs > ONE / SUNRabs(opts.G)) && (opts.scontrol == 0) &&
      (opts.mri_type >= 0) && (opts.mri_type < 5))
  {
    std::cerr
      << "ERROR: hs must be in (0, 1/|G|) for fixed-step explicit MRI, ("
      << opts.hs << " input)\n";
    return -1;
  }
  //   w >= 1.0
  if (opts.w < ONE)
  {
    std::cerr << "ERROR: w must be >= 1.0, (" << opts.w << " input)\n";
    return -1;
  }

  return 0;
}

static void PrintSlowAdaptivity(Options opts)
{
  switch (opts.scontrol)
  {
  case (0):
    std::cout << "    fixed steps, hs = " << opts.hs << std::endl;
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (1):
    std::cout << "    MRI-CC controller based on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (2):
    std::cout << "    MRI-LL controller based on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (3):
    std::cout << "    MRI-PI controller based on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (4):
    std::cout << "    MRI-PID controller based on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (5):
    std::cout
      << "    MRI-HTOL controller (using I for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (6):
    std::cout << "    Decoupled I controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (7):
    std::cout
      << "    MRI-HTOL controller (using PI for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (8):
    std::cout << "    Decoupled PI controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (9):
    std::cout
      << "    MRI-HTOL controller (using PID for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (10):
    std::cout << "    Decoupled PID controller for slow time scale, based on "
                 "order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (11):
    std::cout
      << "    MRI-HTOL controller (using ExpGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (12):
    std::cout << "    Decoupled ExpGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (13):
    std::cout
      << "    MRI-HTOL controller (using ImpGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (14):
    std::cout << "    Decoupled ImpGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (15):
    std::cout
      << "    MRI-HTOL controller (using ImExGus for H) based on order of MRI "
      << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (16):
    std::cout << "    Decoupled ImExGus controller for slow time scale, based "
                 "on order of MRI "
              << ((opts.slow_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  }
}

static void PrintFastAdaptivity(Options opts)
{
  switch (opts.fcontrol)
  {
  case (0):
    std::cout << "    fixed steps, hf = " << opts.hf << std::endl;
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (1):
    std::cout << "    I controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (2):
    std::cout << "    PI controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (3):
    std::cout << "    PID controller for fast time scale, based on order of RK "
              << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (4):
    std::cout
      << "    ExpGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (5):
    std::cout
      << "    ImpGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  case (6):
    std::cout
      << "    ImExGus controller for fast time scale, based on order of RK "
      << ((opts.fast_pq == 1) ? "method\n" : "embedding\n");
    std::cout << "    rtol = " << opts.rtol << ", atol = " << opts.atol << "\n";
    break;
  }
}

//---- end of file ----//
