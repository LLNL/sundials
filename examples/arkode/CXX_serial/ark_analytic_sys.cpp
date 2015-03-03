/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem with analytical 
 * solution,
 *    dy/dt = A*y
 * where A = V*D*Vi, 
 *      V = [1 -1 1; -1 2 1; 0 -1 2];
 *      Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1];
 *      D = [-0.5 0 0; 0 -0.1 0; 0 0 lam];
 * where lam is a large negative number. The analytical solution to
 * this problem is 
 *   Y(t) = V*exp(D*t)*Vi*Y0
 * for t in the interval [0.0, 0.05], with initial condition: 
 * y(0) = [1,1,1]'.
 * 
 * The stiffness of the problem is directly proportional to the 
 * value of "lamda".  The value of lamda should be negative to 
 * result in a well-posed ODE; for values with magnitude larger than
 * 100 the problem becomes quite stiff.
 * 
 * In this example, we choose lamda = -100.
 * 
 * This program solves the problem with the DIRK method,
 * Newton iteration with the ARKDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*/

// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <arkode/arkode.h>            // prototypes for ARKode fcts., consts.
#include <nvector/nvector_serial.h>   // serial N_Vector types, fcts., macros
#include <arkode/arkode_dense.h>      // prototype for ARKDense solver
#include <sundials/sundials_dense.h>  // defs. of DlsMat and DENSE_ELEM
#include <sundials/sundials_types.h>  // def. of type 'realtype'

using namespace std;

// User-supplied Functions Called by the Solver
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to perform matrix-matrix product
static int dense_MM(DlsMat A, DlsMat B, DlsMat C);

// Private function to check function return values
static int check_flag(void *flagvalue, const string funcname, int opt);

// Main Program
int main()
{
  // general problem parameters
  realtype T0 = RCONST(0.0);       // initial time
  realtype Tf = RCONST(0.05);      // final time
  realtype dTout = RCONST(0.005);  // time between outputs
  long int NEQ = 3;                // number of dependent vars.
  realtype reltol = 1.0e-6;        // tolerances
  realtype abstol = 1.0e-10;
  realtype lamda  = -100.0;        // stiffness parameter

  // general problem variables
  int flag;                      // reusable error-checking flag
  N_Vector y = NULL;             // empty vector for storing solution
  void *arkode_mem = NULL;       // empty ARKode memory structure

  // Initial problem output
  cout << "\nAnalytical ODE test problem:\n";
  cout << "    lamda = " << lamda << "\n";
  cout << "   reltol = " << reltol << "\n";
  cout << "   abstol = " << abstol << "\n\n";

  // Initialize data structures
  y = N_VNew_Serial(NEQ);         // Create serial vector solution
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = 1.0;            // Specify initial condition
  NV_Ith_S(y,1) = 1.0;
  NV_Ith_S(y,2) = 1.0;
  arkode_mem = ARKodeCreate();    // Create the solver memory
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;

  /* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  // Set routines
  flag = ARKodeSetUserData(arkode_mem, (void *) &lamda);   // Pass lamda to user functions
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);   // Specify tolerances
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  // Linear solver specification
  flag = ARKDense(arkode_mem, NEQ);              // Specify dense linear solver
  if (check_flag(&flag, "ARKDense", 1)) return 1;
  flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);   // Set Jacobian routine
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return 1;

  // Specify linearly implicit RHS, with non-time-dependent Jacobian
  flag = ARKodeSetLinear(arkode_mem, 0);
  if (check_flag(&flag, "ARKodeSetLinear", 1)) return 1;

  // Open output stream for results, output comment line
  FILE *UFID = fopen("solution.txt","w");
  fprintf(UFID,"# t y1 y2 y3\n");

  // output initial condition to disk 
  fprintf(UFID," %.16e %.16e %.16e %.16e\n", 
	  T0, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));  

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  realtype t = T0;
  realtype tout = T0+dTout;
  cout << "      t        y0        y1        y2\n";
  cout << "   --------------------------------------\n";
  while (Tf - t > 1.0e-15) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);       // call integrator
    if (check_flag(&flag, "ARKode", 1)) break;
    printf("  %8.4f  %8.5f  %8.5f  %8.5f\n",                  // access/print solution
           t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
    fprintf(UFID," %.16e %.16e %.16e %.16e\n", 
	    t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));  
    if (flag >= 0) {                                          // successful solve: update time
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                  // unsuccessful solve: break
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }
  }
  cout << "   --------------------------------------\n";
  fclose(UFID);

  // Print some final statistics
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
  flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKDlsGetNumJacEvals", 1);
  flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKDlsGetNumRhsEvals", 1);

  cout << "\nFinal Solver Statistics:\n";
  cout << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
  cout << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
  cout << "   Total linear solver setups = " << nsetups << "\n";
  cout << "   Total RHS evals for setting up the linear system = " << nfeLS << "\n";
  cout << "   Total number of Jacobian evaluations = " << nje << "\n";
  cout << "   Total number of Newton iterations = " << nni << "\n";
  cout << "   Total number of linear solver convergence failures = " << ncfn << "\n";
  cout << "   Total number of error test failures = " << netf << "\n\n";

  // Clean up and return with successful completion
  N_VDestroy_Serial(y);        // Free y vector
  ARKodeFree(&arkode_mem);     // Free integrator memory
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

// f routine to compute the ODE RHS function f(t,y).
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);
  realtype y2 = NV_Ith_S(y,2);
  realtype yd0, yd1, yd2;

  // fill in the RHS function: f(t,y) = V*D*Vi*y
  yd0 = 0.25*(5.0*y0 + 1.0*y1 - 3.0*y2);     // yd = Vi*y
  yd1 = 0.25*(2.0*y0 + 2.0*y1 - 2.0*y2);
  yd2 = 0.25*(1.0*y0 + 1.0*y1 + 1.0*y2);
  y0  = -0.5*yd0;                            //  y = D*yd
  y1  = -0.1*yd1;
  y2  =  lam*yd2;
  yd0 =  1.0*y0 - 1.0*y1 + 1.0*y2;           // yd = V*y
  yd1 = -1.0*y0 + 2.0*y1 + 1.0*y2;
  yd2 =  0.0*y0 - 1.0*y1 + 2.0*y2;
  NV_Ith_S(ydot,0) = yd0;
  NV_Ith_S(ydot,1) = yd1;
  NV_Ith_S(ydot,2) = yd2;

  return 0;                                  // Return with success
}

// Jacobian routine to compute J(t,y) = df/dy.
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  DlsMat V  = NewDenseMat(3,3);               // create temporary DlsMat objects
  DlsMat D  = NewDenseMat(3,3);
  DlsMat Vi = NewDenseMat(3,3);

  DenseScale(0.0, V);     // initialize temporary matrices to zero
  DenseScale(0.0, D);
  DenseScale(0.0, Vi);

  // Fill in temporary matrices:
  //    V = [1 -1 1; -1 2 1; 0 -1 2]
  DENSE_ELEM(V,0,0) =  1.0;
  DENSE_ELEM(V,0,1) = -1.0;
  DENSE_ELEM(V,0,2) =  1.0;
  DENSE_ELEM(V,1,0) = -1.0;
  DENSE_ELEM(V,1,1) =  2.0;
  DENSE_ELEM(V,1,2) =  1.0;
  DENSE_ELEM(V,2,0) =  0.0;
  DENSE_ELEM(V,2,1) = -1.0;
  DENSE_ELEM(V,2,2) =  2.0;

  //    Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
  DENSE_ELEM(Vi,0,0) =  0.25*5.0;
  DENSE_ELEM(Vi,0,1) =  0.25*1.0;
  DENSE_ELEM(Vi,0,2) = -0.25*3.0;
  DENSE_ELEM(Vi,1,0) =  0.25*2.0;
  DENSE_ELEM(Vi,1,1) =  0.25*2.0;
  DENSE_ELEM(Vi,1,2) = -0.25*2.0;
  DENSE_ELEM(Vi,2,0) =  0.25*1.0;
  DENSE_ELEM(Vi,2,1) =  0.25*1.0;
  DENSE_ELEM(Vi,2,2) =  0.25*1.0;

  //    D = [-0.5 0 0; 0 -0.1 0; 0 0 lam]
  DENSE_ELEM(D,0,0) = -0.5;
  DENSE_ELEM(D,1,1) = -0.1;
  DENSE_ELEM(D,2,2) = lam;

  // Compute J = V*D*Vi
  if (dense_MM(D,Vi,J) != 0) {     // J = D*Vi
    cerr << "matmul error\n";
    return 1;
  }
  if (dense_MM(V,J,D) != 0) {      // D = V*J [= V*D*Vi]
    cerr << "matmul error\n";
    return 1;
  }
  DenseCopy(D, J);                 // J = D [= V*D*Vi]

  return 0;                        // Return with success
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

// DlsMat matrix-multiply utility routine: C = A*B.
static int dense_MM(DlsMat A, DlsMat B, DlsMat C)
{
  // check for legal dimensions
  if ((A->N != B->M) || (C->M != A->M) || (C->N != B->N)) {
    cerr << "\n matmul error: dimension mismatch\n\n";
    return 1;
  }

  realtype **adata = A->cols;     // access data and extents
  realtype **bdata = B->cols;
  realtype **cdata = C->cols;
  long int m = C->M;
  long int n = C->N;
  long int l = A->N;
  int i, j, k;
  DenseScale(0.0, C);             // initialize output

  // perform multiply (not optimal, but fine for 3x3 matrices)
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      for (k=0; k<l; k++)
     cdata[i][j] += adata[i][k] * bdata[k][j];

  return 0;                       // Return with success
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer  
*/
static int check_flag(void *flagvalue, const string funcname, int opt)
{
  int *errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    cerr << "\nSUNDIALS_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  // Check if flag < 0
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      cerr << "\nSUNDIALS_ERROR: " << funcname << " failed with flag = " << *errflag << "\n\n";
      return 1; 
    }
  }
  
  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL) {
    cerr << "\nMEMORY_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  return 0;
}



//---- end of file ----
