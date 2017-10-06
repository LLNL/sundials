/*
 * -----------------------------------------------------------------
 * $Revision: 4882 $
 * $Date: 2016-09-01 16:05:08 -0700 (Thu, 01 Sep 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on PETSc TS example 15 and a SUNDIALS example by 
 * Allan Taylor, Alan Hindmarsh and Radu Serban
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, PETSc vector, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver IDASpgmr.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MX x MY
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = MX * MY. Here MX = MY = 10.
 *
 * The system is actually implemented on submeshes, processor by
 * processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
 * processors.
 *
 * The system is solved with IDA using the Krylov linear solver
 * IDASPGMR. The preconditioner uses the diagonal elements of the
 * Jacobian only. Routines for preconditioning, required by
 * IDASPGMR, are supplied here. The constraints u >= 0 are posed
 * for all components. Local error testing on the boundary values
 * is suppressed. Output is taken at t = 0, .01, .02, .04,
 * ..., 10.24.
 * 
 * The example uses PETSc vector and data management functions to 
 * generate mesh and handle MPI communication.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <ida/ida_spgmr.h>
#include <nvector/nvector_petsc.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <mpi.h>
#include <petscdm.h>
#include <petscdmda.h>

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define NOUT         11             /* Number of output times */

#define NPEX         2              /* No. PEs in x direction of PE array */
#define NPEY         2              /* No. PEs in y direction of PE array */
                                    /* Total no. PEs = NPEX*NPEY */
#define MXSUB        5              /* No. x points per subgrid */
#define MYSUB        5              /* No. y points per subgrid */

#define MX           (NPEX*MXSUB)   /* MX = number of x mesh points */
#define MY           (NPEY*MYSUB)   /* MY = number of y mesh points */
                                    /* Spatial mesh is MX by MY */

typedef struct {  
  N_Vector pp;    /* vector of diagonal preconditioner elements */
  DM       da;    /* PETSc data management object */
} *UserData;

/* User-supplied residual function */

int resHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr, void *user_data);

/* User-supplied preconditioner routines */

int PsetupHeat(realtype tt, 
               N_Vector yy, N_Vector yp, N_Vector rr, 
               realtype c_j, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int PsolveHeat(realtype tt, 
               N_Vector uu, N_Vector up, N_Vector rr, 
               N_Vector rvec, N_Vector zvec,
               realtype c_j, realtype delta, void *user_data, 
               N_Vector tmp);

/* Private function to check function return values */

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, void *user_data);

static void PrintHeader(long int Neq, realtype rtol, realtype atol);

static void PrintOutput(int id, void *mem, realtype t, N_Vector uu);

static void PrintFinalStats(void *mem);

static int check_flag(void *flagvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  MPI_Comm comm;
  void *mem;
  UserData data;
  int iout, thispe, ier, npes;
  long int Neq;
  realtype rtol, atol, t0, t1, tout, tret;
  N_Vector uu, up, constraints, id, res;
  PetscErrorCode ierr;                  /* PETSc error code  */
  Vec uvec;

  mem = NULL;
  data = NULL;
  uu = up = constraints = id = res = NULL;

  /* Get processor number and total number of pe's. */

  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &thispe);
  
  if (npes != NPEX*NPEY) {
    if (thispe == 0)
      fprintf(stderr, 
              "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n", 
              npes,NPEX*NPEY);
    MPI_Finalize();
    return(1);
  }
  
  /* Set global vector length Neq. */
  Neq = MX * MY;
  
  /* Initialize PETSc */
  ierr = PetscInitializeNoArguments();
  CHKERRQ(ierr);

  /* Allocate and initialize the data structure and N-vectors. */
  data = (UserData) malloc(sizeof *data);
  if(check_flag((void *)data, "malloc", 2, thispe)) 
    MPI_Abort(comm, 1);
  data->pp = NULL;
  data->da = NULL;
  
  /* Create the object that will manage the communication of 2D data */
  ierr = DMDACreate2d(comm, 
                      DM_BOUNDARY_NONE,  /* NONE, PERIODIC, GHOSTED */
                      DM_BOUNDARY_NONE,
                      DMDA_STENCIL_STAR, /* STAR, BOX */
                      MX,                
                      MY,                
                      NPEX,              /* Set numbers or use PETSC_DECIDE */
                      NPEY,              
                      1,                 /* degrees of freedom per node */
                      1,                 /* stencil width */
                      NULL,              /* number of nodes per cell in x */
                      NULL,              /* number of nodes per cell in y */
                      &(data->da));
  CHKERRQ(ierr);

  /* Create PETSc vector */
  ierr = DMCreateGlobalVector(data->da, &uvec);
  CHKERRQ(ierr);

  /* Make N_Vector wrapper for uvec */
  uu = N_VMake_Petsc(&uvec);
  if(check_flag((void *)uu, "N_VNew_Petsc", 0, thispe)) 
    MPI_Abort(comm, 1);

  up = N_VClone(uu);
  if(check_flag((void *)up, "N_VNew_Petsc", 0, thispe)) 
    MPI_Abort(comm, 1);

  res = N_VClone(uu);
  if(check_flag((void *)res, "N_VNew_Petsc", 0, thispe)) 
    MPI_Abort(comm, 1);

  constraints = N_VClone(uu);
  if(check_flag((void *)constraints, "N_VNew_Petsc", 0, thispe)) 
    MPI_Abort(comm, 1);

  id = N_VClone(uu);
  if(check_flag((void *)id, "N_VNew_Petsc", 0, thispe)) 
    MPI_Abort(comm, 1);

  /* An N-vector to hold preconditioner. */
  data->pp = N_VClone(uu); 
  if(check_flag((void *)data->pp, "N_VNew_Petsc", 0, thispe)) 
    MPI_Abort(comm, 1);

  /* Initialize the uu, up, id, and res profiles. */

  SetInitialProfile(uu, up, id, res, data);
  
  /* Set constraints to all 1's for nonnegative solution values. */

  N_VConst(ONE, constraints);
  
  t0 = ZERO; t1 = RCONST(0.01);
  
  /* Scalar relative and absolute tolerance. */

  rtol = ZERO;
  atol = RCONST(1.0e-3);

  /* Call IDACreate and IDAMalloc to initialize solution. */

  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0, thispe)) MPI_Abort(comm, 1);

  ier = IDASetUserData(mem, data);
  if(check_flag(&ier, "IDASetUserData", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASetSuppressAlg(mem, TRUE);
  if(check_flag(&ier, "IDASetSuppressAlg", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASetId(mem, id);
  if(check_flag(&ier, "IDASetId", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASetConstraints(mem, constraints);
  if(check_flag(&ier, "IDASetConstraints", 1, thispe)) MPI_Abort(comm, 1);
  N_VDestroy_Petsc(constraints);  

  ier = IDAInit(mem, resHeat, t0, uu, up);
  if(check_flag(&ier, "IDAInit", 1, thispe)) MPI_Abort(comm, 1);
  
  ier = IDASStolerances(mem, rtol, atol);
  if(check_flag(&ier, "IDASStolerances", 1, thispe)) MPI_Abort(comm, 1);

  /* Call IDASpgmr to specify the linear solver. */

  ier = IDASpgmr(mem, 0);
  if(check_flag(&ier, "IDASpgmr", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASpilsSetPreconditioner(mem, PsetupHeat, PsolveHeat);
  if(check_flag(&ier, "IDASpilsSetPreconditioner", 1, thispe)) MPI_Abort(comm, 1);

  /* Print output heading (on processor 0 only) and intial solution  */
  
  if (thispe == 0) PrintHeader(Neq, rtol, atol);
  PrintOutput(thispe, mem, t0, uu); 
  
  /* Loop over tout, call IDASolve, print output. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    ier = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_flag(&ier, "IDASolve", 1, thispe)) MPI_Abort(comm, 1);

    PrintOutput(thispe, mem, tret, uu);

  }
  
  /* Print remaining counters. */

  if (thispe == 0) PrintFinalStats(mem);

  /* Free memory */

  IDAFree(&mem);

  N_VDestroy_Petsc(id);
  N_VDestroy_Petsc(res);
  N_VDestroy_Petsc(up);
  N_VDestroy_Petsc(uu);

  N_VDestroy_Petsc(data->pp);
  ierr = DMDestroy(&data->da);
  CHKERRQ(ierr);
  free(data);
  
  ierr = VecDestroy(&uvec);
  CHKERRQ(ierr);

  ierr = PetscFinalize();
  CHKERRQ(ierr);

  MPI_Finalize();

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resHeat: heat equation system residual function                       
 * This uses 5-point central differencing on the interior points, and    
 * includes algebraic equations for the boundary values.                 
 * So for each interior point, the residual component has the form       
 *    res_i = u'_i - (central difference)_i                              
 * while for each boundary point, it is res_i = u_i. 
 *                    
 */

int resHeat(realtype tt, 
            N_Vector uu, N_Vector up, N_Vector rr, 
            void *user_data)
{
  PetscErrorCode ierr;
  UserData       data = (UserData) user_data;
  DM             da   = (DM) data->da;
  PetscInt       i, j, Mx, My, xs, ys, xm, ym;
  PetscReal      hx, hy, sx, sy;
  PetscScalar    u, uxx, uyy, **uarray, **f, **udot;
  Vec localU;
  Vec *U    = N_VGetVector_Petsc(uu);
  Vec *Udot = N_VGetVector_Petsc(up);
  Vec *F    = N_VGetVector_Petsc(rr);

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(da, &localU);
  CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,
                     PETSC_IGNORE,
                     &Mx,          /* Get global grid x size */
                     &My,          /* Get global grid y size */
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE);

  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  ierr = DMGlobalToLocalBegin(da, *U, INSERT_VALUES, localU);
  CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da, *U, INSERT_VALUES, localU);
  CHKERRQ(ierr);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArrayRead(da, localU, &uarray);
  CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, *F, &f);
  CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, *Udot, &udot);
  CHKERRQ(ierr);

  /* Get local grid boundaries */
  ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      /* Boundary conditions */
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i] = uarray[j][i]; /* F = U */
      } else { /* Interior */
        u = uarray[j][i];
        /* 5-point stencil */
        uxx = (-2.0*u + uarray[j][i-1] + uarray[j][i+1]);
        uyy = (-2.0*u + uarray[j-1][i] + uarray[j+1][i]);
        f[j][i] = udot[j][i] - (uxx*sx + uyy*sy);
      }
    }
  }

  /* Restore vectors */
  ierr = DMDAVecRestoreArrayRead(da, localU, &uarray);
  CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, *F, &f);
  CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, *Udot, &udot);
  CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da, &localU);
  CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


/*
 * PsetupHeat: setup for diagonal preconditioner for heatsk.    
 *                                                                 
 * The optional user-supplied functions PsetupHeat and          
 * PsolveHeat together must define the left preconditoner        
 * matrix P approximating the system Jacobian matrix               
 *                   J = dF/du + cj*dF/du'                         
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   
 * systems P z = r.   This is done in this case by keeping only    
 * the diagonal elements of the J matrix above, storing them as    
 * inverses in a vector pp, when computed in PsetupHeat, for    
 * subsequent use in PsolveHeat.                                 
 *                                                                 
 * In this instance, only cj and data (user data structure, with    
 * pp etc.) are used from the PsetupHeat argument list.         
 *
 */

int PsetupHeat(realtype tt, 
               N_Vector yy, N_Vector yp, N_Vector rr, 
               realtype c_j, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  PetscErrorCode ierr;
  PetscInt    i, j, Mx, My, xs, ys, xm, ym;
  PetscReal   hx,hy,sx,sy;
  PetscScalar pelinv;
  PetscScalar **ppv;
  UserData data = (UserData) user_data;
  DM da = (DM) data->da;
  Vec *ppvec = N_VGetVector_Petsc((data->pp));
  
  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,
                     PETSC_IGNORE,
                     &Mx,          /* Get global grid x size */
                     &My,          /* Get global grid y size */
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE);
  ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  CHKERRQ(ierr);

  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);
  
  pelinv = ONE/(2.0*(sx + sy) + c_j);

  ierr = DMDAVecGetArray(da, *ppvec, &ppv);
  CHKERRQ(ierr);

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        ppv[j][i] = ONE;    /* Boundary */
      } else {
        ppv[j][i] = pelinv; /* Interior */
      }
    }
  }

  ierr = DMDAVecRestoreArray(da, *ppvec, &ppv);
  CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


/*
 * PsolveHeat: solve preconditioner linear system.              
 * This routine multiplies the input vector rvec by the vector pp 
 * containing the inverse diagonal Jacobian elements (previously  
 * computed in PsetupHeat), returning the result in zvec.      
 */

int PsolveHeat(realtype tt, 
               N_Vector uu, N_Vector up, N_Vector rr, 
               N_Vector rvec, N_Vector zvec,
               realtype c_j, realtype delta, void *user_data, 
               N_Vector tmp)
{
  UserData data;

  data = (UserData) user_data;
  
  N_VProd(data->pp, rvec, zvec);

  return(0);

}


/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */


/*
 * SetInitialProfile sets the initial values for the problem. 
 */

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, void *user_data)
{
  UserData       data = (UserData) user_data;
  DM             da   = data->da;
  PetscErrorCode ierr;
  PetscInt       i, j, xs, ys, xm, ym, Mx, My;
  PetscScalar    **u;
  PetscScalar    **iddat;
  PetscReal      hx, hy, x, y;
  Vec *U     = N_VGetVector_Petsc(uu);
  Vec *idvec = N_VGetVector_Petsc(id);

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,
                     PETSC_IGNORE,
                     &Mx,
                     &My,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE);

  hx = 1.0/(PetscReal)(Mx-1);
  hy = 1.0/(PetscReal)(My-1);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArray(da, *U, &u);
  CHKERRQ(ierr);

  /* Get pointers to differentiable variable IDs */
  ierr = DMDAVecGetArray(da, *idvec, &iddat);
  CHKERRQ(ierr);

  /* Get local grid boundaries */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);
  CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (j=ys; j<ys+ym; j++) {
    y = j*hy;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx;
      u[j][i] = 16.0 * x*(1.0 - x) * y*(1.0 - y);
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        iddat[j][i] = 0.0; /* algebraic variables on the boundary */
      } else { 
        iddat[j][i] = 1.0; /* differential variables in the interior */
      }
    }
  }

  /* Restore vectors */
  ierr = DMDAVecRestoreArray(da, *U, &u);
  CHKERRQ(ierr);

   /* Restore vectors */
  ierr = DMDAVecRestoreArray(da, *idvec, &iddat);
  CHKERRQ(ierr);

  /* Initialize up. */
  N_VConst(ZERO, up);    /* Initially set up = 0. */
  
  /* resHeat sets res to negative of ODE RHS values at interior points. */
  resHeat(ZERO, uu, up, res, data);
  
  /* Copy -res into up to get correct initial up values. */
  N_VScale(-ONE, res, up);
  
  PetscFunctionReturn(0);
}


/*
 * Print first lines of output and table heading
 */

static void PrintHeader(long int Neq, realtype rtol, realtype atol)
{ 
  printf("\nidaHeat2D_kry_petsc: Heat equation, parallel example problem for IDA\n");
  printf("            Discretized heat equation on 2D unit square.\n");
  printf("            Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("            Mesh dimensions: %d x %d", MX, MY);
  printf("        Total system size: %ld\n\n", Neq);
  printf("Subgrid dimensions: %d x %d", MXSUB, MYSUB);
  printf("        Processor array: %d x %d\n", NPEX, NPEY);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");
  printf("SUPPRESSALG = TRUE to suppress local error testing on ");
  printf("all boundary components. \n");
  printf("Linear solver: IDASPGMR  ");
  printf("Preconditioner: diagonal elements only.\n"); 
  printf("This example uses PETSc vector.\n");
  
  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nli   nre   nreLS    h      npe nps\n");
  printf("----------------------------------------------------------------------\n");
}

/*
 * PrintOutput: print max norm of solution and current solver statistics
 */

static void PrintOutput(int id, void *mem, realtype t, N_Vector uu)
{
  realtype hused, umax;
  long int nst, nni, nje, nre, nreLS, nli, npe, nps;
  int kused, ier;

  umax = N_VMaxNorm(uu);

  if (id == 0) {

    ier = IDAGetLastOrder(mem, &kused);
    check_flag(&ier, "IDAGetLastOrder", 1, id);
    ier = IDAGetNumSteps(mem, &nst);
    check_flag(&ier, "IDAGetNumSteps", 1, id);
    ier = IDAGetNumNonlinSolvIters(mem, &nni);
    check_flag(&ier, "IDAGetNumNonlinSolvIters", 1, id);
    ier = IDAGetNumResEvals(mem, &nre);
    check_flag(&ier, "IDAGetNumResEvals", 1, id);
    ier = IDAGetLastStep(mem, &hused);
    check_flag(&ier, "IDAGetLastStep", 1, id);
    ier = IDASpilsGetNumJtimesEvals(mem, &nje);
    check_flag(&ier, "IDASpilsGetNumJtimesEvals", 1, id);
    ier = IDASpilsGetNumLinIters(mem, &nli);
    check_flag(&ier, "IDASpilsGetNumLinIters", 1, id);
    ier = IDASpilsGetNumResEvals(mem, &nreLS);
    check_flag(&ier, "IDASpilsGetNumResEvals", 1, id);
    ier = IDASpilsGetNumPrecEvals(mem, &npe);
    check_flag(&ier, "IDASpilsGetPrecEvals", 1, id);
    ier = IDASpilsGetNumPrecSolves(mem, &nps);
    check_flag(&ier, "IDASpilsGetNumPrecSolves", 1, id);

#if defined(SUNDIALS_EXTENDED_PRECISION)  
    printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2Le  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#elif defined(SUNDIALS_DOUBLE_PRECISION)  
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#else
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#endif

  }
}

/*
 * Print some final integrator statistics
 */

static void PrintFinalStats(void *mem)
{
  long int netf, ncfn, ncfl;

  IDAGetNumErrTestFails(mem, &netf);
  IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  IDASpilsGetNumConvFails(mem, &ncfl);

  printf("\nError test failures            = %ld\n", netf);
  printf("Nonlinear convergence failures = %ld\n", ncfn);
  printf("Linear convergence failures    = %ld\n", ncfl);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt, int id)
{
  int *errflag;

  if (opt == 0 && flagvalue == NULL) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n", 
            id, funcname);
    return(1); 
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR(%d): %s() failed with flag = %d\n\n", 
              id, funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n", 
            id, funcname);
    return(1); 
  }

  return(0);
}
