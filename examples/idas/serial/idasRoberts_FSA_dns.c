/* -----------------------------------------------------------------
 * Programmer(s): Cosmin Petra and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * This simple example problem for IDA, due to Robertson,
 * is from chemical kinetics, and consists of the following three
 * equations:
 *
 *      dy1/dt = -p1*y1 + p2*y2*y3
 *      dy2/dt = p1*y1 - p2*y2*y3 - p3*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7
 *
 * Optionally, IDAS can compute sensitivities with respect to the
 * problem parameters p1, p2, and p3.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of two sensitivity methods (SIMULTANEOUS and STAGGERED can be
 * used and sensitivities may be included in the error test or not
 *(error control set on SUNTRUE or SUNFALSE, respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % idasRoberts_FSA_dns -nosensi
 * If sensitivities are to be computed:
 *    % idasRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <idas/idas.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* defs. of SUNRabs, SUNRexp, etc.      */

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */

/* Problem Constants */

#define NEQ   3             /* number of equations  */
#define T0    RCONST(0.0)   /* initial time */
#define T1    RCONST(0.4)   /* first output time */
#define TMULT RCONST(10.0)  /* output time factor */
#define NOUT  12            /* number of output times */

#define NP    3             /* number of problem parameters */
#define NS    3             /* number of sensitivities computed */

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)

/* Type : UserData */

typedef struct {
  realtype p[3];           /* problem parameters */
  realtype coef;
} *UserData;

/* Prototypes of functions by IDAS */

static int res(realtype t, N_Vector y, N_Vector yp, N_Vector resval, void *user_data);

static int resS(int Ns, realtype t,
                N_Vector y, N_Vector yp, N_Vector resval,
                N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int rhsQ(realtype tres, N_Vector yy, N_Vector yp,
                 N_Vector rrQ, void *user_data);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth,
                        booleantype *err_con);
static void WrongArgs(char *name);

static void PrintIC(N_Vector y, N_Vector yp);
static void PrintSensIC(N_Vector y, N_Vector yp, N_Vector* yS, N_Vector* ypS);

static void PrintOutput(void *ida_mem, realtype t, N_Vector u);
static void PrintSensOutput(N_Vector *uS);

static int check_retval(void *returnvalue, const char *funcname, int opt);
/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  SUNContext ctx;
  void *ida_mem;
  SUNMatrix A;
  SUNLinearSolver LS;
  UserData data;
  realtype reltol, t, tout;
  N_Vector y, yp, abstol, id;
  int iout, retval;
  FILE* FID;
  char fname[256];

  realtype pbar[NS];
  int is;
  N_Vector *yS, *ypS;
  booleantype sensi, err_con;
  int sensi_meth;

  N_Vector yQ, *yQS;

  ida_mem = NULL;
  A       = NULL;
  LS      = NULL;
  data    = NULL;
  y       = NULL;
  yS      = NULL;
  ypS     = NULL;
  yQS     = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* Create the SUNDIALS context object for this simulation */
  retval = SUNContext_Create(NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  if (check_retval((void *)data, "malloc", 2)) return(1);
  data->p[0] = RCONST(0.040);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);
  data->coef = RCONST(0.5);

  /* Initial conditions */
  y = N_VNew_Serial(NEQ, ctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  Ith(y,1) = ONE;
  Ith(y,2) = ZERO;
  Ith(y,3) = ZERO;

  yp = N_VClone(y);
  if(check_retval((void *)yp, "N_VNew_Serial", 0)) return(1);

  /* These initial conditions are NOT consistent. See IDACalcIC below. */
  Ith(yp,1) = RCONST(0.1);
  Ith(yp,2) = ZERO;
  Ith(yp,3) = ZERO;

  /* Create IDAS object */
  ida_mem = IDACreate(ctx);
  if (check_retval((void *)ida_mem, "IDACreate", 0)) return(1);

  /* Allocate space for IDAS */
  retval = IDAInit(ida_mem, res, T0, y, yp);
  if (check_retval(&retval, "IDAInit", 1)) return(1);

  /* Specify scalar relative tol. and vector absolute tol. */
  reltol = RCONST(1.0e-6);
  abstol = N_VClone(y);
  Ith(abstol,1) = RCONST(1.0e-8);
  Ith(abstol,2) = RCONST(1.0e-14);
  Ith(abstol,3) = RCONST(1.0e-6);
  retval = IDASVtolerances(ida_mem, reltol, abstol);
  if (check_retval(&retval, "IDASVtolerances", 1)) return(1);

  /* Set ID vector */
  id = N_VClone(y);
  Ith(id,1) = 1.0;
  Ith(id,2) = 1.0;
  Ith(id,3) = 0.0;
  retval = IDASetId(ida_mem, id);
  if (check_retval(&retval, "IDASetId", 1)) return(1);

  /* Attach user data */
  retval = IDASetUserData(ida_mem, data);
  if (check_retval(&retval, "IDASetUserData", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ, ctx);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(y, A, ctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  printf("\n3-species chemical kinetics problem\n");

  /* Sensitivity-related settings */
  if (sensi) {

    pbar[0] = data->p[0];
    pbar[1] = data->p[1];
    pbar[2] = data->p[2];

    yS = N_VCloneVectorArray(NS, y);
    if (check_retval((void *)yS, "N_VCloneVectorArray", 0)) return(1);
    for (is=0;is<NS;is++) N_VConst(ZERO, yS[is]);

    ypS = N_VCloneVectorArray(NS, y);
    if (check_retval((void *)ypS, "N_VCloneVectorArray", 0)) return(1);
    for (is=0;is<NS;is++) N_VConst(ZERO, ypS[is]);

    /*
    * Only non-zero sensitivity I.C. are ypS[0]:
    * - Ith(ypS[0],1) = -ONE;
    * - Ith(ypS[0],2) =  ONE;
    *
    * They are not set. IDACalcIC also computes consistent IC for sensitivities.
    */

    retval = IDASensInit(ida_mem, NS, sensi_meth, resS, yS, ypS);
    if(check_retval(&retval, "IDASensInit", 1)) return(1);

    retval = IDASensEEtolerances(ida_mem);
    if(check_retval(&retval, "IDASensEEtolerances", 1)) return(1);

    retval = IDASetSensErrCon(ida_mem, err_con);
    if (check_retval(&retval, "IDASetSensErrCon", 1)) return(1);

    retval = IDASetSensParams(ida_mem, data->p, pbar, NULL);
    if (check_retval(&retval, "IDASetSensParams", 1)) return(1);

    printf("Sensitivity: YES ");
    if(sensi_meth == IDA_SIMULTANEOUS)
      printf("( SIMULTANEOUS +");
    else
      printf("( STAGGERED +");
    if(err_con) printf(" FULL ERROR CONTROL )");
    else        printf(" PARTIAL ERROR CONTROL )");

  } else {

    printf("Sensitivity: NO ");

  }

  /*----------------------------------------------------------
   *               Q U A D R A T U R E S
   * ---------------------------------------------------------*/
  yQ = N_VNew_Serial(2, ctx);

  Ith(yQ,1) = 0;
  Ith(yQ,2) = 0;

  IDAQuadInit(ida_mem, rhsQ, yQ);

  if (sensi) {
    yQS = N_VCloneVectorArray(NS, yQ);
    for (is=0;is<NS;is++) N_VConst(ZERO, yQS[is]);

    IDAQuadSensInit(ida_mem, NULL, yQS);
  }

  /* Call IDACalcIC to compute consistent initial conditions. If sensitivity is
     enabled, this function also try to find consistent IC for the sensitivities. */

  retval = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, T1);;
  if (check_retval(&retval, "IDACalcIC", 1)) return(1);

  retval = IDAGetConsistentIC(ida_mem, y, yp);
  if (check_retval(&retval, "IDAGetConsistentIC", 1)) return(1);

  PrintIC(y, yp);

  if(sensi) {
      IDAGetSensConsistentIC(ida_mem, yS, ypS);
      PrintSensIC(y, yp, yS, ypS);
    }

  /* In loop over output points, call IDA, print results, test for error */

  printf("\n\n");
  printf("===========================================");
  printf("============================\n");
  printf("     T     Q       H      NST           y1");
  printf("           y2           y3    \n");
  printf("===========================================");
  printf("============================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {

    retval = IDASolve(ida_mem, tout, &t, y, yp, IDA_NORMAL);
    if (check_retval(&retval, "IDASolve", 1)) break;

    PrintOutput(ida_mem, t, y);

    if (sensi) {
      retval = IDAGetSens(ida_mem, &t, yS);
      if (check_retval(&retval, "IDAGetSens", 1)) break;
      PrintSensOutput(yS);
    }
    printf("-----------------------------------------");
    printf("------------------------------\n");

  }

  printf("\nQuadrature:\n");
  IDAGetQuad(ida_mem, &t, yQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("G:      %10.4Le\n", Ith(yQ,1));
#else
  printf("G:      %10.4e\n", Ith(yQ,1));
#endif

  if(sensi) {
    IDAGetQuadSens(ida_mem, &t, yQS);
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("\nSensitivities at t=%Lg:\n",t);
    printf("dG/dp1: %11.4Le\n", Ith(yQS[0], 1));
    printf("dG/dp1: %11.4Le\n", Ith(yQS[1], 1));
    printf("dG/dp1: %11.4Le\n", Ith(yQS[2], 1));
#else
    printf("\nSensitivities at t=%g:\n",t);
    printf("dG/dp1: %11.4e\n", Ith(yQS[0], 1));
    printf("dG/dp1: %11.4e\n", Ith(yQS[1], 1));
    printf("dG/dp1: %11.4e\n", Ith(yQS[2], 1));
#endif
  }

  /* Print final statistics to the screen */
  printf("\nFinal Statistics:\n");
  retval = IDAPrintAllStats(ida_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

  /* Print final statistics to a file in CSV format */
  strcpy(fname, "idasRoberts_FSA_dns_stats");
  if (sensi)
  {
    if(sensi_meth == IDA_SIMULTANEOUS)
      strcat(fname, "_-sensi_sim");
    else
      strcat(fname, "_-sensi_stg");
    if(err_con)
      strcat(fname, "_t");
    else
      strcat(fname, "_f");
  }
  strcat(fname, ".csv");
  FID = fopen(fname, "w");
  retval = IDAPrintAllStats(ida_mem, FID, SUN_OUTPUTFORMAT_CSV);
  fclose(FID);

  /* Free memory */
  N_VDestroy(y);
  N_VDestroy(yp);
  N_VDestroy(abstol);
  N_VDestroy(id);
  N_VDestroy(yQ);
  if (sensi) {
    N_VDestroyVectorArray(yS,  NS);
    N_VDestroyVectorArray(ypS, NS);
    N_VDestroyVectorArray(yQS, NS);
  }
  free(data);
  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNContext_Free(&ctx);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDAS
 *--------------------------------------------------------------------
 */

/*
 * Residual routine. Compute F(t,y,y',p).
 */
static int res(realtype t, N_Vector yy, N_Vector yp, N_Vector resval, void *user_data)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype yp1, yp2;

  data = (UserData) user_data;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);

  yp1 = Ith(yp,1);
  yp2 = Ith(yp,2);

  Ith(resval,1) = yp1 + p1*y1 - p2*y2*y3;
  Ith(resval,2) = yp2 - p1*y1 + p2*y2*y3 + p3*y2*y2;
  Ith(resval,3) = y1 + y2 + y3 - ONE;

  return(0);
}


/*
 * resS routine. Compute sensitivity r.h.s.
 */

static int resS(int Ns, realtype t,
                N_Vector yy, N_Vector yp, N_Vector resval,
                N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                void *user_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype sd1, sd2;
  realtype rs1, rs2, rs3;
  int is;

  data = (UserData) user_data;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);

  for (is=0; is<NS; is++) {

    s1 = Ith(yyS[is],1);
    s2 = Ith(yyS[is],2);
    s3 = Ith(yyS[is],3);

    sd1 = Ith(ypS[is],1);
    sd2 = Ith(ypS[is],2);

    rs1 = sd1 + p1*s1 - p2*y3*s2 - p2*y2*s3;
    rs2 = sd2 - p1*s1 + p2*y3*s2 + p2*y2*s3 + 2*p3*y2*s2;
    rs3 = s1 + s2 + s3;

    switch (is) {
    case 0:
      rs1 += y1;
      rs2 -= y1;
      break;
    case 1:
      rs1 -= y2*y3;
      rs2 += y2*y3;
      break;
    case 2:
      rs2 += y2*y2;
      break;
    }

    Ith(resvalS[is],1) = rs1;
    Ith(resvalS[is],2) = rs2;
    Ith(resvalS[is],3) = rs3;

  }

  return(0);
}

static int rhsQ(realtype t, N_Vector y, N_Vector yp,
              N_Vector ypQ, void* user_data)
{
  UserData data;

  data = (UserData) user_data;

  Ith(ypQ,1) = Ith(y,3);

  Ith(ypQ,2) = data->coef*( Ith(y,1)*Ith(y,1)+
                            Ith(y,2)*Ith(y,2)+
                            Ith(y,3)*Ith(y,3) );

  return(0);
}


/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to idasfwddenx.
 */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = SUNFALSE;
  *sensi_meth = -1;
  *err_con = SUNFALSE;

  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = SUNFALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = SUNTRUE;
  else
    WrongArgs(argv[0]);

  if (*sensi) {

    if (argc != 4)
      WrongArgs(argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = IDA_SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = IDA_STAGGERED;
    else
      WrongArgs(argv[0]);

    if (strcmp(argv[3],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs(argv[0]);
  }

}

static void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim or stg\n");
    printf("         err_con    = t or f\n");

    exit(0);
}


static void PrintIC(N_Vector y, N_Vector yp)
{
  realtype* data;

  data = N_VGetArrayPointer(y);
  printf("\n\nConsistent IC:\n");
  printf("\ty = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", data[0], data[1], data[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#endif

  data = N_VGetArrayPointer(yp);
  printf("\typ= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", data[0], data[1], data[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#endif

}
static void PrintSensIC(N_Vector y, N_Vector yp, N_Vector* yS, N_Vector* ypS)
{
  realtype *sdata;

  sdata = N_VGetArrayPointer(yS[0]);
  printf("                  Sensitivity 1  ");

  printf("\n\ts1 = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  sdata = N_VGetArrayPointer(ypS[0]);
  printf("\ts1'= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif


  printf("                  Sensitivity 2  ");
  sdata = N_VGetArrayPointer(yS[1]);
  printf("\n\ts2 = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  sdata = N_VGetArrayPointer(ypS[1]);
  printf("\ts2'= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif


  printf("                  Sensitivity 3  ");
  sdata = N_VGetArrayPointer(yS[2]);
  printf("\n\ts3 = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  sdata = N_VGetArrayPointer(ypS[2]);
  printf("\ts3'= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif


}

/*
 * Print current t, step count, order, stepsize, and solution.
 */

static void PrintOutput(void *ida_mem, realtype t, N_Vector u)
{
  long int nst;
  int qu, retval;
  realtype hu, *udata;

  udata = N_VGetArrayPointer(u);

  retval = IDAGetNumSteps(ida_mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastOrder(ida_mem, &qu);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetLastStep(ida_mem, &hu);
  check_retval(&retval, "IDAGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#else
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#endif

  printf("                  Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", udata[0], udata[1], udata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#endif

}

/*
 * Print sensitivities.
*/

static void PrintSensOutput(N_Vector *uS)
{
  realtype *sdata;

  sdata = N_VGetArrayPointer(uS[0]);
  printf("                  Sensitivity 1  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

  sdata = N_VGetArrayPointer(uS[1]);
  printf("                  Sensitivity 2  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

  sdata = N_VGetArrayPointer(uS[2]);
  printf("                  Sensitivity 3  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if retval < 0 */
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1);
    }
  } else if (opt == 2 && returnvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  }

  return(0);
}
