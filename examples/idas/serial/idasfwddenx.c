/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2007-06-05 21:03:55 $
 * -----------------------------------------------------------------
 * Programmer(s): Cosmin Petra and Radu Serban @ LLNL
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
 *(error control set on TRUE or FALSE, respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % idasfwddenx -nosensi
 * If sensitivities are to be computed:
 *    % idasfwddenx -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <idas/idas.h>           /* prototypes for IDAS fcts. and consts. */
#include <idas/idas_dense.h>
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

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
} *UserData;

/* Prototypes of functions by IDAS */

static int res(realtype t, N_Vector y, N_Vector yp, N_Vector resval, void *user_data);

static int resS(int Ns, realtype t, 
                N_Vector y, N_Vector yp, N_Vector resval,
                N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                void *user_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth, 
                        booleantype *err_con);
static void WrongArgs(char *name);
static void PrintOutput(void *ida_mem, realtype t, N_Vector u);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(void *ida_mem, booleantype sensi);
static int check_flag(void *flagvalue, char *funcname, int opt);
/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *ida_mem;
  UserData data;
  realtype reltol, t, tout;
  N_Vector y, yp, abstol;
  int iout, flag;

  realtype pbar[NS];
  int is; 
  N_Vector *yS, *ypS;
  booleantype sensi, err_con;
  int sensi_meth;

  ida_mem = NULL;
  data    = NULL;
  y       =  NULL;
  yS      = NULL;
  ypS     = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2)) return(1);
  data->p[0] = RCONST(0.04);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);

  /* Initial conditions */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  Ith(y,1) = ONE;
  Ith(y,2) = ZERO;
  Ith(y,3) = ZERO;

  yp = N_VNew_Serial(NEQ);
  if(check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);

  Ith(yp,1) = RCONST(-0.04);
  Ith(yp,2) = RCONST(0.04);
  Ith(yp,3) = ZERO;  

  /* Create IDAS object */
  ida_mem = IDACreate();
  if (check_flag((void *)ida_mem, "IDACreate", 0)) return(1);

  /* Allocate space for IDAS */
  flag = IDAInit(ida_mem, res, T0, y, yp);
  if (check_flag(&flag, "IDAInit", 1)) return(1);

  /* Specify scalar relative tol. and vector absolute tol. */
  reltol = RCONST(1.0e-4);
  abstol = N_VNew_Serial(NEQ);
  Ith(abstol,1) = RCONST(1.0e-8);
  Ith(abstol,2) = RCONST(1.0e-14);
  Ith(abstol,3) = RCONST(1.0e-6);
  flag = IDASVtolerances(ida_mem, reltol, abstol);
  if (check_flag(&flag, "IDASVtolerances", 1)) return(1);

  /* Attach user data */
  flag = IDASetUserData(ida_mem, data);
  if (check_flag(&flag, "IDASetUserData", 1)) return(1);

  /* Attach linear solver */
  flag = IDADense(ida_mem, NEQ);
  if (check_flag(&flag, "IDADense", 1)) return(1);

  printf("\n3-species chemical kinetics problem\n");

  /* Sensitivity-related settings */
  if (sensi) {

    pbar[0] = data->p[0];
    pbar[1] = data->p[1];
    pbar[2] = data->p[2];

    yS = N_VCloneVectorArray_Serial(NS, y);
    if (check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0)) return(1);
    for (is=0;is<NS;is++) N_VConst(ZERO, yS[is]);
    
    ypS = N_VCloneVectorArray_Serial(NS, y);
    if (check_flag((void *)ypS, "N_VCloneVectorArray_Serial", 0)) return(1);
    for (is=0;is<NS;is++) N_VConst(ZERO, ypS[is]);

    /* Only non-zero sensitivity I.C. are ypS[0] */
    Ith(ypS[0],1) = -ONE;
    Ith(ypS[0],2) =  ONE;

    flag = IDASensInit(ida_mem, NS, sensi_meth, resS, yS, ypS);
    if(check_flag(&flag, "IDASensInit", 1)) return(1);

    flag = IDASensEEtolerances(ida_mem);
    if(check_flag(&flag, "IDASensEEtolerances", 1)) return(1);

    flag = IDASetSensErrCon(ida_mem, err_con);
    if (check_flag(&flag, "IDASetSensErrCon", 1)) return(1);

    flag = IDASetSensParams(ida_mem, data->p, pbar, NULL);
    if (check_flag(&flag, "IDASetSensParams", 1)) return(1);

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
  
  /* In loop over output points, call IDA, print results, test for error */
  
  printf("\n\n");
  printf("===========================================");
  printf("============================\n");
  printf("     T     Q       H      NST           y1");
  printf("           y2           y3    \n");
  printf("===========================================");
  printf("============================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {

    flag = IDASolve(ida_mem, tout, &t, y, yp, IDA_NORMAL);
    if (check_flag(&flag, "IDASolve", 1)) break;

    PrintOutput(ida_mem, t, y);

    if (sensi) {
      flag = IDAGetSens(ida_mem, &t, yS);
      if (check_flag(&flag, "IDAGetSens", 1)) break;
      PrintOutputS(yS);
    } 
    printf("-----------------------------------------");
    printf("------------------------------\n");

  }

  /* Print final statistics */
  PrintFinalStats(ida_mem, sensi);

  /* Free memory */
  N_VDestroy_Serial(y);
  if (sensi) {
    N_VDestroyVectorArray_Serial(yS, NS);
  }
  free(data);
  IDAFree(&ida_mem);

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
  realtype yp1, yp2, yp3;

  data = (UserData) user_data;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);

  yp1 = Ith(yp,1);
  yp2 = Ith(yp,2);
  yp3 = Ith(yp,3);

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
  realtype yp1, yp2, yp3;
  realtype s1, s2, s3;
  realtype sd1, sd2, sd3;
  realtype rs1, rs2, rs3;
  int is;

  data = (UserData) user_data;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);

  yp1 = Ith(yp,1);
  yp2 = Ith(yp,2);
  yp3 = Ith(yp,3);

  for (is=0; is<NS; is++) {

    s1 = Ith(yyS[is],1);
    s2 = Ith(yyS[is],2);
    s3 = Ith(yyS[is],3);

    sd1 = Ith(ypS[is],1);
    sd2 = Ith(ypS[is],2);
    sd3 = Ith(ypS[is],3);

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
  *sensi = FALSE;
  *sensi_meth = -1;
  *err_con = FALSE;

  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = FALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = TRUE;
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
      *err_con = TRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = FALSE;
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

/*
 * Print current t, step count, order, stepsize, and solution.
 */

static void PrintOutput(void *ida_mem, realtype t, N_Vector u)
{
  long int nst;
  int qu, flag;
  realtype hu, *udata;
  
  udata = NV_DATA_S(u);

  flag = IDAGetNumSteps(ida_mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = IDAGetLastOrder(ida_mem, &qu);
  check_flag(&flag, "IDAGetLastOrder", 1);
  flag = IDAGetLastStep(ida_mem, &hu);
  check_flag(&flag, "IDAGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.3le %2d  %8.3le %5ld\n", t, qu, hu, nst);
#else
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#endif

  printf("                  Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", udata[0], udata[1], udata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4le %12.4le %12.4le \n", udata[0], udata[1], udata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#endif

}

/* 
 * Print sensitivities.
*/

static void PrintOutputS(N_Vector *uS)
{
  realtype *sdata;

  sdata = NV_DATA_S(uS[0]);
  printf("                  Sensitivity 1  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4le %12.4le %12.4le \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  
  sdata = NV_DATA_S(uS[1]);
  printf("                  Sensitivity 2  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4le %12.4le %12.4le \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

  sdata = NV_DATA_S(uS[2]);
  printf("                  Sensitivity 3  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4le %12.4le %12.4le \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
}

/* 
 * Print some final statistics from the IDAS memory.
 */

static void PrintFinalStats(void *ida_mem, booleantype sensi)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nje, nfeLS;
  int flag;

  flag = IDAGetNumSteps(ida_mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = IDAGetNumResEvals(ida_mem, &nfe);
  check_flag(&flag, "IDAGetNumRhsEvals", 1);
  flag = IDAGetNumLinSolvSetups(ida_mem, &nsetups);
  check_flag(&flag, "IDAGetNumLinSolvSetups", 1);
  flag = IDAGetNumErrTestFails(ida_mem, &netf);
  check_flag(&flag, "IDAGetNumErrTestFails", 1);
  flag = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  check_flag(&flag, "IDAGetNumNonlinSolvIters", 1);
  flag = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  check_flag(&flag, "IDAGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    flag = IDAGetNumSensResEvals(ida_mem, &nfSe);
    check_flag(&flag, "IDAGetNumSensRhsEvals", 1);
    flag = IDAGetNumResEvalsSens(ida_mem, &nfeS);
    check_flag(&flag, "IDAGetNumRhsEvalsSens", 1);
    flag = IDAGetNumSensLinSolvSetups(ida_mem, &nsetupsS);
    check_flag(&flag, "IDAGetNumSensLinSolvSetups", 1);
    flag = IDAGetNumSensErrTestFails(ida_mem, &netfS);
    check_flag(&flag, "IDAGetNumSensErrTestFails", 1);
    flag = IDAGetNumSensNonlinSolvIters(ida_mem, &nniS);
    check_flag(&flag, "IDAGetNumSensNonlinSolvIters", 1);
    flag = IDAGetNumSensNonlinSolvConvFails(ida_mem, &ncfnS);
    check_flag(&flag, "IDAGetNumSensNonlinSolvConvFails", 1);
  }

  flag = IDADlsGetNumJacEvals(ida_mem, &nje);
  check_flag(&flag, "IDADlsGetNumJacEvals", 1);
  flag = IDADlsGetNumResEvals(ida_mem, &nfeLS);
  check_flag(&flag, "IDADlsGetNumResEvals", 1);

  printf("\nFinal Statistics\n\n");
  printf("nst     = %5ld\n\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

  printf("\n");
  printf("nje    = %5ld    nfeLS     = %5ld\n", nje, nfeLS);

}

/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
