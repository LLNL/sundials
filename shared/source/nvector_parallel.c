/****************************************************************
 *                                                              *
 * File          : nvector.c                                    *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh,           *
 *                 Radu Serban, and Allan G. Taylor, LLNL       *
 * Version of    : 30 October 2001                              *
 *--------------------------------------------------------------*
 *                                                              *
 * This is the implementation file for an MPI NVECTOR           *
 * package. It contains the implementation of the N_Vector      *
 * kernels listed in nvector.h.                                 *
 * It uses MPI for message-passing.                             *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "nvector.h"
#include "llnltyps.h"
#include "llnlmath.h" 


#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error message */

#define BAD_N1  "PVecInitMPI-- Sum of local vector lengths differs from "
#define BAD_N2  "input global length. \n\n"
#define BAD_N    BAD_N1 BAD_N2

/* Private Helper Prototypes */

static real PVecAllReduce(real d, int op, machEnvType machEnv);
/* Reduction operations add/max/min over the processor group */
static void VCopy(N_Vector x, N_Vector z); /* z=x */
static void VSum(N_Vector x, N_Vector y, N_Vector z); /* z=x+y */
static void VDiff(N_Vector x, N_Vector y, N_Vector z); /* z=x-y */
static void VNeg(N_Vector x, N_Vector z); /* z=-x */
/* z=c(x+y) */
static void VScaleSum(real c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff(real c, N_Vector x, N_Vector y, N_Vector z); 
static void VLin1(real a, N_Vector x, N_Vector y, N_Vector z); /* z=ax+y */
static void VLin2(real a, N_Vector x, N_Vector y, N_Vector z); /* z=ax-y */
static void Vaxpy(real a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy(real a, N_Vector x); /* x <- ax */


/********************* Exported Functions ************************/


void *PVecInitMPI(MPI_Comm comm,  integer local_vec_length, 
		  integer global_vec_length, int *argc, char ***argv)
{
  int initflag, initerr;
  integer n, Nsum;
  machEnvType env;

  /* Create structure env and begin loading it */
  env = (machEnvType) malloc (sizeof *env);
  if (env == NULL) return(NULL);

  env->local_vec_length = local_vec_length;
  env->global_vec_length = global_vec_length;

  MPI_Initialized(&initflag);
  if (!initflag) {
    initerr = MPI_Init(argc,argv);
    if (initerr != MPI_SUCCESS) return(NULL);
    }
  env->init_by_user = initflag;

  env->comm = comm;

  /* If this PE is inactive, return now */
  if (local_vec_length <= 0) return((void *)env);

  /* Compute global length as sum of local lengths */
  n = local_vec_length;
  MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  env->global_vec_length = Nsum;

  /* Check input global length against computed value */
  if (Nsum != global_vec_length) {
    printf(BAD_N);
    PVecFreeMPI(env);
    return(NULL);
    } 

  /* Return the pointer env */
  return((void *)env);
}


void PVecFreeMPI(void *machEnv)
{
  machEnvType env;

  env = (machEnvType) machEnv;
  if (env == NULL) return;

  if (!(env->init_by_user)) MPI_Finalize();

  free(machEnv);
}

 
N_Vector N_VNew(integer N, machEnvType machEnv)
{
  N_Vector v;
  int N_local, N_global;

  if (N <= 0) return(NULL);
  if (machEnv == NULL) return(NULL);

  N_local = machEnv->local_vec_length;
  N_global = machEnv->global_vec_length;

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  v->data = (real *) malloc(N_local * sizeof(real));
  if (v->data == NULL) {
    free(v);
    return(NULL);
  }

  v->length = N_local;
  v->global_length = N_global;
  v->machEnv = machEnv;
  
  return(v);
}


N_Vector *N_VNew_S(integer Ns, integer N, machEnvType machEnv)
{
    N_Vector *vs;
    integer is, j;

    if (Ns <= 0 || N <= 0) return(NULL);
    if (machEnv == NULL) return(NULL);

    vs = (N_Vector *) malloc(Ns * sizeof(N_Vector *));
    if (vs == NULL) return(NULL);

    for (is=0; is<Ns; is++) {
        vs[is] = N_VNew(N,machEnv);
        if (vs[is] == NULL) {
            for (j=0; j<is; j++) N_VFree(vs[j]);
            free(vs);
            return(NULL);
        }
    }
    
    return(vs);
}


void N_VFree(N_Vector x)
{
  free(x->data);
  free(x);
}


void N_VFree_S(integer Ns, N_Vector *vs)
{
    integer is;
    
    for (is=0; is<Ns; is++) N_VFree(vs[is]);
    free(vs);
}


void N_VLinearSum(real a, N_Vector x, real b, N_Vector y, N_Vector z)
{
  integer i, N;
  real c, *xd, *yd, *zd;
  N_Vector v1, v2;
  boole test;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++) 
    *zd++ = a * (*xd++) + b * (*yd++);
}


void N_VConst(real c, N_Vector z)
{
  integer i, N;
  real *zd;

  N = z->length;
  zd = z->data;

  for (i=0; i < N; i++) 
    *zd++ = c;
}


void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = (*xd++) * (*yd++);
}


void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = (*xd++) / (*yd++);
}


void N_VScale(real c, N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd;

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VScaleBy(c, x);
    return;
  }

  if (c == ONE) {
    VCopy(x, z);
  } else if (c == -ONE) {
    VNeg(x, z);
  } else {
    N = x->length;
    xd = x->data;
    zd = z->data;
    for (i=0; i < N; i++) *zd++ = c * (*xd++);
  }
}


void N_VAbs(N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; i++, xd++, zd++)
    *zd = ABS(*xd);
}


void N_VInv(N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = ONE / (*xd++);
}


void N_VAddConst(N_Vector x, real b, N_Vector z)
{
  integer i, N;
  real *xd, *zd;
  
  N = x->length;
  xd = x->data;
  zd = z->data;
  
  for (i=0; i < N; i++) *zd++ = (*xd++) + b; 
}


real N_VDotProd(N_Vector x, N_Vector y)
{
  integer i, loclen;
  real sum = ZERO, *xd, *yd, gsum;
  machEnvType machenv;

  loclen = x->length;
  xd = x->data;
  yd = y->data;
  machenv = x->machEnv;

  for (i=0; i < loclen; i++) sum += xd[i] * yd[i];

  gsum = PVecAllReduce(sum, 1, machenv);
  return(gsum);
}


real N_VMaxNorm(N_Vector x)
{
  integer i, N;
  real max = ZERO, *xd, gmax;
  machEnvType machenv;

  N = x->length;
  xd = x->data;
  machenv = x->machEnv;

  for (i=0; i < N; i++, xd++) {
    if (ABS(*xd) > max) max = ABS(*xd);
  }
   
  gmax = PVecAllReduce(max, 2, machenv);
  return(gmax);
}


real N_VWrmsNorm(N_Vector x, N_Vector w)
{
  integer i, N, N_global;
  real sum = ZERO, prodi, *xd, *wd, gsum;
  machEnvType machenv;

  N = x->length;
  N_global = x->global_length;
  xd = x->data;
  wd = w->data;
  machenv = x->machEnv;

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  gsum = PVecAllReduce(sum, 1, machenv);
  return(RSqrt(gsum / N_global));
}

real N_VMin(N_Vector x)
{
  integer i, N;
  real min, *xd, gmin;
  machEnvType machenv;

  N = x->length;
  xd = x->data;
  min = xd[0];
  machenv = x->machEnv;

  xd++;
  for (i=1; i < N; i++, xd++) {
    if ((*xd) < min) min = *xd;
  }

  gmin = PVecAllReduce(min, 3, machenv);
  return(gmin);
}


real N_VWL2Norm(N_Vector x, N_Vector w)
{
  integer i, N;
  real sum = ZERO, prodi, *xd, *wd, gsum;
  machEnvType machenv;

  N = x->length;
  xd = x->data;
  wd = w->data;
  machenv = x->machEnv;

  for (i=0; i < N; i++) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  gsum = PVecAllReduce(sum, 1, machenv);
  return(RSqrt(gsum));
}


real N_VL1Norm(N_Vector x)
{
  integer i, N;
  real sum = ZERO, gsum, *xd;
  machEnvType machenv;


  N = x->length;
  xd = x->data;
  machenv = x->machEnv;

  for (i=0; i<N; i++,xd++) sum += ABS(*xd);

  gsum = PVecAllReduce(sum, 1, machenv);

  return(gsum);
}


void N_VOneMask(N_Vector x)
{
  integer i, N;
  real *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i<N; i++,xd++) {
    if (*xd != ZERO) *xd = ONE;
  }
}


void N_VCompare(real c, N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd;
  
  N = x->length;
  xd = x->data;
  zd = z->data;
  
  for (i=0; i < N; i++, xd++, zd++) {
    *zd = (ABS(*xd) >= c) ? ONE : ZERO;
  }
}


boole N_VInvTest(N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd, val, gval;
  machEnvType machenv;

  N = x->length;
  xd = x->data;
  zd = z->data;
  machenv = x->machEnv;

  val = ONE;
  for (i=0; i < N; i++) {
    if (*xd == ZERO) 
      val = ZERO;
    else
      *zd++ = ONE / (*xd++);
  }

  gval = PVecAllReduce(val, 3, machenv);
  if (gval == ZERO)
    return(FALSE);
  else
    return(TRUE);
}


boole N_VConstrProdPos(N_Vector c, N_Vector x)
{
  integer i, N;
  boole test;
  real  *xd, *cd;
  machEnvType machenv;

  N =  x->length;
  xd = x->data;
  cd = c->data;
  machenv = x->machEnv;

  test = TRUE;
  for (i=0; i < N; i++, xd++,cd++) {
    if (*cd != ZERO) {
      if ( (*xd)*(*cd) <= ZERO) {
	test = FALSE;
        break;
      }
    }
  }
  return((boole)PVecAllReduce((real)test, 3, machenv));
}


boole N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
  integer i, N;
  boole test;
  real *cd, *xd, *md;
  machEnvType machenv;
 
  N = x->length;
  cd = c->data;
  xd = x->data;
  md = m->data;
  machenv = x->machEnv;
  
  test = TRUE;

  for (i=0; i<N; i++, cd++, xd++, md++) {
    *md = ZERO;
    if (*cd == ZERO) continue;
    if (*cd > ONEPT5 || (*cd) < -ONEPT5) {
      if ( (*xd)*(*cd) <= ZERO) {test = FALSE; *md = ONE; }
      continue;
    }
    if ( (*cd) > HALF || (*cd) < -HALF) {
      if ( (*xd)*(*cd) < ZERO ) {test = FALSE; *md = ONE; }
    }
  }
  return((boole)PVecAllReduce((real)test, 3, machenv));
}


real N_VMinQuotient(N_Vector num, N_Vector denom)
{
  boole notEvenOnce;
  integer i, N;
  real *nd, *dd, min;
  machEnvType machenv;

  N = num->length;
  nd = num->data;
  dd = denom->data;
  machenv = num->machEnv;
  notEvenOnce = TRUE;

  for (i=0; i<N; i++, nd++, dd++) {
    if (*dd == ZERO) continue;
    else {
      if (notEvenOnce) {
	min = *nd / *dd ;
	notEvenOnce = FALSE;
      }
      else 
	min = MIN(min, (*nd)/(*dd));
    }
  }
  if (notEvenOnce) min = 1.e99;

  return(PVecAllReduce(min, 3, machenv));
}

 
void N_VPrint(N_Vector x)
{
  integer i, N;
  real *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i < N; i++) printf("%g\n", *xd++);

  printf("\n");
}


/***************** Private Helper Functions **********************/

static real PVecAllReduce(real d, int op, machEnvType machEnv)
{
  /* This function does a global reduction.  The operation is
       sum if op = 1,
       max if op = 2,
       min if op = 3.
     The operation is over all processors in the group defined by
     the parameters within machEnv. */

  MPI_Comm comm;
  real out;

  comm = machEnv->comm;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);
}


static void VCopy(N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = *xd++; 
}


static void VSum(N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = (*xd++) + (*yd++);
}


static void VDiff(N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;
 
  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = (*xd++) - (*yd++);
}


static void VNeg(N_Vector x, N_Vector z)
{
  integer i, N;
  real *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = -(*xd++);
}


static void VScaleSum(real c, N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = c * ((*xd++) + (*yd++));
}


static void VScaleDiff(real c, N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = c * ((*xd++) - (*yd++));
}


static void VLin1(real a, N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = a * (*xd++) + (*yd++);
}


static void VLin2(real a, N_Vector x, N_Vector y, N_Vector z)
{
  integer i, N;
  real *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; i++)
    *zd++ = a * (*xd++) - (*yd++);
}

static void Vaxpy(real a, N_Vector x, N_Vector y)
{
  integer i, N;
  real *xd, *yd;

  N = x->length;
  xd = x->data;
  yd = y->data;

  if (a == ONE) {
    for (i=0; i < N; i++)
      *yd++ += (*xd++);
    return;
  }

  if (a == -ONE) {
    for (i=0; i < N; i++)
      *yd++ -= (*xd++);
    return;
  }    

  for (i=0; i < N; i++)
    *yd++ += a * (*xd++);
}

static void VScaleBy(real a, N_Vector x)
{
  integer i, N;
  real *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i < N; i++)
    *xd++ *= a;
}
