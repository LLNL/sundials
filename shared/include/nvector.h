/*******************************************************************
 *                                                                 *
 * File          : nvector.h                                       *
 * Programmers   : Radu Serban, LLNL                               *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the header file for a generic NVECTOR package.          *
 * It defines the N_Vector and M_Env structures:                   *
 *   M_Env has an implementation-dependent 'content' field         *
 *         which contains the data needed to generate a new        *
 *         nvector in that implementation and an 'ops' filed       *
 *         which is a structure listing operations acting on       *
 *         such nvectors.                                          *
 *   N_Vector has an implementation-dependent 'content' field      *
 *         which contains the description and actual data of       *
 *         the nvector and a 'menv' field which points to the      *
 *         M_Env structure used in creating the nvector.           * 
 *                                                                 *
 * Part I of this file contains type declarations for the          *
 * the following structures: _generic_M_Env, _generic_N_Vector,    *
 * and _generic_N_Vector_Ops, as well as references to pointers    *
 * to such structures (M_Env and N_Vector).                        *
 *                                                                 *
 * Part II of this file contains the prototypes for the vector     *
 * kernels which operate on N_Vector.                              * 
 *                                                                 *
 * A particular implementation of an NVECTOR package must then     *
 * specify the 'content' fields of M_Env and N_Vector, define      *
 * the propotypes for kernel operations on those N_Vectors         *
 * (NOTE: kernel routine names must be unique to that              *
 * implementation), and finally provide an initialization          *
 * routine (which generates an M_Env with that particular          *
 * 'content' field and links the defined vector kernel routines    *
 * into the 'ops' field).                                          *
 *                                                                 *
 *******************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif
  
#ifndef included_nvector_h
#define included_nvector_h
  
#include "sundialstypes.h"  /* definition of types */
  
/****************************************************************
 * Generic definitions of machine environment and N_Vector      *
 ****************************************************************/

/* Forward reference for pointer to N_Vector_Ops object */
typedef struct _generic_N_Vector_Ops *N_Vector_Ops;

/* Forward reference for pointer to M_Env object */
typedef struct _generic_M_Env *M_Env;
  
/* Forward reference for pointer to N_Vector object */
typedef struct _generic_N_Vector *N_Vector;

/* Define array of N_Vectors */
typedef N_Vector *N_Vector_S;

/* Structure containing function pointers to vector operations  */  
struct _generic_N_Vector_Ops {
  N_Vector    (*nvnew)(integertype, M_Env);
  N_Vector_S  (*nvnewS)(integertype, integertype, M_Env);
  void        (*nvfree)(N_Vector);
  void        (*nvfreeS)(integertype, N_Vector_S);
  N_Vector    (*nvmake)(integertype, realtype *, M_Env);
  void        (*nvdispose)(N_Vector);
  realtype*   (*nvgetdata)(N_Vector);
  void        (*nvsetdata)(realtype *, N_Vector);
  void        (*nvlinearsum)(realtype, N_Vector, realtype, N_Vector, N_Vector); 
  void        (*nvconst)(realtype, N_Vector);
  void        (*nvprod)(N_Vector, N_Vector, N_Vector);
  void        (*nvdiv)(N_Vector, N_Vector, N_Vector);
  void        (*nvscale)(realtype, N_Vector, N_Vector);
  void        (*nvabs)(N_Vector, N_Vector);
  void        (*nvinv)(N_Vector, N_Vector);
  void        (*nvaddconst)(N_Vector, realtype, N_Vector);
  realtype    (*nvdotprod)(N_Vector, N_Vector);
  realtype    (*nvmaxnorm)(N_Vector);
  realtype    (*nvwrmsnorm)(N_Vector, N_Vector);
  realtype    (*nvmin)(N_Vector);
  realtype    (*nvwl2norm)(N_Vector, N_Vector);
  realtype    (*nvl1norm)(N_Vector);
  void        (*nvonemask)(N_Vector);
  void        (*nvcompare)(realtype, N_Vector, N_Vector);
  booleantype (*nvinvtest)(N_Vector, N_Vector);
  booleantype (*nvconstrprodpos)(N_Vector, N_Vector);
  booleantype (*nvconstrmask)(N_Vector, N_Vector, N_Vector);
  realtype    (*nvminquotient)(N_Vector, N_Vector);
  void        (*nvprint)(N_Vector);
};
  
/* A machine environment is a structure with an implementation
   dependent 'content' representation (used to generate a new vector
   in that implementation), a set of operations defined in the above 
   structure, and an ID tag */
struct _generic_M_Env {
  void *content;
  struct _generic_N_Vector_Ops *ops;
  char tag[8];
};

/* A vector is a structure with an implementation dependent content
   representation and a pointer to the machine environment 
   corresponding to that implementation */
struct _generic_N_Vector {
  void *content;
  struct _generic_M_Env *menv;
};
  
/****************************************************************
 * Functions exported by nvector                                *
 ****************************************************************/

/*--------------------------------------------------------------*
 * Function : N_VNew                                            *
 * Usage    : v = N_VNew(n, machEnv);                           *
 *--------------------------------------------------------------*
 * Returns a new N_Vector of length n. The parameter machEnv    *
 * is a pointer to machine environment-specific information.    *
 * If there is not enough memory for a new N_Vector, then       *
 * N_VNew returns NULL.                                         *
 *--------------------------------------------------------------*/
  
N_Vector N_VNew(integertype n, M_Env machEnv);

/*--------------------------------------------------------------*
 * Function : N_VNew_S                                          *
 * Usage    : v = N_VNew_S(ns, n, machEnv);                     *
 *--------------------------------------------------------------*
 * Returns an array of ns new N_Vectors of length n. The        *
 * parameter machEnv is a pointer to machine environment        *
 * specific information.                                        *
 * If there is not enough memory for a new array of N_Vectors   *
 * or for one of the components, then N_VNew_S returns NULL.    *
 *--------------------------------------------------------------*/

N_Vector_S N_VNew_S(integertype ns, integertype n, M_Env machEnv);

/*--------------------------------------------------------------*
 * Function : N_VFree                                           *
 * Usage    : N_VFree(v);                                       *
 *--------------------------------------------------------------*
 * Frees the N_Vector v. It is illegal to use v after the call  *
 * N_VFree(v).                                                  *
 *--------------------------------------------------------------*/

void N_VFree(N_Vector v);

/*--------------------------------------------------------------*
 * Function : N_VFree_S                                         *
 * Usage    : N_VFree_S(ns, vs);                                *
 *--------------------------------------------------------------*
 * Frees the array of ns N_Vectors vs.                          *
 * It is illegal to use vs after the call N_VFree_S(Ns,vs).     *
 *--------------------------------------------------------------*/

void N_VFree_S(integertype ns, N_Vector_S vs);

/*--------------------------------------------------------------*
 * Function : N_VMake                                           *
 * Usage    : v = N_VMake(n, v_data, machEnv);                  *
 *--------------------------------------------------------------*
 * Creates an N_Vector with component array data allocated by   *
 * the user.                                                    *
 *--------------------------------------------------------------*/

N_Vector N_VMake(integertype n, realtype *v_data, M_Env machEnv);

/*--------------------------------------------------------------*
 * Function : N_VDispose                                        *
 * Usage    : N_VDispose(v);                                    *
 *--------------------------------------------------------------*
 * Destroys an N_Vector with component array data allocated by  *
 * the user.                                                    *
 *--------------------------------------------------------------*/

void N_VDispose(N_Vector v);

/*--------------------------------------------------------------*
 * Function : N_VGetData                                        *
 * Usage    : v_data = N_VGetData(v);                           *
 *--------------------------------------------------------------*
 * Extracts the data component array from the N_Vector v.       *
 * Note: this routine is used in the solver-specific interfaces *
 *       to the dense and banded linear solvers, as well as the *
 *       interfaces to the banded preconditioners provided with *
 *       SUNDIALS. It needs not be implemented by a user        *
 *       defined NVECTOR module, if these linear solvers are not*
 *       used.                                                  *
 *--------------------------------------------------------------*/

realtype *N_VGetData(N_Vector v);

/*--------------------------------------------------------------*
 * Function : N_VSetData                                        *
 * Usage    : N_VSetData(v_data, v);                            *
 *--------------------------------------------------------------*
 * Attaches the data component array v_data to the N_Vector v.  *
 * Note: this routine is used in the solver-specific interfaces *
 *       to the dense and banded linear solvers, as well as the *
 *       interfaces to the banded preconditioners provided with *
 *       SUNDIALS. It needs not be implemented by a user        *
 *       defined NVECTOR module, if these linear solvers are not*
 *       used.                                                  *
 *--------------------------------------------------------------*/

void N_VSetData(realtype *v_data, N_Vector v);

/*--------------------------------------------------------------*
 * Function  : N_VLinearSum                                     *
 * Operation : z = a x + b y                                    *
 *--------------------------------------------------------------*/

void N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, 
                  N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VConst                                         *
 * Operation : z[i] = c for i=0, 1, ..., N-1                    *
 *--------------------------------------------------------------*/

void N_VConst(realtype c, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VProd                                          *
 * Operation : z[i] = x[i] * y[i] for i=0, 1, ..., N-1          *
 *--------------------------------------------------------------*/

void N_VProd(N_Vector x, N_Vector y, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VDiv                                           *
 * Operation : z[i] = x[i] / y[i] for i=0, 1, ..., N-1          *
 *--------------------------------------------------------------*/

void N_VDiv(N_Vector x, N_Vector y, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VScale                                         *
 * Operation : z = c x                                          *
 *--------------------------------------------------------------*/

void N_VScale(realtype c, N_Vector x, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VAbs                                           *
 * Operation : z[i] = |x[i]|,   for i=0, 1, ..., N-1            *
 *--------------------------------------------------------------*/

void N_VAbs(N_Vector x, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VInv                                           *
 * Operation : z[i] = 1.0 / x[i] for i = 0, 1, ..., N-1         *
 *--------------------------------------------------------------*
 * This routine does not check for division by 0. It should be  *
 * called only with an N_Vector x which is guaranteed to have   *
 * all non-zero components.                                     *
 *--------------------------------------------------------------*/

void N_VInv(N_Vector x, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VAddConst                                      *
 * Operation : z[i] = x[i] + b   for i = 0, 1, ..., N-1         *
 *--------------------------------------------------------------*/

void N_VAddConst(N_Vector x, realtype b, N_Vector z);

/*--------------------------------------------------------------*
 * Function : N_VDotProd                                        *
 * Usage    : dotprod = N_VDotProd(x, y);                       *
 *--------------------------------------------------------------*
 * Returns the value of the ordinary dot product of x and y:    *
 * -> sum (i=0 to N-1) {x[i] * y[i]}                            *
 * Returns 0.0 if N <= 0.                                       *
 *--------------------------------------------------------------*/

realtype N_VDotProd(N_Vector x, N_Vector y);

/*--------------------------------------------------------------*
 * Function : N_VMaxNorm                                        *
 * Usage    : maxnorm = N_VMaxNorm(x);                          *
 *--------------------------------------------------------------*
 * Returns the maximum norm of x:                               *
 * -> max (i=0 to N-1) |x[i]|                                   *
 * Returns 0.0 if N <= 0.                                       *
 *--------------------------------------------------------------*/

realtype N_VMaxNorm(N_Vector x);

/*--------------------------------------------------------------*
 * Function : N_VWrmsNorm                                       *
 * Usage    : wrmsnorm = N_VWrmsNorm(x, w);                     *
 *--------------------------------------------------------------*
 * Returns the weighted root mean square norm of x with         *
 * weight vector w:                                             *
 * -> sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) / N]           *
 * Returns 0.0 if N <= 0.                                       *
 *--------------------------------------------------------------*/

realtype N_VWrmsNorm(N_Vector x, N_Vector w);

/*--------------------------------------------------------------*
 * Function : N_VMin                                            *
 * Usage    : min = N_VMin(x);                                  *
 *--------------------------------------------------------------*
 * Returns the smallest element of x:                           *
 * -> min (i=0 to N-1) x[i]                                     *
 * Returns 0.0 if N <= 0.                                       *
 *--------------------------------------------------------------*/

realtype N_VMin(N_Vector x);

/*--------------------------------------------------------------*
 * Function : N_VWL2Norm                                        *
 * Usage    : wl2norm = N_VWL2Norm(x, w);                       *
 *--------------------------------------------------------------*
 * Returns the weighted Euclidean L2 norm of x with             *
 * weight vector w:                                             *
 * -> sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) ]              *
 * Returns 0.0 if N <= 0.                                       *
 *--------------------------------------------------------------*/

realtype N_VWL2Norm(N_Vector x, N_Vector w);

/*--------------------------------------------------------------*
 * Function : N_VL1Norm                                         *
 * Usage    : l1norm = N_VL1Norm(x);                            *
 *--------------------------------------------------------------*
 * Returns the L1 norm of x:                                    *
 * -> sum (i=0 to N-1) {ABS(x[i])}                              *
 * Returns 0.0 if N <= 0.                                       *
 *--------------------------------------------------------------*/

realtype N_VL1Norm(N_Vector x);

/*--------------------------------------------------------------*
 * Function  : N_VOneMask                                       *
 * Operation : x[i] = 1.0 if |x[i]| != 0.  i = 0, 1, ..., N-1   *
 *                    0.0 otherwise                             *
 *--------------------------------------------------------------*/

void N_VOneMask(N_Vector x);

/*--------------------------------------------------------------*
 * Function  : N_VCompare                                       *
 * Operation : z[i] = 1.0 if |x[i]| >= c   i = 0, 1, ..., N-1   *
 *                    0.0 otherwise                             *
 *--------------------------------------------------------------*/

void N_VCompare(realtype c, N_Vector x, N_Vector z);

/*--------------------------------------------------------------*
 * Function  : N_VInvTest                                       *
 * Operation : z[i] = 1.0 / x[i] with a test for x[i]==0.0      *
 *             before inverting x[i].                           *
 *--------------------------------------------------------------*
 * This routine returns TRUE if all components of x are         *
 * non-zero (successful inversion) and returns FALSE            *
 * otherwise.                                                   *
 *--------------------------------------------------------------*/

booleantype N_VInvTest(N_Vector x, N_Vector z);

/*--------------------------------------------------------------*
 * Function : N_VConstrProdPos                                  *
 * Usage    : booltest = N_VConstrProdPos(c,x);                 *
 *--------------------------------------------------------------*
 * Returns a boolean equal to                                   *
 *   FALSE if some c[i] != 0.0 and x[i]*c[i] <= 0.0,  or        *
 *   TRUE otherwise.                                            *
 *                                                              *
 * This routine is used for constraint checking.                *
 *--------------------------------------------------------------*/

booleantype N_VConstrProdPos(N_Vector c, N_Vector x);

/*--------------------------------------------------------------*
 * Function  : N_VConstrMask                                    *
 * Operation : m[i] = 1.0 if constraint test fails for x[i]     *
 *             m[i] = 0.0 if constraint test passes for x[i]    *
 * where the constraint tests are as follows:                   *
 *    If c[i] = 2.0,  then x[i] must be > 0.0.                  *
 *    If c[i] = 1.0,  then x[i] must be >= 0.0.                 *
 *    If c[i] = -1.0, then x[i] must be <= 0.0.                 *
 *    If c[i] = -2.0, then x[i] must be < 0.0.                  *
 *--------------------------------------------------------------*
 * This routine returns a boolean FALSE if any element failed   *
 * the constraint test, TRUE if all passed.  It also sets a     *
 * mask vector m, with elements equal to 1.0 where the          *
 * corresponding constraint test failed, and equal to 0.0       *
 * where the constraint test passed.                            *
 * This routine is specialized in that it is used only for      *
 * constraint checking.                                         *
 *--------------------------------------------------------------*/

booleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m);

/*--------------------------------------------------------------*
 * Function  : N_VMinQuotient                                   *
 * Operation : minq  = min ( num[i]/denom[i]) over all i such   *
 *             that   denom[i] != 0.                            *
 *--------------------------------------------------------------*
 * This routine returns the minimum of the quotients obtained   *
 * by term-wise dividing num[i] by denom[i]. A zero element     *
 * in denom will be skipped. If no such quotients are found,    *
 * then the large value 1.e99 is returned.                      *
 *--------------------------------------------------------------*/

realtype N_VMinQuotient(N_Vector num, N_Vector denom);

/*--------------------------------------------------------------*
 * Function : N_VPrint                                          *
 * Usage    : N_VPrint(x);                                      *
 *--------------------------------------------------------------*
 * Prints the N_Vector x to stdout.                             *
 * This routine is provided as an aid in debugging code which   *
 * uses this vector package.                                    *
 *--------------------------------------------------------------*/

void N_VPrint(N_Vector x);


#endif

#ifdef __cplusplus
}
#endif
