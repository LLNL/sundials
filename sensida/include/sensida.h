/******************************************************************
 *                                                                *
 * File          : sensida.h                                      *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 8 December 2000                                *
 *----------------------------------------------------------------*
 * This is the header file for the implementation of computing    *
 * sensitivity information as given in the file:                  *
 * ../source/sensida.c                                            *
 *                                                                *
 ******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _sensida_h
#define _sensida_h

#include "ida.h"
#include "llnltyps.h" 
#include "sensidaspgmr.h"

/******************************************************************
 *                                                                *
 * Function : SensIDAMalloc                                       *
 *----------------------------------------------------------------*
 * SensIDAMalloc initializes a block of sensitivity data.         *
 * It also allocates and initializes memory for a problem to      *
 * to be solved by IDA.                                           *
 *                                                                *
 * Neq     is the total number of differential and algebriac      *
 *         equations in the DAE system.                           *
 *         (In the parallel case, Neq is the global system size,  *
 *         not the local size.)                                   *
 *                                                                *
 * res     is the residual function F in F(t,y,y',p) = 0.         *          
 *                                                                *
 * rdata   is the data memory block (supplied by user) for res.   *
 *         It is passed as a void pointer and is to be cast before*
 *         use in res.                                            *
 *                                                                *
 * t0      is the initial value of t, the independent variable.   *
 *                                                                *
 * y0      is the initial condition vector y(t0).                 *
 *                                                                *
 * yp0     is the initial condition vector y'(t0)                 *
 *                                                                *
 * itol    is the type of tolerances to be used.                  *
 *            The legal values are:                               *
 *               SS (scalar relative and absolute  tolerances),   *
 *               SV (scalar relative tolerance and vector         *
 *                   absolute tolerance).                         *
 *                                                                *
 * rtol    is a pointer to the relative tolerance scalar.         *
 *                                                                *
 * atol    is a pointer (void) to the absolute tolerance scalar   *
 *            or an N_Vector tolerance.                           *
 * (ewt)                                                          *
 *         Both rtol and atol are used to compute the error weight*
 *         vector, ewt. The error test required of a correction   *
 *         delta is that the weighted-RMS norm of delta be less   *
 *         than or equal to 1.0. Other convergence tests use the  *
 *         same norm. The weighting vector used in this norm is   *
 *         ewt. The components of ewt are defined by              *
 *         ewt[i] = 1.0/(rtol*yy[i] + atol[i]). Here, yy is the   *
 *         current approximate solution.  See the routine         *
 *         N_VWrmsNorm for the norm used in this error test.      *
 *                                                                *
 * id      is an N_Vector, required conditionally, which states a *
 *         given element to be either algebraic or differential.  *
 *         A value of 1.0 indicates a differential variable while *
 *         a 0.0 indicates an algebraic variable. 'id' is required*
 *         if optional input SUPPRESSALG is set, or if IDACalcIC  *
 *         is to be called with icopt = CALC_YA_YDP_INIT.         *
 *         Otherwise, 'id' may be NULL.                           *
 *                                                                *
 * constraints  is an N_Vector defining inequality constraints    *
 *         for each component of the solution vector y. If a given*
 *         element of this vector has values +2 or -2, then the   *
 *         corresponding component of y will be constrained to be *
 *         > 0.0 or < 0.0, respectively, while if it is +1 or -1, *
 *         the y component is constrained to be >= 0.0 or <= 0.0, *
 *         respectively. If a component of constraints is 0.0,    *
 *         then no constraint is imposed on the corresponding     *
 *         component of y. The presence of a non-NULL constraints *
 *         vector that is not 0.0 (ZERO) in all components will   *
 *         cause constraint checking to be performed.             *
 *                                                                *
 * errfp   is the file pointer for an error file where all IDA    *
 *            warning and error messages will be written. This    *
 *            parameter can be stdout (standard output), stderr   *
 *            (standard error), a file pointer (corresponding to  *
 *            a user error file opened for writing) returned by   *
 *            fopen, or NULL. If the user passes NULL, then all   *
 *            messages will be written to standard output.        *
 *                                                                *
 * optIn   is a flag (boole) indicating whether there are any     *
 *            optional inputs from the user in the arrays         *
 *            iopt and iopt.                                      *
 *            Pass FALSE to indicate no optional inputs and TRUE  *
 *            to indicate that optional inputs are present.       *
 *                                                                *
 * iopt    is the user-allocated array (of size OPT_SIZE given    *
 *            later) that will hold optional integer inputs and   *
 *            outputs.  The user can pass NULL if he/she does not *
 *            wish to use optional integer inputs or outputs.     *
 *            If optIn is TRUE, the user should preset to 0 those *
 *            locations for which default values are to be used.  *
 *                                                                *
 * ropt    is the user-allocated array (of size OPT_SIZE given    *
 *            later) that will hold optional real inputs and      *
 *            outputs.  The user can pass NULL if he/she does not *
 *            wish to use optional real inputs or outputs.        *
 *            If optIn is TRUE, the user should preset to 0.0 the *
 *            locations for which default values are to be used.  *
 *                                                                *
 * machEnv is a pointer to a block of machine environment-specific*
 *            information.  You may pass NULL for this argument   *
 *            in a serial computer environment.                   *
 *                                                                *
 * p       is a pointer to the array of real parameters given in  *
 *            the DAE system F(t,y,y',p) = 0. This pointer must   *
 *            be stored in the user's rdata data block so that    *
 *            the parameters can be accessed through rdata and    *
 *            this pointer every time res is called.              *
 *                                                                *
 * pbar    is a pointer to an array of nonzero, real values that  *
 *            are used to scale the Ns sensitivity vectors.       *
 *         NOTE: pbar can have the same values as p (if all the   *
 *            parameters are nonzero); however, they cannot be    *
 *            the same array.                                     *
 *                                                                *
 * rhomax  is a real value used for selecting a finite difference *
 *            formula to estimate the residual of the scaled      *
 *            sensitivity equation for the DAE F(t,y,y',p):       *
 *                                                                *
 *            (dF/dy)*w_i + (dF/dy')*w'_i + pbar_i*(dF/dp_i),     *
 *                                                                *
 *            for 1 <= i <= Ns.                                   *
 *                                                                * 
 *          See the following discussion for details concerning   *
 *          the finite difference formulas used for estimating    *
 *          this residual.                                        *
 *----------------------------------------------------------------*
 *** The residual of the scaled sensitivity equations is estimated*
 *** using difference formulas involving the function F(t,y,y',p) *
 *** and forward (or backward) perturbations, such as             *
 *** y + deltay*w_i, y' + deltay*w'_i, and p + deltap*pbar_i.     *
 *** The vectors y, y' and p can be perturbed simultaneously if   *
 *** deltay and deltap have the same value.                       *
 ***                                                              *
 *** The heuristics for selecting deltay and deltap are based on  *
 *** the machine unit roundoff (uround), the user-specified       *
 *** relative error tolerance (rtol), and a weighted root-mean-   *
 *** square norm of the scaled sensitivity vector (wrmsnorm(w_i)):*
 ***                                                              *
 ***   deltap = sqrt(MAX(rtol, uround)),                          *
 ***                                                              *
 ***   deltay = 1/(MAX(wrmsnorm(w_i), 1/deltap)).                 *
 ***                                                              *
 *** Additionally, this subroutine computes                       *
 ***                                                              *
 ***   ratio = MAX(deltay/deltap, deltap/deltay)                  *
 ***                                                              *
 *** to assess if deltay and deltap have comparable magnitudes.   *
 *----------------------------------------------------------------*
 *                                                                *
 * The values and meanings for rhomax are as follows:             *
 *                                                                *
 * rhomax = 0: Choose delta = MIN(deltay, deltap) and use         *
 *             1 centered difference to estimate the residual of  *
 *             the scaled sensitivity equations.                  *
 *             Cost: 2 evaluations of F(t,y,y',p).                *
 *             Accuracy: 2nd order.                               *
 *                                                                *
 * 0 < rhomax < 1: Use centered differences to estimate the sum   *
 *                               (dF/dy)*w_i + (dF/dy')*w'_i,     *
 *                 and the term                                   *
 *                               pbar_i*(dF/dp_i).                *
 *             Cost: 4 evaluations of F(t,y,y',p).                *
 *             Accuracy: 2nd order.                               *
 *                                                                *
 * rhomax >= 1: If (ratio > rhomax), use centered differences to  *
 *              estimate the sum (dF/dy)*w_i + (dF/dy')*w'_i,     *
 *              and the term pbar_i*(dF/dp_i).                    *
 *              Otherwise, choose delta = MIN(deltay, deltap) and *
 *              use 1 centered difference to estimate the residual*
 *              of the scaled sensitivity equations.              *
 *             Cost: 2 or 4 evaluations of F(t,y,y',p).           *
 *             Accuracy: 2nd order.                               *
 *                                                                *
 * -1 < rhomax < 0: Use forward differences to estimate the sum   *
 *                  (dF/dy)*w_i + (dF/dy')*w'_i, and the term     *
 *                  pbar_i*(dF/dp_i).                             *
 *             Cost: 2 evaluations of F(t,y,y',p).                *
 *             Accuracy: 1st order.                               *
 *                                                                *
 * rhomax <= -1: If (ratio > abs(rhomax)), use forward differences*
 *               to estimate the sum (dF/dy)*w_i + (dF/dy')*w'_i, *
 *               and the term pbar_i*(dF/dp_i).                   *
 *               Otherwise, choose delta = MIN(deltay, deltap)    *
 *               and use 1 forward difference to estimate the     *
 *               residual of the scaled sensitivity equations.    *
 *             Cost: 1 or 2 evaluations of F(t,y,y',p).           *
 *             Accuracy: 1st order.                               *
 *                                                                *
 *----------------------------------------------------------------*
 *                                                                *
 * plist   is a pointer to an array of integers. These integers   *
 *             identify the order in which the parameters are to  *
 *             be studied.                                        *
 *         For plist equal to NULL: w_i is the scaled sensitivity *
 *             of y with respect to p_i, for i = 1,...,Ns.        *
 *         For plist NOT equal to NULL: w_i is the scaled         *
 *             sensitivity of y with respect to p_j, where        *
 *             j = plist[i-1] for i = 1,...,Ns.                   *
 *         For example: let Ns = 2, plist[0] = 5, plist[1] = 3.   *
 *         w_1 is the scaled sensitivity of y with respect to     *
 *             parameter 5.                                       *
 *         w_2 is the scaled sensitivity of y with respect to     *
 *             parameter 3.                                       *
 *                                                                *
 * Note: The tolerance values may be changed in between calls to  *
 *       IDASolve for the same problem. These values refer to     *
 *       (*rtol) and either (*atol), for a scalar absolute        *
 *       tolerance, or the components of atol, for a vector       *
 *       absolute tolerance.                                      *
 *                                                                * 
 * If successful, SensIDAMalloc returns a pointer to initialized  *
 * problem memory. This pointer should be passed to IDA. If       *
 * an initialization error occurs, SensIDAMalloc prints an error  *
 * message to the file specified by errfp and returns NULL.       *
 *                                                                *
 ******************************************************************/

void *SensIDAMalloc(integer Ny, integer Ns, integer Ntotal, ResFn res,
		    void *rdata, real t0, N_Vector y0, N_Vector yp0, 
		    int itol, real *rtol, void *atol, N_Vector id, 
		    N_Vector constraints, FILE *errfp, boole optIn, 
		    long int iopt[], real ropt[], void *machEnv,
		    real p[], real pbar[], void *plist, real rhomax); 


/******************************************************************
 *                                                                *
 * Function : SensIDAFree                                         *
 *----------------------------------------------------------------*
 * SensIDAFree frees the problem memory ida_mem allocated by      *
 * SensIDAMalloc. Its only argument is the pointer ida_mem        *
 * returned by SensIDAMalloc.                                     *
 *                                                                *
 ******************************************************************/

void SensIDAFree(void *ida_mem);


/***************************************************************
 *                                                             *
 * Function : SensSetVecAtol                                   *
 *-------------------------------------------------------------*
 * SensSetVecAtol initializes the vector absolute error        *
 * tolerances of the Ns scaled sensitivity vectors to be the   *
 * same as the vector absolute error tolerances for the DAE    *
 * problem F(t,y,y',p) = 0.                                    *
 *                                                             *
 * atol is the pointer to the vector of absolute error         *
 *         tolerances for the DAE problem F(t,y,y',p) = 0.     *
 *                                                             *
 * Ns     is the number of sensitivity vectors to be solved.   *
 *                                                             *
 ***************************************************************/

void SensSetVecAtol(N_Vector atol, integer Ns);


/***************************************************************
 *                                                             *
 * Function : SensSetId                                        *
 *-------------------------------------------------------------*
 * This routine identifies each component of the Ns sensitivity*
 * vectors as a differential or algebraic variable.            *
 *                                                             *
 * id     is a pointer to a vector of size Ntotal = Ny*(1+Ns). *
 *                                                             *
 * Ns     is the number of sensitivity vectors to be solved.   *
 *                                                             *
 ***************************************************************/

void SensSetId(N_Vector id, integer Ns);


/***************************************************************
 *                                                             *
 * Function : SensInitZero                                     *
 *-------------------------------------------------------------*
 * This routine initializes the trailing Ny*Ns components of a *
 * vector to zero.                                             *
 *                                                             *
 * z      is a pointer to a vector of size Ntotal = Ny*(1+Ns). *
 *                                                             *
 * Ns     is the number of sensitivity vectors to be solved.   *
 ***************************************************************/

void SensInitZero(N_Vector z, integer Ns);


/***************************************************************
 *                                                             *
 * Type : SensData                                             *
 *-------------------------------------------------------------*
 * The type SensData is a type for the block of data required  *
 * for computing the sensitivity of the DAE solution with      *
 * respect to its parameters. This block of data is created    *
 * by a user call to SensIDAMalloc and is not normally seen    *
 * by the user.                                                *
 *                                                             *
 * ida_mem is the pointer to CVODE memory returned by          *
 *           SensIDAMalloc.                                    *
 *                                                             *
 * res_ptr is the pointer to the function in F(t,y,y',p).      *
 *                                                             *
 * rdata  is the pointer to user data that will be passed to   *
 *           the user's F function every time F is called.     *
 *                                                             *
 * p      is a pointer to an array of real parameters.         *
 *                                                             *
 * pbar   is a pointer to an array of nonzero, real parameters *
 *           that are used to scale the sensitivity equations. *
 *                                                             *
 * rhomax is a real value used for selecting a finite          *
 *           difference formula to estimate and scale the      *
 *           sensitivity of the DAE residual with respect to   *
 *           various parameters.                               *
 *                                                             *
 * Ny     is the total number of differential and algebraic    *
 *           variables in the DAE system F(t,y,y',p) = 0.      *
 *                                                             *
 * Ns     is the number of sensitivity vectors to be solved.   *
 *                                                             *
 * restemp is temporary storage (of size Ny) for estimating    *
 *           sensitivity derivatives using finite differences. *
 *                                                             *
 * ytemp is temporary storage (of size Ny) for estimating      *
 *           sensitivity derivatives using finite differences. *
 *                                                             *
 * yptemp is temporary storage (of size Ny) for estimating     *
 *           sensitivity derivatives using finite differences. *
 *                                                             *
 ***************************************************************/

typedef struct {
  void *ida_mem;
  int (*res_ptr)(integer, real, N_Vector, N_Vector, N_Vector, void *);
  void *res_data;
  real *p;
  real *pbar;
  int *plist;
  real rhomax;
  integer Ny;
  integer Ns; 
  N_Vector restemp;
  N_Vector ytemp;
  N_Vector yptemp;
} *SensData;


#endif

#ifdef __cplusplus
}
#endif
