/******************************************************************
 *                                                                *
 * File          : senspvode.h                                    *
 * Programmers   : Steven L. Lee and Alan C. Hindmarsh @ LLNL     *
 * Last Modified : 05 February 2001                               *
 *----------------------------------------------------------------*
 * This is the header file for the implementation of computing    *
 * sensitivity information as given in the file:                  *
 * ../source/senspvode.c                                          *
 *                                                                *
 ******************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _senspvode_h
#define _senspvode_h

#include "cvode.h"
#include "llnltyps.h" 
#include "senscvspgmr.h"

/******************************************************************
 *                                                                *
 * Function : SensCVodeMalloc                                     *
 *----------------------------------------------------------------*
 * SensCVodeMalloc initializes a block of sensitivity data.       *
 * It also allocates and initializes memory for a problem         *
 * to be solved by CVODE.                                         *
 *                                                                *
 * Ny      is the number of ODEs contained in y' = f(t,y,p).      *
 *                                                                *
 * Ns      is the number of sensitivity vectors to be solved.     *
 *                                                                *
 * Ntotal  is the total number of ODEs to be solved by CVODE.     *
 *         Normally Ntotal = Ny*(1+Ns), but this is not required. *
 *                                                                *
 * f       is the right hand side function in y' = f(t,y,p).      *          
 *                                                                *
 * t0      is the initial value of t.                             *
 *                                                                *
 * y0      is a vector of length Ntotal that contains the initial *
 *            values for the Ny state variables and Ny*Ns         *
 *            sensitivity variables at time t = t0.               *
 *                                                                *
 * lmm     is the type of linear multistep method to be used.     *
 *            The legal values are ADAMS and BDF (see previous    *
 *            description).                                       *
 *                                                                *
 * iter    is the type of iteration used to solve the nonlinear   *
 *            system that arises during each internal time step.  *
 *            The legal values are FUNCTIONAL and NEWTON.         *
 *                                                                *
 * itol    is the type of tolerances to be used.                  *
 *            The legal values are:                               *
 *               SS (scalar relative and absolute  tolerances),   *
 *               SV (scalar relative tolerance and vector         *
 *                   absolute tolerance).                         *
 *                                                                *
 * reltol  is a pointer to the relative tolerance scalar.         *
 *                                                                *
 * abstol  is a pointer to the absolute tolerance scalar or       *
 *            an N_Vector tolerance of absolute tolerances.       *
 *                                                                *
 * The parameters itol, reltol, and abstol define a vector of     *
 * error weights, ewt, with components                            *
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol)   (if itol = SS), or  *
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (if itol = SV).  *
 * This vector is used in all error and convergence tests, which  *
 * use a weighted RMS norm on all error-like vectors v.           *
 *                                                                *
 * f_data  is a pointer to a block of user data that will be      *
 *            passed to the user's f function every time f is     *
 *            called. This data block must contain a pointer to   *
 *            the real array p, described below.                  *
 *                                                                *
 * errfp   is the file pointer for an error file where all CVODE  *
 *            warning and error messages will be written. This    *
 *            parameter can be stdout (standard output), stderr   *
 *            (standard error), a file pointer (corresponding to  *
 *            a user error file opened for writing) returned by   *
 *            fopen, or NULL. If the user passes NULL, then all   *
 *            messages will be written to standard output.        *
 *                                                                *
 * optIn   is a flag indicating whether there are any optional    *
 *            inputs from the user in the arrays iopt and ropt.   *
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
 * machEnv is a pointer to machine environment-specific           *
 *            information.                                        *
 *                                                                *
 * p       is a pointer to the array of real parameters given in  *
 *            the ODE y' = f(t,y,p). This pointer must be stored  *
 *            in the user's f_data data block so that the         *
 *            parameters can be accessed through f_data and this  *
 *            pointer every time f is called.                     *
 *                                                                *
 * pbar    is a pointer to an array of nonzero, real values that  *
 *            are used to scale the Ns sensitivity vectors.       *
 *         NOTE: pbar can have the same values as p (if all the   *
 *            parameters are nonzero); however, they cannot be    *
 *            the same array.                                     *
 *                                                                *
 * rhomax  is a real value used for selecting a finite difference *
 *            formula to estimate the scaled sensitivity          *
 *            derivative                                          *
 *                                                                *
 *               w'_i = (df/dy)*w_i + pbar_i*(df/dp_i),           *
 *                                                                * 
 *            where 1 <= i <= Ns.                                 *
 *            See the following discussion for details concerning *
 *            the finite difference formulas for w'_i.            *
 *----------------------------------------------------------------*
 *** The scaled sensitivity derivative w'_i is estimated using    *
 *** difference formulas involving the function f(t,y,p) and      *
 *** forward (or backward) perturbations, such as y + deltay*w_i  *
 *** and p + deltap*pbar_i. Vectors y and p can be perturbed      *
 *** separately, or simultaneously if the same value is used for  *
 *** deltay and deltap.                                           *
 ***                                                              *
 *** The heuristics for selecting deltay and deltap are based on  *
 *** the machine unit roundoff (uround), the user-specified       *
 *** relative error tolerance (rtol), and a weighted root-mean-   *
 *** square norm of the scaled sensitivity vector (wrmsnorm(w_i)):*
 ***                                                              *
 ***   deltap = sqrt(MAX(rtol, uround),                           *
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
 *   The values and meanings for rhomax are as follows:           *
 *                                                                *
 *   rhomax = 0: Choose delta = MIN(deltay, deltap) and use       *
 *               1 centered difference to estimate w'_i.          *
 *               Cost: 2 evaluations of f(t,y,p).                 *
 *               Accuracy: 2nd order.                             *
 *                                                                *
 *   0 < rhomax < 1: Use centered differences on each term:       *
 *                   (df/dy)*w_i, pbar_i*(df/dp_i).               *
 *                   Cost: 4 evaluations of f(t,y,p).             *
 *                   Accuracy: 2nd order.                         *
 *                                                                *
 *   rhomax >= 1: If  (ratio > rhomax), use centered differences  *
 *                on each term: (df/dy)*w_i, pbar_i*(df/dp_i).    *
 *                Otherwise, choose delta = MIN(deltay, deltap)   *
 *                and use 1 centered difference to estimate w'_i. *
 *                Cost: 2 or 4 evaluations of f(t,y,p).           *
 *                Accuracy: 2nd order.                            *
 *                                                                *
 *   -1 < rhomax < 0: Use forward differences on each term:       *
 *                    (df/dy)*w_i, pbar_i*(df/dp_i).              *
 *                    Cost: 2 evaluations of f(t,y,p).            *
 *                    Accuracy: 1st order.                        *
 *                                                                *
 *   rhomax <= -1: If  (ratio > abs(rhomax)), use forward         *
 *                 differences on each term: (df/dy)*w_i,         *
 *                 pbar_i*(df/dp_i). Otherwise, choose            *
 *                 delta = MIN(deltay, deltap) and use 1 forward  *
 *                 difference to estimate w'_i.                   *
 *                 Cost: 1 or 2 evaluations of f(t,y,p).          *
 *                 Accuracy: 1st order accuracy.                  *
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
 *       CVode for the same problem. These values refer to        *
 *       (*reltol) and either (*abstol), for a scalar absolute    *
 *       tolerance, or the components of abstol, for a vector     *
 *       absolute tolerance.                                      *
 *                                                                * 
 * If successful, SensCVodeMalloc returns a pointer to the        *
 * initialized problem memory. This pointer should be passed to   *
 * CVode. If an initialization error occurs, SensCVodeMalloc      *
 * prints an error message to the file specified by errfp and     *
 * returns NULL.                                                  *
 *                                                                *
 ******************************************************************/

void *SensCVodeMalloc(int Ny, int Ns, int Ntotal, RhsFn f, real t0, 
		      N_Vector y0, int lmm, int iter, int itol, 
		      real *reltol, void *abstol, void *f_data, 
		      FILE *errfp, boole optIn, long int iopt[], 
		      real ropt[], void *machEnv, real p[], 
		      real pbar[], void *plist, real rhomax);


/******************************************************************
 *                                                                *
 * Function : SensCVReInit                                        *
 *----------------------------------------------------------------*
 * SensCVReInit reinitializes CVode for the solution of a         *
 * problem, where a prior call to SensCVodeMalloc has been made   *
 * with the same problem size parameters Ny, Ns, and Ntotal.      *
 * SensCVReInit performs the same input checking and              *
 * initializations that SensCVodeMalloc does (except for Ny, Ns,  * 
 * and Ntotal). But it does no memory allocation, assuming that   *
 * the existing internal memory is sufficient for the new problem.*
 *                                                                *
 * The use of SensCVReInit requires that the maximum method order,*
 * maxord, is no larger for the new problem than for the problem  *
 * specified in the last call to SensCVodeMalloc. This condition  *
 * is automatically fulfilled if the multistep method parameter   *
 * lmm is unchanged (or changed from ADAMS to BDF) and the        *
 * default value for maxord is specified.                         *
 *                                                                *
 * The first argument to SensCVReInit is:                         *
 *                                                                *
 * cvode_mem = pointer to CVODE memory, from SensCVodeMalloc.     *
 *                                                                *
 * All the remaining arguments to SensCVReInit have names and     *
 * meanings identical to those of SensCVodeMalloc. Note that the  *
 * problem size parameters (Ny, Ns, Ntotal) are not passed as     *
 * arguments to SensCVReInit, as they are assumed unchanged since *
 * the SensCVodeMalloc call.                                      *
 *                                                                *
 * The return value of SensCVReInit is equal to SUCCESS = 0 if    *
 * there were no errors; otherwise it is a negative int equal to: *
 *  SensCVREI_NO_MEM     indicating cvode_mem was NULL, or        *
 *  SensCVREI_ILL_INPUT  indicating an input argument was illegal *
 *                      (including an attempt to increase maxord).*
 * In case of an error return, an error message is also printed.  *
 *                                                                *
 * Note: the reported workspace sizes iopt[LENRW] and iopt[LENIW] *
 * are unchanged from the values computed by SensCVodeMalloc, and *
 * so may be larger than would be computed for the new problem.   *
 ******************************************************************/

int SensCVReInit(void *cvode_mem, RhsFn f, real t0, N_Vector y0,
                 int lmm, int iter, int itol, real *reltol, void *abstol, 
		 void *f_data, FILE *errfp, boole optIn, long int iopt[], 
		 real ropt[], void *machEnv, real p[], real pbar[],
		 void *plist, real rhomax);


/* SensCVReInit return values: */

/* SUCCESS = 0  (Defined under CVode return values, but listed
                 here also for completeness)                      */
enum {SensCVREI_NO_MEM = -1, SensCVREI_ILL_INPUT = -2};


/******************************************************************
 *                                                                *
 * Function : SensCVodeFree                                       *
 *----------------------------------------------------------------*
 * SensCVodeFree frees the problem memory cvode_mem allocated by  *
 * SensCVodeMalloc. Its only argument is the pointer cvode_mem    *
 * returned by SensCVodeMalloc.                                   *
 *                                                                *
 ******************************************************************/

void SensCVodeFree(void *cvode_mem);


/***************************************************************
 *                                                             *
 * Function : SensSetVecAtol                                   *
 *-------------------------------------------------------------*
 * SensSetVecAtol initializes the vector absolute error        *
 * tolerances of the Ns scaled sensitivity vectors to be the   *
 * same as the vector absolute error tolerances for the ODEs   *
 * contained in y' = f(t,y,p).                                 *
 *                                                             *
 * atol is the pointer to the vector of absolute error         *
 *         tolerances for the CVODE problem y' = f(t,y,p).     *
 *                                                             *
 * Ns     is the number of sensitivity vectors to be solved.   *
 *                                                             *
 ***************************************************************/

void SensSetVecAtol(N_Vector atol, integer Ns);


/***************************************************************
 *                                                             *
 * Type : SensData                                             *
 *-------------------------------------------------------------*
 * The type SensData is a type for the block of data required  *
 * for computing the sensitivity of the ODE solution with      *
 * respect to its parameters. This block of data is created    *
 * by a user call to SensCVodeMalloc and is not normally seen  *
 * by the user.                                                *
 *                                                             *
 * cv_mem is the pointer to CVODE memory returned by           *
 *           SensCVodeMalloc.                                  *
 *                                                             *
 * f_ptr  is the pointer to the right hand side function in    *
 *           y' = f(t,y,p).                                    *
 *                                                             *
 * f_data is the pointer to user data that will be passed to   *
 *           the user's f function every time f is called.     *
 *                                                             *
 * p      is a pointer to an array of real parameters.         *
 *                                                             *
 * pbar   is a pointer to an array of nonzero, real parameters *
 *           that are used to scale the sensitivity ODEs.      *
 *                                                             *
 * rhomax is a real value used for selecting a finite          *
 *           difference formula to estimate the scaled         *
 *           sensitivity derivative w'.                        *
 *                                                             *
 * plist  is a pointer to an array of integers. These integers *
 *           identify the order in which parameters are to be  *
 *           studied.                                          *
 *                                                             *
 * Ny     is the number of ODEs contained in y' = f(t,y,p).    *
 *                                                             *
 * Ns     is the number of sensitivity vectors to be solved.   *
 *                                                             *
 * ftemp  is temporary storage for estimating sensitivity      *
 *           derivatives using finite differences.             *
 *                                                             *
 * ytemp  is temporary storage for estimating sensitivity      *
 *           derivatives using finite differences.             *
 *                                                             *
 ***************************************************************/

typedef struct {
  void *cv_mem;
  void (*f_ptr)(integer, real, N_Vector, N_Vector, void *);
  void *f_data;
  real *p;
  real *pbar;
  void *plist;
  real rhomax;
  integer Ny;
  integer Ns; 
  N_Vector ftemp;
  N_Vector ytemp;
} *SensData;


#endif

#ifdef __cplusplus
}
#endif
