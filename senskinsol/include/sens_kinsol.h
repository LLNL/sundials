
/**********************************************************************
 *                                                                    *
 * File        : sens_kinsol.h                                        *
 * Programmers : Keith E. Grant and Alan C. Hindmarsh @ LLNL          *
 *                                                                    *
 * Version of  : 07 Aug  2001: Added snormtol to SensKinLinInit       *
 *--------------------------------------------------------------------*
 * This is the header file for the sensitivity parameter structures   *
 * used by SensKINSOL                                                 *
 *                                                                    *
 **********************************************************************/

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _sens_kinparams_h
#define _sens_kinparams_h

#include <stdio.h>
#include "llnltyps.h"
#include "nvector.h"
#include "kinsol.h"


/**********************************************************************
 *                                                                    *
 * Types : struct SensKINParams                                       *
 *--------------------------------------------------------------------*
 * The type SensKINParams is type pointer to a structure containing   *
 * fields for parameter sensitivity analysis. An instance of this     *
 * structure will be allocated and initialized by the call to         *
 * SensKINMalloc which will in turn call KINMalloc and set a pointer  *
 * to the KINSOL memory structure. A pointer to the parameter array   *
 * (p) passed to SensKINMalloc should also be placed (by the user)    *
 * into the fdata (i.e. UserData) structure so that it can be passed  *
 * to the user object function. It is by this mechanism that both the *
 * object function and the parameter variation routines of SensKINSOL *
 * share access to the function parameters.                           *
 *                                                                    *
 **********************************************************************/

typedef struct
{
   integer np;       /* The number of sensitivity parameters            */
   integer diffrhs;  /* Selects RHS sensitivity differencing option     */
   real    diffscl;  /* Scale factor for the RHS difference increment   */
   real   *p;        /* Array of sensitvity parameters                  */
   real   *pbar;     /* Array of nominal magnitudes of sens. parameters */
   void   *kin_mem;  /* Pointer to the KINSOL memory structure          */

   /* Linear solver function for sensitivity problems                   */
   int (*lsolve) (KINMem kin_mem, N_Vector vv, N_Vector ww,
      real *res_norm);

} *SensKINParams;


/* Define Sensitivity RHS difference options */
enum {DIFF_FRWRD, DIFF_BKWRD, DIFF_CNTRD};

/* Define SensKINSOL return codes that are in     */
/* addition to the underlying KINSOL return codes */
enum { SENSKIN_NO_MEM=KINSOL_LSOLV_NO_MEM-1 };


/**********************************************************************
 * Function : SensKINMalloc                                           *
 *                                                                    *
 * SensKINMalloc allocates and initializes a sensitivity parameter    *
 * structure                                                          *
 *--------------------------------------------------------------------*
 *                                                                    *
 * Arguments...                                                       *
 *                                                                    *
 * Neq      size of vectors being handled by the current memory       *
 *          allocation call to SensKINMalloc                          *
 *                                                                    *
 * np       The number of sensitivity parameters (np >= 0)            *
 *                                                                    *
 * diffrhs  Determines the form of the RHS difference equation        *
 *             = DIFF_FRWRD = 0 : Forward difference                  *
 *             = DIFF_BKWRD = 1 : Backward difference                 *
 *             = DIFF_CNTRD = 2 : Centered difference                 *
 *                                                                    *
 * diffscl  Positive scale factor for the RHS difference increment.   *
 *          The nominal value for diffscl is 1.0                      *
 *                                                                    *
 * p        Array of sensitivity parameters of length np. The user    *
 *          should also place a pointer to this array within the      *
 *          f_data structure passed by KINSol to the object function. *
 *                                                                    *
 * pbar     Array of sensitivity parameter nominal magnitudes > 0     *
 *                                                                    *
 * msgfp    Pointer to the file for error message output              *
 *                                                                    *
 * machEnv  A pointer to machine environment-specific information.    *
 *          Pass NULL for the sequential case (see nvector.h).        *
 *                                                                    *
 *                                                                    *
 * Return values...                                                   *
 *                                                                    *
 * SensKINMalloc =   Pointer to sensitivity memory. NULL on error     *
 *                                                                    *
 **********************************************************************/

void *SensKINMalloc(integer Neq, integer np, integer diffrhs,
   real diffscl, real *p, real *pbar, FILE *msgfp, void *machEnv);


/**********************************************************************
 *                                                                    *
 * Function: SensKINSol                                               *
 *                                                                    *
 * This routine encapsulates a call to KINSol. The arguments are      *
 * formally the same as for KINSol except that a pointer to a         *
 * sensitivity memory structure is passed to SensKINsol instead of a  *
 * pointer to KINSol memory. An additional underlying difference is   *
 * that the user's f_data structure must contain a pointer to the     *
 * same parameter array p[] that was passed to SensKINMalloc          *
 * -------------------------------------------------------------------*
 *                                                                    *
 * sens_params    Pointer to a sensitivity memory structure that      *
 *                contains a pointer to the KINSol memory             *
 *                                                                    *
 * ...            See kinsol.h and the KINSol user's guide            *
 *                                                                    *
 **********************************************************************/

 int SensKINSol(void *sens_params, integer Neq, N_Vector uu,
         SysFn func, int globalstrategy, N_Vector uscale,
         N_Vector fscale, real fnormtol, real scsteptol,
         N_Vector constraints, boole optIn, long int iopt[],
         real ropt[], void *f_data);


/**********************************************************************
 * Function: SensKINInit                                              *
 *                                                                    *
 * This routine initializes KINSol system variables needed by         *
 * SensKINSol without calling the nonlinear solver. SensKINInit is    *
 * to be used when the solution uu of the nonlinear system F(uu) = 0  *
 * is already accurately known. At the cost of including some unused  *
 * variables (as noted below), the calling sequence to SensKINInit    *
 * has been kept the same as that to SensKINSol                       *
 *                                                                    *
 * Arguments...                                                       *
 *                                                                    *
 *    sens_params       Pointer to the sensitivity memory structure   *
 *                                                                    *
 *    Neq               The number of equations in the algebraic      *
 *                      system                                        *
 *                                                                    *
 *    uu                Known solution for the system func(uu)=0.     *
 *                                                                    *
 *    func              Objective function for the system func(uu)=0. *
 *                                                                    *
 *    globalstrategy    Unused, but set to zero in KINSOL memory      *
 *                                                                    *
 *    uscale            An array (type N_Vector) of diagonal elements *
 *                      of the scaling matrix for uu. The elements of *
 *                      uscale must be positive values. The scaling   *
 *                      matrix uscale should be chosen so that        *
 *                      uscale * uu (as a matrix multiplication)      *
 *                      should have all its components with roughly   *
 *                      the same magnitude when uu is close to a root *
 *                      of func.                                      *
 *                                                                    *
 *    fscale            An array (type N_Vector) of diagonal elements *
 *                      of the scaling matrix for func. the elements  *
 *                      of fscale must be positive values. The        *
 *                      scaling matrix fscale should be chosen so     *
 *                      that fscale * func(uu) (as a matrix           *
 *                      multiplication) should have all its           *
 *                      components with roughly the same magnitude    *
 *                      when uu is NOT too near a root of func.       *
 *                                                                    *
 *    fnormtol          A real (scalar) value containing the stopping *
 *                      tolerance on maxnorm( fscale * func(uu) ). If *
 *                      fnormtol is input as 0., then a default value *
 *                      of (uround) to the 1/3 power will be used.    *
 *                      uround is the unit roundoff for the machine   *
 *                      in use for the calculation. (see UnitRoundoff *
 *                      in the llnlmath module                        *
 *                                                                    *
 *    scsteptol         Unused by SensKINInit                         *
 *                                                                    *
 *    constraints       Unused by SensKINInit                         *
 *                                                                    *
 *    optIn             is a flag indicating whether optional inputs  *
 *                      from the user in the arrays iopt and ropt are *
 *                      to be used. pass FALSE to ignore all optional *
 *                      inputs and TRUE to use all optional inputs    *
 *                      that are present. Either choice does NOT      *
 *                      affect outputs in other positions of iopt     *
 *                      or ropt.                                      *
 *                                                                    *
 *    iopt              The user-allocated array (of size OPT_SIZE)   *
 *                      that will hold optional integer inputs and    *
 *                      outputs. The user can pass NULL if he/she     *
 *                      does not wish to use optional integer inputs  *
 *                      or outputs. If optIn is TRUE, the user should *
 *                      preset to 0 those locations for which default *
 *                      values are to be used. See Optional Inputs    *
 *                      and Outputs in kinsol.h                       *
 *                                                                    *
 *    ropt              The user-allocated array (of size OPT_SIZE)   *
 *                      that will hold optional real inputs and       *
 *                      outputs. The user can pass NULL if he/she     *
 *                      does not wish to use optional real inputs or  *
 *                      outputs. If optIn is TRUE, the user should    *
 *                      preset to 0.0 the optional input locations    *
 *                      for which default values are to be used. See  *
 *                      Optional Inputs and Outputs in kinsol.h       *
 *                                                                    *
 *    f_data            A pointer to work space for use by the user-  *
 *                      supplied function func. The space allocated   *
 *                      to f_data is allocated by the user's program  *
 *                      before the call to SensKINMalloc.             *
 *                                                                    *
 * Return values...                                                   *
 *                                                                    *
 * SensKINInit = -4   Sensitivity memory allocation failuse           *
 *             = -3   KINSol linear solver memory allocation failure  *
 *             = -2   KINSol input error (see error output file)      *
 *             = -1   KINSol memory allocation failure                *
 *             = +1   Success flag for F(uu)==0 satisfied             *
 *                                                                    *
 * Note: The above return values are compatible with those returned   *
 * by SensKINSol. If uu is not a solution, this will be returned as   *
 * an input error (-2) with an appropriate error message.             *
 *                                                                    *
 **********************************************************************/

int SensKINInit (void *sens_params, integer Neq, N_Vector uu,
        SysFn func, int globalstrategy, N_Vector uscale,
        N_Vector fscale, real fnormtol, real scsteptol,
        N_Vector constraints, boole optIn, long int iopt[],
        real ropt[], void *f_data);


/******************************************************************
 *                                                                *
 * Function : SensKINLinInit                                      *
 *                                                                *
 * This routine calls the linear solver setup routine to insure   *
 * that the preconditioning is current. It also allows the user   *
 * to reset the linear solver norm tolerance                      *
 *                                                                *
 *----------------------------------------------------------------*
 *                                                                *
 * Arguments...                                                   *
 *                                                                *
 * sens_params  Pointer to the sensitivity memory structure       *
 *                                                                *
 * snormtol     A real (scalar) value containing the stopping     *
 *              tolerance for the linear solver. Passing a zero   *
 *              value will result in a default value of           *
 *              sqrt(fnormtol), fnormtol as used for the          *
 *              nonlinear solution.                               *
 *                                                                *
 * Return Values...                                               *
 *                                                                *
 * SensKINLINInit = 0  Success                                    *
 *                = 1  Pointer to sensitivity memory NULL         *
 *                = 2  Pointer to KINSOL memory NULL              *
 *                = 3  Error in Preconditioning Setup             *
 *                                                                *
 ******************************************************************/

int SensKINLinInit (void *sens_params, real snormtol);


/******************************************************************
 *                                                                *
 * Function : SensKINDiff                                         *
 *                                                                *
 * What SensKINDIFF allows the user to do is to change the RHS    *
 * difference option and difference increment scaling factor      *
 * between the solutions for individual parameters. Use of        *
 * SensKINDiff is optional, since the user specifices values for  *
 * both diffrhs and diffscl in the call to SensKINMalloc          *
 *----------------------------------------------------------------*
 *                                                                *
 * Arguments...                                                   *
 *                                                                *
 * sens     Pointer to the SensKINSOL sensitivity structure       *
 *                                                                *
 *                                                                *
 * diffrhs  Determines the form of the RHS difference equation    *
 *             = DIFF_FRWRD = 0 : Forward difference              *
 *             = DIFF_BKWRD = 1 : Backward difference             *
 *             = DIFF_CNTRD = 2 : Centered difference             *
 *                                                                *
 * diffscl  Scale factor (positive) for the RHS difference        *
 *          increment. The nominal value for diffscl is 1.0       *
 *                                                                *
 * Return Values...                                               *
 *                                                                *
 * SensKINDiff = 0   Success                                      *
 *             = 1   Pointer to sensitivity memory NULL           *
 *             = 2   Difference scaling factor <= 0               *
 *             = 3   Differencing option not valid                *
 *                                                                *
 ******************************************************************/

int SensKINDiff (void *sens_params, integer diffrhs, real diffscl);


/******************************************************************
 *                                                                *
 * Function : SensKINLinSolve                                     *
 *                                                                *
 * This routine solves the linear system J w = v, where w is the  *
 * scaled sensitivity vector for the current parameter (ip) and   *
 * v = - pbar*partial(F)/partial(p).                              *
 *----------------------------------------------------------------*
 *                                                                *
 * Arguments...                                                   *
 *                                                                *
 * sens_params Pointer to the sensitivity memory structure        *
 *                                                                *
 * ip          Index of the current sensitivity parameter         *
 *                                                                *
 * ww          The sensitivity vector for parameter ip (returned) *
 *                                                                *
 *                                                                *
 * Return Values...                                               *
 *                                                                *
 * SensKINLinSolve:                                               *
 *                  = -10  Error in evaluating sensitivity RHS    *
 *                  =  -9  Pointer to sensitivity memory is NULL  *
 *                  =  -8  N_Vector for RHS is NULL               *
 *                  =  -7  N-Vector for solution is NULL          *
 *                  =  -6  Pointer to linear solver mem. is NULL  *
 *                  =  -5  Pointer to KINSOL memory is NULL       *
 *                  =  -4  Other unrecoverable failure            *
 *                  =  -3  Unrecoverable preconditioner error     *
 *                  =  -2  Unrecoverable ATimes routine error     *
 *                  =  -1  NULL pointer found by Krylov solver    *
 *                  =   0  A linear solver success story          *
 *                  =   1  Recoverable failure to converge        *
 *                  =   2  Recoverable preconditioner error       *
 *                  =   3  A singular matrix was encountered      *
 *                  =   4  Other recoverable failure              *
 *                  =   5  Out of range sens. parameter index     *
 *                                                                *
 ******************************************************************/

int SensKINLinSolve (void *sens_params, integer ip, N_Vector ww);


/**********************************************************************
 *                                                                    *
 * Function : SensKINFree                                             *
 *                                                                    *
 * This routine deallocates the sensitivity structure memory and      *
 * calls KINFree to deallocate the KINSOL memory structure            *
 *--------------------------------------------------------------------*
 *                                                                    *
 * Arguments...                                                       *
 *                                                                    *
 * sens     Pointer to the sensitivity structure to free              *
 *                                                                    *
 **********************************************************************/

void SensKINFree (void *sens_params);

#endif

#ifdef __cplusplus
}
#endif
