
/**********************************************************************
 *                                                                    *
 * File        : sens_kinsol.c                                        *
 * Programmers : Keith E. Grant and Alan C. Hindmarsh @ LLNL          *
 *                                                                    *
 * Version of  : 07 Aug  2001: Added snormtol to SensKinLinInit       *
 *--------------------------------------------------------------------*
 * This is the file of routines associated with the sensitivity       *
 * parameter structures used by SensKINSOL                            *
 *                                                                    *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llnltyps.h"
#include "llnlmath.h"
#include "nvector.h"
#include "kinsol.h"
#include "sens_kinsol.h"


/* Error Messages */

#define SENSKINSOL         "\nSensKINSOL-- "
#define MSG_SENS_MEMFAIL   SENSKINSOL "Sens memory allocation failed\n\n"
#define MSG_SENS_NULL      SENSKINSOL "Pointer to sens is NULL\n\n"
#define MSG_SENS_NOP       SENSKINSOL "Pointer to p values is NULL\n\n"
#define MSG_SENS_NOPBAR    SENSKINSOL "Pointer to pbar values is NULL\n\n"
#define MSG_SENS_BPBAR     SENSKINSOL "pbar <= 0 is illegal\n\n"
#define MSG_SENS_DIFFRHS   SENSKINSOL "diffrhs = %d is illegal\n\n"
#define MSG_SENS_DIFFSCL   SENSKINSOL "diffscl <= 0 is illegal\n\n"


/* Constants */

#define ZERO   RCONST(0.0)    /* real 0.0  */
#define ONE    RCONST(1.0)    /* real 1.0  */
#define POINT5 RCONST(0.5)    /* real 0.5  */
#define TWO    RCONST(2.0)    /* real 2.0  */


/* Private function prototypes */


/**********************************************************************
 *                                                                    *
 * Function : SensKINRHS                                              *
 *                                                                    *
 * This routine generates the linear sensitivity equation RHS by      *
 * using a difference quotient approximation for pbar*partial(F) /    *
 * partial(p).                                                        *
 *                                                                    *
***********************************************************************/

static int SensKINRHS (void *sens_params, integer ip, N_Vector vv);


/**********************************************************************
 * Function: SensKINInitialStop                                       *
 *                                                                    *
 * This routine checks the initial guess (uu) to see if the system    *
 * func(uu) = 0. is satisfied to the level of fnormtol                *
 *                                                                    *
 **********************************************************************/

static boole SensKINInitialStop (KINMem kin_mem);

/* Sensitivity Function Definitions */

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
   real diffscl, real *p, real *pbar, FILE *msgfp, void *machEnv)
{

   SensKINParams sens;
   integer ip;
   boole bad_pbar;
   FILE *fp;

   fp = (msgfp == NULL) ? stdout : msgfp;


   /* Allocate the sensitivity parameter structure */

   sens = (SensKINParams) malloc(sizeof(*sens));

   if ( sens == NULL) {
      fprintf(fp, MSG_SENS_MEMFAIL);
      return (NULL);
   }


   /* Allocate KINSOL memory */

   sens->kin_mem = KINMalloc(Neq, fp, machEnv);
   if ( sens->kin_mem == NULL ) {
      free (sens);
      return (NULL);
   }


   /* Return to the sensitivity structure setup */

   sens->np = MAX(np,0);
   sens->diffrhs = diffrhs;
   sens->diffscl = diffscl;
   sens->lsolve = NULL;
   sens->p       = p;
   sens->pbar    = pbar;


   /* Allow for the case with zero parameters */

   if  ( sens->np == 0 ) return (sens);


   /* Check that pointer to p values is not NULL */

   if ( p == NULL) {
      fprintf(fp, MSG_SENS_NOP);
      KINFree (sens->kin_mem);
      free(sens);
      return (NULL);
   }


   /* Check that pointer to pbar nominal magnitudes is not NULL */

   if ( pbar == NULL) {
      fprintf(fp, MSG_SENS_NOPBAR);
      KINFree (sens->kin_mem);
      free(sens);
      return (NULL);
   }


   /* Check that the value for diffrhs is valid */
   if ( diffrhs != DIFF_FRWRD &&
        diffrhs != DIFF_BKWRD &&
        diffrhs != DIFF_CNTRD ) {

      fprintf(fp, MSG_SENS_DIFFRHS, diffrhs);
      KINFree (sens->kin_mem);
      free(sens);
      return (NULL);
   }


   /* Check that diffscl is positive */
   if ( diffscl <= ZERO ) {
      fprintf(fp, MSG_SENS_DIFFSCL);
      KINFree (sens->kin_mem);
      free(sens);
      return (NULL);
   }


   /* Check that the values of pbar are positive */

   bad_pbar = FALSE;
   for (ip=0; ip<np; ip++) bad_pbar = bad_pbar || pbar[ip] <= ZERO;
   if ( bad_pbar ) {
      fprintf(fp, MSG_SENS_BPBAR);
      KINFree (sens->kin_mem);
      free(sens);
      return (NULL);
   }

   return ( (void *) sens);
}


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
         real ropt[], void *f_data)
 {

   SensKINParams sens;
   integer ret;

   sens = (SensKINParams) sens_params;
   if ( sens == NULL) {
      fprintf(stdout, MSG_SENS_NULL);
      return (SENSKIN_NO_MEM);
   }

   ret = KINSol(sens->kin_mem, Neq, uu, func, globalstrategy,
            uscale, fscale, fnormtol,scsteptol, constraints,
            optIn, iopt, ropt, f_data);

   return (ret);

 }

#define HALF   RCONST(0.5)                    /* real 0.5   */
#define ONETHIRD RCONST(.3333333333333333)    /* real 1/3   */

#define PRINTFL_DEFAULT 0

#define f1norm          (kin_mem->kin_f1norm)
#define fnorm           (kin_mem->kin_fnorm)
#define fval            (kin_mem->kin_fval)
#define nbcf            (kin_mem->kin_nbcf)
#define nbktrk          (kin_mem->kin_nbktrk)
#define nfe             (kin_mem->kin_nfe)
#define nni             (kin_mem->kin_nni)
#define nnilpre         (kin_mem->kin_nnilpre)
#define printfl         (kin_mem->kin_printfl)
#define sqrt_relfunc    (kin_mem->kin_sqrt_relfunc)
#define uround          (kin_mem->kin_uround)


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
        real ropt[], void *f_data)
{
   SensKINParams sens;
   KINMem kin_mem;
   FILE *fp;
   boole ioptBad;
   boole ioptExists;
   boole roptExists;


   /* SensKINInit is unique in having to check both sensitivity */
   /* argument errors and KINSOL argument errors. Here we set   */
   /* up the error necessary error messages for this capability */

   const char myself[]="SensKINInit: ";
   const char emsg_sens_null[] = "\n%s Null pointer to sensitivity memory\n\n";
   const char emsg_kinm_null[] = "\n%s Null pointer to KINSol memory\n\n";
   const char emsg_linm_null[] = "\n%s Null pointer to lsolver memory\n\n";
   const char emsg_neq_bad[]   = "\n%s Neq > 0 is required instead of %d\n\n";
   const char emsg_uu_null[]   = "\n%s Null pointer to uu\n\n";
   const char emsg_func_null[] = "\n%s Null pointer to object function\n\n";
   const char emsg_uscl_null[] = "\n%s Null pointer to uscale\n\n";
   const char emsg_fscl_null[] = "\n%s Null pointer to fscale\n\n";
   const char emsg_uscl_npos[] = "\n%s One or more uscale values <= 0\n\n";
   const char emsg_fscl_npos[] = "\n%s One or more fscale values <= 0\n\n";
   const char emsg_fntol_neg[] = "\n%s fnormtol >= 0 required; was %12.4e\n\n";
   const char emsg_optin_bad[] = "\n%s Optin=%d is illegal; must be %d or %d\n\n";
   const char emsg_opt_null[]  = "\n%s Either iopt or ropt must be nonnull\n\n";
   const char emsg_iopt_out[]  = "\n%s 0 <= iopt[PRINTFL] <= 3 required\n\n";
   const char emsg_uu_nsoln[]  = "\n%s UU was not a solution to F(UU)==0\n\n";


   /* Check the pointer to the sensitivity memory structure */

   sens = (SensKINParams) sens_params;
   if ( sens == NULL) {
      fprintf(stderr, emsg_sens_null, myself);
      return (SENSKIN_NO_MEM);
   }


   /* Check the pointer from sensitivity memory to KINSOL */
   /* memory. This is initialized by SensKINMalloc.       */

   kin_mem = (KINMem) sens->kin_mem;
   if (kin_mem == NULL) {
     fprintf(stderr, emsg_kinm_null, myself);
     return(KINSOL_NO_MEM);
   }


   /* We now have access to the file pointer for error      */
   /* messages that was passed to SensKINMalloc. Subsequent */
   /* error message will use that file rather than stderr.  */

   fp = kin_mem->kin_msgfp;


   /* Check that the pointer to linear solver memory was set */

   if(kin_mem->kin_lmem == NULL) {
      fprintf(fp, emsg_linm_null, myself);
      return(KINSOL_LSOLV_NO_MEM);
   }


   if (Neq <= 0) {
     fprintf(fp, emsg_neq_bad, myself, Neq);
     return(KINSOL_INPUT_ERROR);
   }

   if (uu==NULL) {
      fprintf(fp, emsg_uu_null, myself);
      return(KINSOL_INPUT_ERROR);
   }

   if (func == NULL) {
      fprintf(fp, emsg_func_null, myself);
      return(KINSOL_INPUT_ERROR);
   }


  if (uscale == NULL)  {
    fprintf(fp, emsg_uscl_null, myself);
    return(KINSOL_INPUT_ERROR);
  }

  if (N_VMin(uscale) <= ZERO){
    fprintf(fp, emsg_uscl_npos, myself);
    return(KINSOL_INPUT_ERROR);
  }

  if (fscale == NULL)  {
    fprintf(fp, emsg_fscl_null, myself);
    return(KINSOL_INPUT_ERROR);
  }

  if (N_VMin(fscale) <= ZERO){
    fprintf(fp, emsg_fscl_npos, myself);
    return(KINSOL_INPUT_ERROR);
  }

  if (fnormtol < ZERO) {
    fprintf(fp, emsg_fntol_neg, myself, fnormtol);
    return(KINSOL_INPUT_ERROR);
  }

  if ((optIn != FALSE) && (optIn != TRUE)) {
    fprintf(fp, emsg_optin_bad, myself, optIn, FALSE, TRUE);
    return(KINSOL_INPUT_ERROR);
  }

  if ((optIn) && (iopt == NULL) && (ropt == NULL)) {
    fprintf(fp, emsg_opt_null, myself);
    return(KINSOL_INPUT_ERROR);
  }

  kin_mem->kin_ioptExists = ioptExists = (iopt != NULL);
  kin_mem->kin_roptExists = roptExists = (ropt != NULL);


  /* Copy the input parameters into KINSol state */

  kin_mem->kin_uu = uu;
  kin_mem->kin_func = func;
  kin_mem->kin_f_data = f_data;
  kin_mem->kin_globalstrategy = 0;
  kin_mem->kin_fnormtol = fnormtol;
  kin_mem->kin_scsteptol = ZERO;
  kin_mem->kin_iopt = iopt;
  kin_mem->kin_ropt = ropt;
  kin_mem->kin_constraints = NULL;
  kin_mem->kin_uscale = uscale;
  kin_mem->kin_fscale = fscale;
  kin_mem->kin_precondflag = FALSE;  /* set to the correct state in KINSpgmr */


  /* Check/set the value of the function norm tolerance */

  if(fnormtol <= ZERO )  fnormtol = RPowerR(uround,ONETHIRD);


  /* Handle the remaining optional inputs */

  printfl = PRINTFL_DEFAULT;

  if(ioptExists && optIn) {
    if ( iopt[PRINTFL] > 0 && iopt[PRINTFL] < 4 ) printfl = iopt[PRINTFL];
    ioptBad = (iopt[PRINTFL] < 0  || iopt[PRINTFL] > 3 );
    if(ioptBad){
       fprintf(fp, emsg_iopt_out, myself);
      return(KINSOL_INPUT_ERROR);
    }
  }

  if ( printfl > 0 ) fprintf(fp,"fnormtol  used: %12.3g \n",fnormtol);

  sqrt_relfunc = RSqrt(uround);



  /* Initialize all the counters */
  nfe = nnilpre = nni  = nbcf = nbktrk = 0;

  /* Initialize optional output locations in iopt, ropt */

  if (ioptExists) {
    iopt[NFE]    = iopt[NNI] = iopt[NNI] = 0;
    iopt[NBKTRK] = 0;
  }

  if (roptExists) {
    ropt[FNORM] = ZERO;
  }


  /* SensKINInitialStop returns "TRUE" if the system func(uu) = 0 is */
  /* satisfied by the initial guess uu                               */


  if ( SensKINInitialStop(kin_mem) ) {

    /* initialize the L2 norms of f for the linear iteration steps  */
    fnorm = N_VWL2Norm(fval,fscale);
    f1norm = HALF * fnorm * fnorm;
    if(printfl >0)fprintf(fp,
      "KINSolInit nni= %4ld  fnorm= %26.16g  nfe=%6ld \n",nni, fnorm, nfe);

    /* Problem has been successfully initialized */

    return (KINSOL_SUCCESS);

  } else {

    /* UU was not a solution to F(UU)==0 */
    fprintf(fp, emsg_uu_nsoln, myself);
    return (KINSOL_INPUT_ERROR);

  }

 }


/* Define Readability Constants */

#define func            (kin_mem->kin_func)
#define f_data          (kin_mem->kin_f_data)
#define fnormtol        (kin_mem->kin_fnormtol)
#define fscale          (kin_mem->kin_fscale)
#define iopt            (kin_mem->kin_iopt)
#define ioptExists      (kin_mem->kin_ioptExists)
#define lsetup          (kin_mem->kin_lsetup)
#define msgfp           (kin_mem->kin_msgfp)
#define Neq             (kin_mem->kin_Neq)
#define res_norm        (kin_mem->kin_res_norm)
#define uu              (kin_mem->kin_uu)
#define vtemp1          (kin_mem->kin_vtemp1)
#define vtemp2          (kin_mem->kin_vtemp2)


/******************************************************************
 *                                                                *
 * Function : SensKINLinInit                                      *
 *                                                                *
 * This routine calls the linear solver setup routine to insure   *
 * that the preconditioning is current. It also allows the user   *
 * to reset the linear solver norm tolerance.                     *
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

int SensKINLinInit (void *sens_params, real snormtol)
{
   SensKINParams sens;
   KINMem kin_mem;
   int ret;

   enum {OKAY, BAD_SENS, BAD_KINMEM, ERR_LSETUP};

   sens = (SensKINParams) sens_params;
   if ( sens == NULL ) {
      return (BAD_SENS);
   }

   kin_mem = (KINMem) sens->kin_mem;
   if ( kin_mem == NULL ) {
      return (BAD_KINMEM);
   }

   /* Update preconditioning setup via the KINSOL interface routine */
   if ( lsetup != NULL ) {
      ret = lsetup (kin_mem);
      if ( ret != 0 ) {
         return (ERR_LSETUP);
      }
   }

   fnormtol = ( snormtol > 0.0 ) ? snormtol : sqrt(fnormtol);

   return (OKAY);

}


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

int SensKINDiff (void *sens_params, integer diffrhs, real diffscl)
{
   SensKINParams sens;
   KINMem kin_mem;

   enum {OKAY, BAD_SENS, BAD_DIFFSCL, BAD_DIFFRHS};

   sens = (SensKINParams) sens_params;

   /* Check that the sens pointer is valid */
   if ( sens == NULL ) {
      fprintf(stdout, MSG_SENS_NULL);
      return (BAD_SENS);
   }

   kin_mem = (KINMem) sens->kin_mem;

   /* Check that the value for diffrhs is valid */
   if ( diffrhs != DIFF_FRWRD &&
        diffrhs != DIFF_BKWRD &&
        diffrhs != DIFF_CNTRD ) {

      fprintf(msgfp, MSG_SENS_DIFFRHS, diffrhs);
      free(sens);
      return (BAD_DIFFRHS);
   }
   sens->diffrhs = diffrhs;


   /* Check that diffscl is positive */
   if ( diffscl <= ZERO ) {
      fprintf(msgfp, MSG_SENS_DIFFSCL);
      free(sens);
      return (BAD_DIFFSCL);
   }
   sens->diffscl = diffscl;

   return (OKAY);
}


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

int SensKINLinSolve (void *sens_params, integer ip, N_Vector ww)
{

   SensKINParams sens;
   KINMem kin_mem;
   int ret;
   N_Vector vv;

   /* Definitions for return values from the Sensitivity Krylov      */
   /* Solver (SKS) beyond the generic Krylov returns. These return   */
   /* must be kept compatible with values within the generic Krylov  */
   /* solver set (e.g., kinspgmr.h) and those added by a             */
   /* particular linear solver wrapper.                              */

   enum {
      SKS_UNRECOV_RHSERR=-10,   /* Error in evaluating sens. RHS     */
      SKS_UNRECOV_SENSNULL,     /* Pointer to sensitiv. mem is NULL  */
      SKS_UNRECOV_RHSNULL,      /* N_Vector for RHS is NULL          */
      SKS_UNRECOV_SOLNNULL,     /* N_Vector for solution is NULL     */
      SKS_UNRECOV_KLSMEMNULL,   /* Pointer to lin solver mem is NULL */
      SKS_UNRECOV_KINMEMNULL,   /* Pointer to KINSOL memory is NULL  */
      SKS_RECOV_BADPARAMIX=5    /* Out of range sens. param. index   */
   };

   /* Verify that the pointer to sensitivity memory is nonNULL */
   sens = (SensKINParams) sens_params;
   if ( sens == NULL )  return (SKS_UNRECOV_SENSNULL);


   /* Verify that the pointer to KINSOL memory is nonNULL */
   kin_mem = (KINMem) sens->kin_mem;
   if ( kin_mem == NULL ) return (SKS_UNRECOV_KINMEMNULL);


   /* Verify that the pointer to the solution N_Vector is nonNULL */
   if ( ww == NULL ) return (SKS_UNRECOV_SOLNNULL);

   /* Verify that the current sens. parameter index is in range */
   if ( ip < 0 || ip >= sens->np ) return (SKS_RECOV_BADPARAMIX);

   /* As does KINSOL, use unew as work space for vv, the RHS vector */
   vv = kin_mem->kin_unew;
   if (vv == NULL ) return (SKS_UNRECOV_RHSNULL);

   /* Provide information for optional diagnostics */
   if ( ioptExists && iopt[PRINTFL] > 3 ) {
      fprintf(msgfp,
         "Solving sensitivity problem for parameter %d\n", ip);
   }


   /* Calculate the RHS of the linear sensitivity equation */
   ret = SensKINRHS (sens_params, ip, vv);
   if ( ret != 0) return (SKS_UNRECOV_RHSERR);


   /* Call the linear solver wrapper to solve J w = v  */
   ret = sens->lsolve(kin_mem, ww, vv, &res_norm);


   /* Update function evaluation statistics */
   if ( ioptExists ) iopt[NFE] = nfe;

   return (ret);
}


/******************************************************************
 *                                                                *
 * Function : SensKINFree                                         *
 *                                                                *
 * This routine call KINFree to deallocate KINSOL memory and then *
 * deallocates the sensitivity structure memory                   *
 *----------------------------------------------------------------*
 *                                                                *
 * Arguments...                                                   *
 *                                                                *
 * sens     Pointer to the sensitivity structure to free          *
 *                                                                *
 ******************************************************************/

void SensKINFree (void *sens_params)
{
   SensKINParams sens;

   sens = (SensKINParams) sens_params;

   KINFree (sens->kin_mem);

   free (sens);
}


/******************************************************************
 *                   ***** Local Functions *****                  *
 ******************************************************************/


/******************************************************************
 *                                                                *
 * Function : SensKINRHS                                          *
 *                                                                *
 * This routine generates the linear sensitivity equation RHS     *
 * by using a difference quotient approximation for               *
 * -pbar*partial(F)/partial(p).                                   *
 *----------------------------------------------------------------*
 *                                                                *
 * Arguments...                                                   *
 *                                                                *
 * sens_params  Pointer to the sensitivity memory structure       *
 *                                                                *
 * ip          Index of the current sensitivity parameter         *
 *                                                                *
 * vv          Right hand side vector evaluated by SensKINRHS     *
 *                                                                *
 * Return Values...                                               *
 *                                                                *
 * SenKINRHS   = 0  Sucsess                                       *
 *             = 1  Pointer to the Sensitivity was NULL           *
 *             = 2  Pointer to KINSOL memory was NULL             *
 *             = 3  Difference option not valid                   *
 *                                                                *
******************************************************************/

static int SensKINRHS (void *sens_params, integer ip, N_Vector vv)
{

   real sigma, sigma_inv;
   real sign, psave, rhs_norm;
   SensKINParams sens;
   KINMem   kin_mem;

   enum {OKAY, BAD_SENS, BAD_KINMEM, BAD_DIFFRHS};

   sens = (SensKINParams) sens_params;
   if ( sens == NULL ) {
      return (BAD_SENS);
   }

   kin_mem = (KINMem) sens->kin_mem;
   if ( kin_mem == NULL ) {
      return (BAD_KINMEM);
   }


   /* Calculate the difference increment */

   sign = ( sens->p[ip]  >= ZERO) ? ONE : -ONE;
   sigma = sign * sqrt_relfunc *sens->diffscl * sens->pbar[ip];
   sigma_inv = sens->pbar[ip] / sigma;

   if ( ioptExists && iopt[PRINTFL]==3 ) {
      fprintf (msgfp, "   Sensitivity parameter increment = %12.4e\n",
         sigma);
   }


   psave = sens->p[ip];    /* Save the parameter base value */


   /* Implement the chosen difference option */
   switch ( sens->diffrhs ) {

      case DIFF_FRWRD :       /* Forward difference */

         sens->p[ip] += sigma;
         func(Neq, uu, vtemp1, f_data);  nfe++;
         N_VLinearSum(-sigma_inv, vtemp1, sigma_inv, fval, vv);
         break;

      case DIFF_BKWRD :       /* Backward difference */

         sens->p[ip] -= sigma;
         func(Neq, uu, vtemp1, f_data);  nfe++;
         N_VLinearSum(-sigma_inv, vtemp1, sigma_inv, fval, vv);
         break;

      case DIFF_CNTRD :       /* Centered difference */

         sens->p[ip] -= POINT5*sigma;
         func(Neq, uu, vtemp1, f_data);  nfe++;
         sens->p[ip] += sigma;
         func(Neq, uu, vtemp2, f_data);  nfe++;
         N_VLinearSum(-sigma_inv, vtemp2, sigma_inv, vtemp1, vv);
         break;

      default:

         return (BAD_DIFFRHS);  /* Unknown differnce option */
   }

   sens->p[ip] = psave;       /* Restore the base parameter value */

   /*Diagnostic output of the norm of the scaled RHS derivatives */
   if ( ioptExists && iopt[PRINTFL]==3 ) {
      rhs_norm = N_VWL2Norm(vv, fscale);
      fprintf (msgfp,
         "   Sensitivity RHS scaled L2 norm  = %12.4e\n",
         rhs_norm);
   }

   return(OKAY);
}

/**********************************************************************
 * Function: SensKINInitialStop                                       *
 *                                                                    *
 * This routine checks the initial guess (uu) to see if the system    *
 * func(uu) = 0. is satisfied to the level of fnormtol                *
 *                                                                    *
 **********************************************************************/

static boole SensKINInitialStop (KINMem kin_mem)
{

  const char myself[]="SensKINInitialStop: ";
  real fmax;

  func(Neq, uu, fval, f_data);    nfe++;

  /* Compute the scaled max norm */
  N_VAbs (fval, vtemp1);
  N_VProd (fscale, vtemp1, vtemp1);
  fmax = N_VMaxNorm(vtemp1);

  if ( printfl > 1 ) fprintf (msgfp,
        " %s scaled f-norm for valid solution is %12.3g\n",
        myself, fmax);

  return ( fmax <= fnormtol );
}

