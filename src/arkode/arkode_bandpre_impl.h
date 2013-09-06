/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Implementation header file for the ARKBANDPRE module.
---------------------------------------------------------------*/

#ifndef _ARKBANDPRE_IMPL_H
#define _ARKBANDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode_bandpre.h>
#include <sundials/sundials_band.h>
#include <sundials/sundials_direct.h>

/*---------------------------------------------------------------
 Type: ARKBandPrecData
---------------------------------------------------------------*/

typedef struct ARKBandPrecDataRec {

  /* Data set by user in ARKBandPrecInit */
  long int N;
  long int ml, mu;

  /* Data set by ARKBandPrecSetup */
  DlsMat savedJ;
  DlsMat savedP;
  long int *lpivots;

  /* Rhs calls */
  long int nfeBP;

  /* Pointer to arkode_mem */
  void *arkode_mem;

} *ARKBandPrecData;


/*---------------------------------------------------------------
 ARKBANDPRE error messages
---------------------------------------------------------------*/

#define MSGBP_MEM_NULL       "Integrator memory is NULL."
#define MSGBP_LMEM_NULL      "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBP_MEM_FAIL       "A memory request failed."
#define MSGBP_BAD_NVECTOR    "A required vector operation is not implemented."
#define MSGBP_PMEM_NULL      "Band preconditioner memory is NULL. ARKBandPrecInit must be called."
#define MSGBP_RHSFUNC_FAILED "The right-hand side routine failed in an unrecoverable manner."


#ifdef __cplusplus
}
#endif

#endif
