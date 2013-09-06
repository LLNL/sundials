/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Implementation header file for the ARKBBDPRE module.
---------------------------------------------------------------*/

#ifndef _ARKBBDPRE_IMPL_H
#define _ARKBBDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode_bbdpre.h>
#include <sundials/sundials_band.h>


/*---------------------------------------------------------------
 Type: ARKBBDPrecData
---------------------------------------------------------------*/
typedef struct ARKBBDPrecDataRec {

  /* passed by user to ARKBBDPrecAlloc and used by PrecSetup/PrecSolve */
  long int mudq, mldq, mukeep, mlkeep;
  realtype dqrely;
  ARKLocalFn gloc;
  ARKCommFn cfn;

  /* set by ARKBBDPrecSetup and used by ARKBBDPrecSolve */
  DlsMat savedJ;
  DlsMat savedP;
  long int *lpivots;

  /* set by ARKBBDPrecAlloc and used by ARKBBDPrecSetup */
  long int n_local;

  /* available for optional output */
  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* pointer to arkode_mem */
  void *arkode_mem;

} *ARKBBDPrecData;


/*---------------------------------------------------------------
 ARKBBDPRE error messages
---------------------------------------------------------------*/

#define MSGBBD_MEM_NULL    "Integrator memory is NULL."
#define MSGBBD_LMEM_NULL   "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBBD_MEM_FAIL    "A memory request failed."
#define MSGBBD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBD_PMEM_NULL   "BBD peconditioner memory is NULL. ARKBBDPrecInit must be called."
#define MSGBBD_FUNC_FAILED "The gloc or cfn routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
