/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Implementation header file for the generic ARKSLS linear solver 
 * module.
 *---------------------------------------------------------------*/

#ifndef _ARKSPARSE_IMPL_H
#define _ARKSPARSE_IMPL_H

#include "arkode/arkode_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
 ARKSLS solver constants:

 ARKS_MSBJ   maximum number of steps between Jacobian evaluations
 ARKS_DGMAX  maximum change in gamma between Jacobian evaluations
---------------------------------------------------------------*/
#define ARKS_MSBJ  50
#define ARKS_DGMAX RCONST(0.2)

 
/*---------------------------------------------------------------
 Types: ARKSlsMemRec, ARKSlsMem

 ARKSlsMem is pointer to a ARKSlsMemRec structure.
---------------------------------------------------------------*/

typedef struct ARKSlsMemRec {

  ARKSlsSparseJacFn s_Jeval; /* user Jacobian evaluation routine 
                                to be called                   */
  void *s_Jdata;             /* user data passed to s_Jeval    */

  long int s_nje;            /* nje = no. of calls to s_Jeval  */

  long int s_last_flag;      /* last error return flag         */

  int s_first_factorize;     /* flag telling whether the first 
			        factorization needs to happen  */

  int s_nstlj;               /* time step of last J evaluation */

  SlsMat s_A;                /* A = M - gamma * df/dy          */

  SlsMat s_savedJ;           /* saved copy of Jacobian         */

  int sparsetype;            /* matrix type: compressed sparse column or row */

  void *s_solver_data;       /* struct for solver data         */

} *ARKSlsMem;


/*---------------------------------------------------------------
 Types: ARKSlsMassMemRec, ARKSlsMassMem

 ARKSlsMassMem is pointer to a ARKSlsMassMemRec structure.
---------------------------------------------------------------*/

typedef struct ARKSlsMassMemRec {

  ARKSlsSparseMassFn s_Meval; /* user mass matrix evaluation 
                                    routine to be called       */
  void *s_Mdata;              /* user data passed to s_Meval   */

  long int s_nme;             /* nme = no. of calls to s_Meval */

  long int s_last_flag;       /* last error return flag        */

  int s_first_factorize;      /* flag telling whether the first 
			         factorization needs to happen */

  SlsMat s_M;                 /* mass matrix structure         */

  SlsMat s_M_lu;              /* mass matrix for LU decomp     */

  int sparsetype;             /* matrix type: compressed sparse column or row */

  void *s_solver_data;        /* struct for solver data        */

} *ARKSlsMassMem;


/*---------------------------------------------------------------
 Prototypes of internal functions (none)
---------------------------------------------------------------*/
  

/*---------------------------------------------------------------
 Error Messages
---------------------------------------------------------------*/
#define MSGSP_ARKMEM_NULL     "Integrator memory is NULL."
#define MSGSP_BAD_NVECTOR     "A required vector operation is not implemented."
#define MSGSP_MEM_FAIL        "A memory request failed."
#define MSGSP_LMEM_NULL       "Linear solver memory is NULL."
#define MSGSP_MASSMEM_NULL    "Mass matrix solver memory is NULL."
#define MSGSP_ILL_INPUT       "Invalid input detected."
#define MSGSP_JAC_NOSET       "Jacobian evaluation function has not been set."
#define MSGSP_MASS_NOSET      "Mass matrix evaluation function has not been set."
#define MSGSP_JACFUNC_FAILED  "The Jacobian routine failed in an unrecoverable manner."
#define MSGSP_MASSFUNC_FAILED "The mass matrix routine failed in an unrecoverable manner."
#define MSGSP_PACKAGE_FAIL    "A call to an external package failed."

#ifdef __cplusplus
}
#endif

#endif
