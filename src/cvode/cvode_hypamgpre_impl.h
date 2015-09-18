/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Michael Wittman, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Implementation header file for the CVHYPRE_BOOMERAMG module.
 * -----------------------------------------------------------------
 */

#ifndef _CVHYPRE_BOOMERAMG_IMPL_H
#define _CVHYPRE_BOOMERAMG_IMPL_H

#include <cvode/cvode_hypamgpre.h>
#include <sundials/sundials_band.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE.h"
#include "_hypre_parcsr_mv.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Type: CVBoomerAMGData
 * -----------------------------------------------------------------
 */

typedef struct CVBoomerAMGDataRec {

  /* passed by user to CVBoomerAMGAlloc and used by PrecSetup/PrecSolve */

  long int mudq, mldq, mukeep, mlkeep;
  realtype dqrely;
  CVParCsrJacFn jacfn;

  /* set by CVBoomerAMGAlloc and used by CVBoomerAMGSetup */

  long int n_local;

  /* available for optional output */

  long int rpwsize;
  long int ipwsize;
  long int nge;
  
  int ilower, iupper, jlower, jupper, N;

  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
  HYPRE_Solver solver, precond;

  /* pointer to cvode_mem */
  void *cvode_mem;

} *CVBoomerAMGData;

/*
 * -----------------------------------------------------------------
 * CVHYPRE_BOOMERAMG error messages
 * -----------------------------------------------------------------
 */

#define MSGBBD_MEM_NULL    "Integrator memory is NULL."
#define MSGBBD_LMEM_NULL   "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBBD_MEM_FAIL    "A memory request failed."
#define MSGBBD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBD_PMEM_NULL   "BBD peconditioner memory is NULL. CVBoomerAMGInit must be called."
#define MSGBBD_FUNC_FAILED "The gloc or cfn routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
