#ifndef _SUNLINSOL_KLU_INTERNAL_H
#define _SUNLINSOL_KLU_INTERNAL_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <klu.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Interfaces to match 'sunindextype' with the correct KLU types/functions */
#if defined(SUNDIALS_INT64_T)
#define sun_klu_symbolic      klu_l_symbolic
#define sun_klu_numeric       klu_l_numeric
#define sun_klu_common        klu_l_common
#define sun_klu_analyze       klu_l_analyze
#define sun_klu_factor        klu_l_factor
#define sun_klu_refactor      klu_l_refactor
#define sun_klu_rcond         klu_l_rcond
#define sun_klu_condest       klu_l_condest
#define sun_klu_defaults      klu_l_defaults
#define sun_klu_free_symbolic klu_l_free_symbolic
#define sun_klu_free_numeric  klu_l_free_numeric
#elif defined(SUNDIALS_INT32_T)
#define sun_klu_symbolic      klu_symbolic
#define sun_klu_numeric       klu_numeric
#define sun_klu_common        klu_common
#define sun_klu_analyze       klu_analyze
#define sun_klu_factor        klu_factor
#define sun_klu_refactor      klu_refactor
#define sun_klu_rcond         klu_rcond
#define sun_klu_condest       klu_condest
#define sun_klu_defaults      klu_defaults
#define sun_klu_free_symbolic klu_free_symbolic
#define sun_klu_free_numeric  klu_free_numeric
#else  /* incompatible sunindextype for KLU */
#error  Incompatible sunindextype for KLU
#endif

#if defined(SUNDIALS_DOUBLE_PRECISION)
#else
#error  Incompatible realtype for KLU
#endif

/*
 * -----------------------------------------------------------------
 * PART I: KLU implementation of SUNLinearSolver
 *
 * The KLU implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     last_flag -- last error return flag from internal setup/solve
 *     first_factorize -- flag indicating whether the factorization 
 *       has ever been performed
 *     Symbolic -- KLU storage structure for symbolic 
 *       factorization components
 *     Numeric -- KLU storage structure for numeric factorization
 *        components
 *     Common -- storage structure for common KLU solver 
 *        components
 *     klu_solver -- ptr to KLU function to handle CSR/CSC.
 *        We create a typedef for this type of function pointer
 *        to suppress compiler warning messages about sunindextype 
 *        vs internal KLU index types.
 * -----------------------------------------------------------------
 */
typedef sunindextype (*KLUSolveFn)(sun_klu_symbolic*, sun_klu_numeric*,
                                   sunindextype, sunindextype,
                                   double*, sun_klu_common*);
 
struct _SUNLinearSolverContent_KLU {
  long int         last_flag;
  int              first_factorize;
  sun_klu_symbolic *symbolic;
  sun_klu_numeric  *numeric;
  sun_klu_common   common;
  KLUSolveFn       klu_solver;
};

typedef struct _SUNLinearSolverContent_KLU *SUNLinearSolverContent_KLU;

#ifdef __cplusplus
}
#endif

#include <sunlinsol/sunlinsol_klu.h>
#endif
