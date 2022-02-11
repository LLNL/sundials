/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Implementation header file for ARKODE's relaxation (in time) functionality.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_RELAX_IMPL_H
#define _ARKODE_RELAX_IMPL_H

#include <stdarg.h>
#include <arkode/arkode.h>

/* ==================== *
 * Relaxation Constants *
 * ==================== */

#define ARK_RELAX_LRW   5
#define ARK_RELAX_LIW  12   /* int, ptr, etc */

/* ========================= *
 * Relaxation Data Structure *
 * ========================= */

typedef struct ARKodeRelaxMemRec* ARKodeRelaxMem;

struct ARKodeRelaxMemRec
{
  ARKRelaxFn    rfn;        /* relaxation function               */
  ARKRelaxJacFn rjac;       /* relaxation Jacobian               */
  realtype      est;
  realtype      rcur;       /* current relaxation function value */
  realtype      tol;        /* nonlinear solve tolerance         */
  int           max_iters;  /* nonlinear solve max iterations    */
};

/* ==================== *
 * Relaxation Functions *
 * ==================== */

int arkSetRelaxFn(void* arkode_mem, ARKRelaxFn rfn, ARKRelaxJacFn rjac);

int arkRelax(void* arkode_mem, realtype* gam);

/* int arkRelaxFree(ARKodeRelaxMem* relax_mem); */
/* int arkRelaxPrint(ARKodeRelaMem relax_mem, FILE* outfile); */

#endif
