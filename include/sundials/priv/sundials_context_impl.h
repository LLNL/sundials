/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * This is a 'private' header file and should not be used in user
 * code. It is subject to change without warning.
 * !!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * -----------------------------------------------------------------
 * SUNDIALS context class implementation.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CONTEXT_IMPL_H
#define _SUNDIALS_CONTEXT_IMPL_H

#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

struct SUNContext_
{
  SUNProfiler profiler;
  sunbooleantype own_profiler;
  SUNLogger logger;
  sunbooleantype own_logger;
  SUNErrCode last_err;
  SUNErrHandler err_handler;
  SUNComm comm;
};

#ifdef __cplusplus
}
#endif
#endif
