/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the SUNAdaptController_NOOP module.
 * -----------------------------------------------------------------*/

#ifndef _SUNADAPTCONTROLLER_NOOP_H
#define _SUNADAPTCONTROLLER_NOOP_H

#include <stdio.h>
#include <sundials/sundials_adaptcontroller.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ------------------------------------------
 * No-op implementation of SUNAdaptController
 * ------------------------------------------ */

/* ------------------
 * Exported Functions
 * ------------------ */

SUNDIALS_EXPORT
SUNAdaptController SUNAdaptController_NoOp(SUNContext sunctx);
SUNDIALS_EXPORT
SUNAdaptController_Type SUNAdaptController_GetType_NoOp(SUNAdaptController C);

#ifdef __cplusplus
}
#endif

#endif  /* _SUNADAPTCONTROLLER_PID_H */
