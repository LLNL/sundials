/* -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNAdjointCheckpointScheme class declaration.
 * ----------------------------------------------------------------*/

#ifndef _SUNADJOINT_CHECKPOINTSCHEME_H
#define _SUNADJOINT_CHECKPOINTSCHEME_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef _SUNDIALS_STRUCT_ SUNAdjointCheckpointScheme_Ops_* SUNAdjointCheckpointScheme_Ops;
typedef _SUNDIALS_STRUCT_ SUNAdjointCheckpointScheme_* SUNAdjointCheckpointScheme;

typedef SUNErrCode (*SUNAdjointCheckpointSchemeNeedsSavingFn)(
  SUNAdjointCheckpointScheme, suncountertype step_num, suncountertype stage_num,
  sunrealtype t, sunbooleantype* yes_or_no);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeNeedsDeletingFn)(
  SUNAdjointCheckpointScheme, suncountertype step_num, suncountertype stage_num,
  sunrealtype t, sunbooleantype* yes_or_no);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeInsertVectorFn)(
  SUNAdjointCheckpointScheme, suncountertype step_num, suncountertype stage_num,
  sunrealtype t, N_Vector y);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeLoadVectorFn)(
  SUNAdjointCheckpointScheme, suncountertype step_num, suncountertype stage_num,
  sunbooleantype peek, N_Vector* yout, sunrealtype* tout);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeRemoveVectorFn)(
  SUNAdjointCheckpointScheme, suncountertype step_num, suncountertype stage_num,
  N_Vector* out);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeDestroyFn)(
  SUNAdjointCheckpointScheme* scheme);

typedef SUNErrCode (*SUNAdjointCheckpointSchemeEnableDenseFn)(
  SUNAdjointCheckpointScheme, sunbooleantype on_or_off);

struct SUNAdjointCheckpointScheme_Ops_;

struct SUNAdjointCheckpointScheme_;

/*
 * "static" base class methods
 */

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_NewEmpty(SUNContext sunctx,
                                               SUNAdjointCheckpointScheme*);

/*
 * Base class methods
 */

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetNeedsSavingFn(
  SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeNeedsSavingFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetNeedsDeletingFn(
  SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeNeedsDeletingFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetInsertVectorFn(
  SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeInsertVectorFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetLoadVectorFn(
  SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeLoadVectorFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetRemoveVectorFn(
  SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeRemoveVectorFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetDestroyFn(
  SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeDestroyFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetEnableDenseFn(
  SUNAdjointCheckpointScheme, SUNAdjointCheckpointSchemeEnableDenseFn);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_SetContent(SUNAdjointCheckpointScheme,
                                                 void* content);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_GetContent(SUNAdjointCheckpointScheme,
                                                 void** content);

/*
 * Virtual (overridable) base class methods
 */

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_NeedsSaving(SUNAdjointCheckpointScheme,
                                                  suncountertype step_num,
                                                  suncountertype stage_num,
                                                  sunrealtype t,
                                                  sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_NeedsDeleting(SUNAdjointCheckpointScheme,
                                                    suncountertype step_num,
                                                    suncountertype stage_num,
                                                    sunrealtype t,
                                                    sunbooleantype* yes_or_no);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_InsertVector(SUNAdjointCheckpointScheme,
                                                   suncountertype step_num,
                                                   suncountertype stage_num,
                                                   sunrealtype t, N_Vector state);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_LoadVector(
  SUNAdjointCheckpointScheme, suncountertype step_num, suncountertype stage_num,
  sunbooleantype peek, N_Vector* out, sunrealtype* tout);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_RemoveVector(SUNAdjointCheckpointScheme,
                                                   suncountertype step_num,
                                                   suncountertype stage_num,
                                                   N_Vector* out);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_Destroy(SUNAdjointCheckpointScheme*);

SUNDIALS_EXPORT
SUNErrCode SUNAdjointCheckpointScheme_EnableDense(SUNAdjointCheckpointScheme,
                                                  sunbooleantype on_or_off);

#ifdef __cplusplus
}
#endif

#endif /*_SUNADJOINT_CHECKPOINTSCHEME_H*/
