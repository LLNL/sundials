/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos
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
 * SUNDIALS Fortran 2003 interface utility implementations.
 * -----------------------------------------------------------------*/

#include <string.h>
#include <sundials/priv/sundials_errors_impl.h>
#include <sundials/sundials_errors.h>

/* Create a file pointer with the given file name and mode. */
SUNErrCode SUNDIALSFileOpen(const char* filename, const char* mode, FILE** fp_out)
{
  SUNErrCode err = SUN_SUCCESS;
  FILE* fp       = *fp_out;

  if (filename)
  {
    if (!strcmp(filename, "stdout")) { fp = stdout; }
    else if (!strcmp(filename, "stderr")) { fp = stderr; }
    else { fp = fopen(filename, mode); }
  }

  if (!fp) { err = SUN_ERR_FILE_OPEN; }

  *fp_out = fp;
  return err;
}

/* Close a file pointer with the given file name. */
SUNErrCode SUNDIALSFileClose(FILE** fp_ptr)
{
  if (!fp_ptr) { return SUN_SUCCESS; }
  FILE* fp = *fp_ptr;
  if (fp && (fp != stdout) && (fp != stderr)) { fclose(fp); }
  return SUN_SUCCESS;
}
