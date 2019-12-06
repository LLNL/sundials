/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2002-2018, Southern Methodist University and
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
 * -----------------------------------------------------------------
 * SUNDIALS Fortran 2003 interface utility implementations.
 * -----------------------------------------------------------------*/

#include <sundials/sundials_futils.h>

/* Create a file pointer with the given file name and mode. */
FILE* SUNDIALSFileOpen(const char* filename, const char* mode)
{
  return fopen(filename, mode);
}

/* Close a file pointer with the given file name. */
void SUNDIALSFileClose(FILE* fp)
{
  fclose(fp);
}
