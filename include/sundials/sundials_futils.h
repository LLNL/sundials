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
 * SUNDIALS Fortran 2003 interface utility definitions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_FUTILS_H
#define _SUNDIALS_FUTILS_H

#include <stdio.h>
#include <sundials/sundials_config.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Create a file pointer with the given file name and mode. */
SUNDIALS_EXPORT FILE* SUNDIALSFileOpen(const char* filename, const char* modes);

/* Close a file pointer with the given file name. */
SUNDIALS_EXPORT void SUNDIALSFileClose(FILE* fp);


#ifdef __cplusplus
}
#endif

#endif
