/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:27:52 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This file (companion of nvector.h) contains definitions 
 * needed for the initialization of vector operations in Fortran.
 * -----------------------------------------------------------------
 */


#ifndef _FNVECTOR_H
#define _FNVECTOR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H
#include <sundials/sundials_config.h>
#endif

/* SUNDIALS solver IDs */

#define FCMIX_CVODE   1
#define FCMIX_IDA     2
#define FCMIX_KINSOL  3
#define FCMIX_ARKODE  4

#ifdef __cplusplus
}
#endif

#endif
