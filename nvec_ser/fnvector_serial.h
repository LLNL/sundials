/*
 * -----------------------------------------------------------------
 * $Revision: 1.11.2.1 $
 * $Date: 2005-01-24 21:49:02 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/shared/LICENSE.
 * -----------------------------------------------------------------
 * This file (companion of nvector_serial.h) contains the
 * definitions needed for the initialization of serial
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_SERIAL_H
#define _FNVECTOR_SERIAL_H

#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H
#include "sundials_config.h"
#endif

#if defined(F77_FUNC)

#define FNV_INITS F77_FUNC(fnvinits, FNVINITS)
#define FNV_FREES F77_FUNC(fnvfrees, FNVFREES)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITS fnvinits
#define FNV_FREES fnvfrees

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITS FNVINITS
#define FNV_FREES FNVFREES

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITS fnvinits_
#define FNV_FREES fnvfrees_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITS FNVINITS_
#define FNV_FREES FNVFREES_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FNV_INITS fnvinits__
#define FNV_FREES fnvfrees__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FNV_INITS FNVINITS__
#define FNV_FREES FNVFREES__

#endif

#endif
