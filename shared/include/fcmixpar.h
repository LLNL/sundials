/*******************************************************************
 *                                                                 *
 * File          : fcmixpar.h                                      *
 * Programmers   : Alan C. Hindmarsh @ LLNL                        *
 * Version of    : 19 August 2002                                  *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This file contains machine-dependent definitions of parameters  *
 * that characterize the treatment of externals by the Fortran     *
 * compiler, for use with Fortran/C interface routines.            *
 *                                                                 *
 * Definitions:                                                    *
 *                                                                 *
 *   CRAY          If the system is a Cray, set CRAY to 1, else    *
 *                 set it to 0. This is needed since the Cray has  *
 *                 non-standard include files, and the loader      *
 *                 expects names all upper case.                   *
 *                                                                 *
 *   UNDERSCORE    Set this to 1 if the system's Fortran77 expects *
 *                 C functions to have an underscore at the end of *
 *                 the name, and also the Fortran77 appends an     *
 *                 underscore on all Fortran subroutine names.     *
 *                 Otherwise set it to 0. If CRAY = 1, UNDERSCORE  *
 *                 is ignored.                                     *
 *  Some examples:                                                 *
 *         1 for Sun, DEC Alpha, SGI, IBM with -q extname, Meiko   *
 *         0 for HP, IBM default                                   *
 *                                                                 *
 * Declarations:                                                   *
 *                                                                 *
 *   F2C_machEnv  A global variable used to communicate the        *
 *                machine environment structure among the various  *
 *                routine in a Fortran - C interface. It is        *
 *                declared extern here and must defined in the     *
 *                Fortran - C interface implementation file of the *
 *                NVECTOR module.                                  *
 *                                                                 *
 *******************************************************************/  

#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _fcmixpar_h
#define _fcmixpar_h
  
#define CRAY        0
#define UNDERSCORE  1

#include "nvector.h"

extern M_Env F2C_machEnv;

#endif

#ifdef __cplusplus
}
#endif
