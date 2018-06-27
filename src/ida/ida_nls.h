/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------------------
 * This the header file for the IDA nonlinear solver interface.
 * This file will likely be absorbed into ida.h
 * ---------------------------------------------------------------------------*/

#ifndef _IDA_NLS_H
#define _IDA_NLS_H

#include "ida_impl.h"
#include "sunnls_newton.h" /* relace with sundial_nonlinearsolve.h */

/* exported functions */
int IDASetNonlinearSolver(void *ida_mem);

#endif
