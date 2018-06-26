#ifndef _IDA_NLS_H
#define _IDA_NLS_H

#include "ida_impl.h"
#include "sunnls_newton.h"

#define ONE     RCONST(1.0)    /* real 1.0    */
#define TWENTY  RCONST(20.0)   /* real 20.0   */

/* exported functions */
int IDASetNonlinearSolver(void *ida_mem);

#endif
