/******************************************************************
 *                                                                *
 * File          : fkinuatimes.c                                  *
 * Programmers   : Allan G Taylor and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 17 January 2001                                *
 *----------------------------------------------------------------*
 * Routine used to interface between a Fortran main and a user-   *
 * supplied A-times routine FATIMES (Jacobian A times vector v)   *
 ******************************************************************/

#include "kinsol.h"  /* prototypes of interfaces, global variables  */
#include "fkinsol.h" /* prototypes of interfaces, global variables  */


/*******************************************************************/
/*  C function KINUAtimes is to interface between                   *
 *  KINSpgmr/KINSpgmrATimes and FATIMES, the user-supplied Fortran  * 
 *  ATimes routine FATIMES.                                         *
 *******************************************************************/

int KINUAtimes(void *f_data, N_Vector v, N_Vector z, 
               boole *new_uu,  N_Vector uu)
{
 int retcode;
 real *vdata, *zdata, *uudata;

 vdata      = N_VDATA(v);
 zdata      = N_VDATA(z);
 uudata     = N_VDATA(uu);
 
 F_ATIMES(vdata, zdata, (int *)new_uu, uudata, &retcode);

return(retcode);
}
