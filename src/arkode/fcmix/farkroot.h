/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 *---------------------------------------------------------------
 * This is the Fortran interface include file for the rootfinding
 * feature of ARKODE.
 *--------------------------------------------------------------*/

/*===============================================================
                   FARKROOT Interface Package

 The FARKROOT interface package allows programs written in 
 FORTRAN to use the rootfinding feature of the ARKODE solver 
 module.

 The user-callable functions constituting the FARKROOT package 
 are the following: FARKROOTINIT, FARKROOTINFO, and FARKROOTFREE. 
 The corresponding ARKODE subroutine called by each interface 
 function is given below.

   FARKROOT routine      ARKODE equivalent
   ----------------      -------------------
     FARKROOTINIT    ->   ARKodeRootInit
     FARKROOTINFO    ->   ARKodeGetRootInfo
     FARKROOTFREE    ->   ARKodeRootInit

 FARKROOTFN is a user-supplied subroutine defining the functions 
 whose roots are sought.

 ================================================================

              Usage of the FARKROOT Interface Package

 1. In order to use the rootfinding feature of the ARKODE package 
    the user must define the following subroutine:

    SUBROUTINE FARKROOTFN(T, Y, G, IPAR, RPAR, IER)

    The arguments are:
      T = independent variable value t  [realtype, input]
      Y = dependent variable array y  [realtype, input]
      G = function value array g(t,y)  [realtype, output]
      IPAR = user data array [long int, input/output]
      RPAR = user data array [realtype, input/output]
      IER = return flag [int, output]:
            0 if success
            non-zero if error

 2. After calling FARKMALLOC but prior to calling FARKODE, the 
    user must allocate and initialize memory for the FARKROOT 
    module by making the following call:

    CALL FARKROOTINIT(NRTFN, IER)

    The arguments are:
      NRTFN = total number of root functions  [int, input]
      IER   = return flag [int, output]:
                 0 if success
		-1 if ARKODE memory NULL 
               -11 if memory allocation error

 3. After calling FARKODE, to see whether a root was found, test 
    the FARKODE return flag IER.  The value IER = 2 means one or 
    more roots were found.

 4. If a root was found, and if NRTFN>1, then to determine which 
    root functions G(*) were found to have a root, make the 
    following call:
    
    CALL FARKROOTINFO(NRTFN, INFO, IER)

    The arguments are:
      NRTFN = total number of root functions [int, input]
      INFO = array of length NRTFN [int, output]:
                if G(i) has a root, then INFO(i) = 1
                if G(i) does not have a root, then INFO(i) = 0
      IER = return flag [int, output]
                0 if success
                negative if failure

 5. The total number of calls made to the root function 
    (FARKROOTFN), NGE, can be obtained from IOUT(13).

    Note: if the FARKODE/ARKODE memory block is reinitialized to 
    solve a different problem via a call to FARKREINIT, then the 
    counter NGE is reset to zero.

 6. To free the memory resources allocated by a prior call to 
    FARKROOTINIT make the following call:

    CALL FARKROOTFREE()

 For additional information, see the ARKODE documentation.
===============================================================*/

#ifndef _FARKROOT_H
#define _FARKROOT_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_ROOTINIT SUNDIALS_F77_FUNC(farkrootinit, FARKROOTINIT)
#define FARK_ROOTINFO SUNDIALS_F77_FUNC(farkrootinfo, FARKROOTINFO)
#define FARK_ROOTFREE SUNDIALS_F77_FUNC(farkrootfree, FARKROOTFREE)
#define FARK_ROOTFN   SUNDIALS_F77_FUNC(farkrootfn,   FARKROOTFN)

#else

#define FARK_ROOTINIT farkrootinit_
#define FARK_ROOTINFO farkrootinfo_
#define FARK_ROOTFREE farkrootfree_
#define FARK_ROOTFN   farkrootfn_

#endif

/* Prototypes of exported function */
void FARK_ROOTINIT(int *nrtfn, int *ier);
void FARK_ROOTINFO(int *nrtfn, int *info, int *ier);
void FARK_ROOTFREE(void);

/* Prototype of function called by ARKODE module */
int FARKrootfunc(realtype t, N_Vector y, 
		 realtype *gout, void *user_data);

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
   EOF
===============================================================*/
