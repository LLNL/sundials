/*
 * -----------------------------------------------------------------
 * $Revision: 1.9 $
 * $Date: 2005-01-24 22:28:41 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the Fortran interface include file for the BAND
 * preconditioner (CVBANDPRE).
 * -----------------------------------------------------------------
 */

/*
 * ==============================================================================
 *
 *                  FCVBP Interface Package
 *
 * The FCVBP Interface Package is a package of C functions which,
 * together with the FCVODE Interface Package, support the use of the
 * CVODE solver (serial version) with the CVBANDPRE preconditioner module,
 * for the solution of ODE systems in a mixed Fortran/C setting.  The
 * combination of CVODE and CVBANDPRE solves systems dy/dt = f(t,y) with either
 * the SPGMR (scaled preconditioned GMRES) or SPBCG (scaled preconditioned
 * Bi-CGSTAB) method for the linear systems that arise, and with a banded
 * difference quotient Jacobian-based preconditioner.
 * 
 * The user-callable functions in this package, with the corresponding
 * CVODE and CVBBDPRE functions, are as follows: 
 *   FCVBPINIT  interfaces to CVBandPrecAlloc and CVSpgmr/CVSpbcg
 *   FCVBPOPT   accesses optional outputs
 *   FCVBPFREE  interfaces to CVBandPrecFree
 * 
 * In addition to the Fortran right-hand side function FCVFUN, the
 * user may (optionally) supply a routine FCVJTIMES which is called by 
 * the interface function FCVJtimes of type CVSpgmrJtimesFn or CVSpbcgJtimesFn.
 * (The names of all user-supplied routines here are fixed, in order to
 * maximize portability for the resulting mixed-language program.)
 * 
 * Important note on portability.
 * In this package, the names of the interface functions, and the names of
 * the Fortran user routines called by them, appear as dummy names
 * which are mapped to actual values by a series of definitions in the
 * header file fcvbp.h.
 * 
 * ==============================================================================
 * 
 *               Usage of the FCVODE/FCVBP Interface Packages
 * 
 * The usage of the combined interface packages FCVODE and FCVBP requires
 * calls to seven to ten interface functions, and one or two user-supplied
 * routines which define the problem to be solved and indirectly define
 * the preconditioner.  These function calls and user routines are
 * summarized separately below.
 * 
 * Some details are omitted, and the user is referred to the CVODE user document 
 * for more complete information.
 * 
 * (1) User-supplied right-hand side routine: FCVFUN
 * The user must in all cases supply the following Fortran routine
 *       SUBROUTINE FCVFUN (T, Y, YDOT)
 *       DIMENSION Y(*), YDOT(*)
 * It must set the YDOT array to f(t,y), the right-hand side of the ODE
 * system, as function of T = t and the array Y = y.  Here Y and YDOT
 * are distributed vectors.
 * 
 * (2) Optional user-supplied Jacobian-vector product routine: FCVJTIMES
 * As an option, the user may supply a routine that computes the product
 * of the system Jacobian J = df/dy and a given vector v.  If supplied, it
 * must have the following form:
 *       SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, EWT, WORK, IER)
 *       DIMENSION V(*), FJV(*), Y(*), FY(*), EWT(*), WORK(*)
 * Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
 * compute the product vector Jv, where the vector v is stored in V, and store
 * the product in FJV.  On return, set IER = 0 if FCVJTIMES was successful,
 * and nonzero otherwise.
 * 
 * (3) Initialization:  FNVINITS, FCVMALLOC, FCVBPINIT.
 * 
 * (3.1) To initialize the serial vector specification, the user must make 
 * the following call:
 *        CALL FNVINITS(NEQ, IER)
 * where NEQ is the problem size and IER is a return completion flag.
 * Possible values for IER are 0 = success, -1 = failure.
 * 
 * (3.2) To set various problem and solution parameters and allocate
 * internal memory for CVODE, make the following call:
 *       CALL FCVMALLOC(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL, INOPT,
 *      1               IOPT, ROPT, IER)
 * The arguments are:
 * T0     = initial value of t
 * Y0     = array of initial conditions
 * METH   = basic integration method: 1 = Adams (nonstiff), 2 = BDF (stiff)
 * ITMETH = nonlinear iteration method: 1 = functional iteration, 2 = Newton iter.
 * IATOL  = type for absolute tolerance ATOL: 1 = scalar, 2 = array
 * RTOL   = relative tolerance (scalar)
 * ATOL   = absolute tolerance (scalar or array)
 * INOPT  = optional input flag: 0 = none, 1 = inputs used
 * IOPT   = array of length 40 for integer optional inputs and outputs
 *          (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 * ROPT   = array of length 40 for real optional inputs and outputs
 *          The optional inputs are MAXORD, MXSTEP, MXHNIL, SLDET, H0, HMAX,
 *          HMIN, stored in IOPT(1), IOPT(2), IOPT(3), IOPT(14), ROPT(1),
 *          ROPT(2), ROPT(3), respectively.  If any of these optional inputs
 *          are used, set the others to zero to indicate default values.
 *          The optional outputs are NST, NFE, NSETUPS, NNI, NCFN, NETF, QU, QCUR,
 *          LENRW, LENIW, NOR, HU, HCUR, TCUR, TOLSF, stored in IOPT(4) .. IOPT(13),
 *          IOPT(15), ROPT(4) .. ROPT(7), resp.  See the CVODE manual for details. 
 * IER    = return completion flag.  Values are 0 = success, and -1 = failure.
 *          See printed message for details in case of failure.
 * 
 * (3.3) To allocate memory and initialize data associated with the CVBANDPRE
 * preconditioner, make the following call:
 *       CALL FCVBPINIT(NEQ, MU, ML, IER)
 * 
 * The arguments are:
 * NEQ       = problem size
 * MU, ML    = upper and lower half-bandwidths of the band matrix that 
 *             is retained as an approximation of the Jacobian.
 * IER       = return completion flag: IER=0: success, IER<0: and error occured
 *
 * (3.4A) To specify the SPBCG linear solver with the CVBANDPRE preconditioner,
 * make the following call
 *       CALL FCVBPSPBCG(IPRETYPE, MAXL, DELT, IER)
 * 
 * The arguments are:
 * IPRETYPE  = preconditioner type: 
 *            0 = none
 *            1 = left only
 *            2 = right only
 *            3 = both sides.
 * MAXL      = maximum Krylov subspace dimension; 0 indicates default.
 * DELT      = linear convergence tolerance factor; 0.0 indicates default.
 * IER       = return completion flag: IER=0: success, IER<0: ans error occured
 *
 * (3.4B) To specify the SPGMR linear solver with the CVBANDPRE preconditioner,
 * make the following call
 *       CALL FCVBPSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)
 * 
 * The arguments are:
 * IPRETYPE  = preconditioner type: 
 *            0 = none
 *            1 = left only
 *            2 = right only
 *            3 = both sides.
 * IGSTYPE   = Gram-schmidt process type: 0 = modified G-S, 1 = classical G-S.
 * MAXL      = maximum Krylov subspace dimension; 0 indicates default.
 * DELT      = linear convergence tolerance factor; 0.0 indicates default.
 * IER       = return completion flag: IER=0: success, IER<0: ans error occured
 *
 * (3.5A) To specify whether Bi-CGSTAB should use the supplied FCVJTIMES or the 
 * internal finite difference approximation, make the call
 *        CALL FCVSPBCGSETJAC(FLAG, IER)
 * where FLAG=0 for finite differences approxaimtion or
 *       FLAG=1 to use the supplied routine FCVJTIMES
 *
 * (3.5B) To specify whether GMRES should use the supplied FCVJTIMES or the 
 * internal finite difference approximation, make the call
 *        CALL FCVSPGMRSETJAC(FLAG, IER)
 * where FLAG=0 for finite differences approxaimtion or
 *       FLAG=1 to use the supplied routine FCVJTIMES
 * 
 * (4) The integrator: FCVODE
 * Carrying out the integration is accomplished by making calls as follows:
 *       CALL FCVODE (TOUT, T, Y, ITASK, IER)
 * The arguments are:
 * TOUT  = next value of t at which a solution is desired (input)
 * T     = value of t reached by the solver on output
 * Y     = array containing the computed solution on output
 * ITASK = task indicator: 0 = normal mode (overshoot TOUT and interpolate)
 *         1 = one-step mode (return after each internal step taken)
 * IER   = completion flag: 0 = success, values -1 ... -8 are various
 *         failure modes (see CVODE User Guide).
 * The current values of the optional outputs are available in IOPT and ROPT.
 * 
 * (5) Optional outputs: FCVBPOPT
 * Optional outputs specific to the SPGMR/SPBCG solver are NPE, NLI, NPS, NCFL,
 * LRW, and LIW, stored in IOPT(16) ... IOPT(21), respectively.
 * To obtain the optional outputs associated with the CVBANDPRE module, make
 * the following call:
 *       CALL FCVBPOPT(LENRPW, LENIPW, NFE)
 * The arguments returned are:
 * LENRPW = length of real preconditioner work space, in realtype words.
 *          This size is local to the current processor.
 * LENIPW = length of integer preconditioner work space, in integer words.
 *          This size is local to the current processor.
 * NGE    = number of f(t,y) evaluations for CVBANDPRE
 * 
 * (6) Computing solution derivatives: FCVDKY
 * To obtain a derivative of the solution (optionally), of order up to
 * the current method order, make the following call:
 *       CALL FCVDKY (T, K, DKY)
 * The arguments are:
 * T   = value of t at which solution derivative is desired
 * K   = derivative order (0 .le. K .le. QU)
 * DKY = array containing computed K-th derivative of y on return
 * 
 * (7) Memory freeing: FCVBPFREE, FCVFREE, and FNVFREES
 *   To the free the internal memory created by the calls to FNVINITS,
 * FCVMALLOC, and FCVBPINIT, make the following calls, in this order:
 *       CALL FCVBPFREE
 *       CALL FCVFREE
 *       CALL FNVFREES
 * 
 * ==============================================================================
 */

#ifndef _FCVBP_H
#define _FCVBP_H

/* header files  */

#include "nvector.h"       /* definition of type N_Vector */
#include "sundialstypes.h" /* definition of type realtype */

/* Definitions of interface function names */

#if defined(F77_FUNC)

#define FCV_BPINIT  F77_FUNC(fcvbpinit, FCVBPINIT)
#define FCV_BPSPBCG F77_FUNC(fcvbpspbcg, FCVBPSPBCG)
#define FCV_BPSPGMR F77_FUNC(fcvbpspgmr, FCVBPSPGMR)
#define FCV_BPOPT   F77_FUNC(fcvbpopt, FCVBPOPT)
#define FCV_BPFREE  F77_FUNC(fcvbpfree, FCVBPFREE)

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_LOWER)

#define FCV_BPINIT  fcvbpinit
#define FCV_BPSPBCG fcvbpspbcg
#define FCV_BPSPGMR fcvbpspgmr
#define FCV_BPOPT   fcvbpopt
#define FCV_BPFREE  fcvbpfree

#elif defined(SUNDIALS_UNDERSCORE_NONE) && defined(SUNDIALS_CASE_UPPER)

#define FCV_BPINIT  FCVBPINIT
#define FCV_BPSPBCG FCVBPSPBCG
#define FCV_BPSPGMR FCVBPSPGMR
#define FCV_BPOPT   FCVBPOPT
#define FCV_BPFREE  FCVBPFREE

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_LOWER)

#define FCV_BPINIT  fcvbpinit_
#define FCV_BPSPBCG fcvbpspbcg_
#define FCV_BPSPGMR fcvbpspgmr_
#define FCV_BPOPT   fcvbpopt_
#define FCV_BPFREE  fcvbpfree_

#elif defined(SUNDIALS_UNDERSCORE_ONE) && defined(SUNDIALS_CASE_UPPER)

#define FCV_BPINIT  FCVBPINIT_
#define FCV_BPSPBCG FCVBPSPBCG_
#define FCV_BPSPGMR FCVBPSPGMR_
#define FCV_BPOPT   FCVBPOPT_
#define FCV_BPFREE  FCVBPFREE_

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_LOWER)

#define FCV_BPINIT  fcvbpinit__
#define FCV_BPSPBCG fcvbpspbcg__
#define FCV_BPSPGMR fcvbpspgmr__
#define FCV_BPOPT   fcvbpopt__
#define FCV_BPFREE  fcvbpfree__

#elif defined(SUNDIALS_UNDERSCORE_TWO) && defined(SUNDIALS_CASE_UPPER)

#define FCV_BPINIT  FCVBPINIT__
#define FCV_BPSPBCG FCVBPSPBCG__
#define FCV_BPSPGMR FCVBPSPGMR__
#define FCV_BPOPT   FCVBPOPT__
#define FCV_BPFREE  FCVBPFREE__

#endif

/* Declarations for global variables, shared among various routines */

void *CVBP_Data;

#endif
