/*
 * -----------------------------------------------------------------
 * $Revision: 1.20 $
 * $Date: 2004-07-22 22:54:43 $
 * ----------------------------------------------------------------- 
 * Programmers: Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California
 * Produced at the Lawrence Livermore National Laboratory
 * All rights reserved
 * For details, see sundials/cvode/LICENSE
 * -----------------------------------------------------------------
 * This is the header file for FCVODE, the Fortran interface to
 * the CVODE package.                                            
 * -----------------------------------------------------------------
 */

#ifndef _fcvode_h
#define _fcvode_h

/*
 * ===============================================================================
 *
 *                FCVODE Interface Package
 *
 * The FCVODE Interface Package is a package of C functions which support
 * the use of the CVODE solver, for the solution of ODE systems 
 * dy/dt = f(t,y), in a mixed Fortran/C setting.  While CVODE is written
 * in C, it is assumed here that the user's calling program and
 * user-supplied problem-defining routines are written in Fortran. 
 * This package provides the necessary interface to CVODE for both the
 * serial and the parallel NVECTOR implementations.
 * 
 * The user-callable functions, with the corresponding CVODE functions,
 * are as follows:
 * 
 *   FNVINITS and FNVINITP initialize serial and parallel vector
 *              computations, respectively
 * 
 *   FCVMALLOC  interfaces to CVodeCreate, CVodeSet*, and CVodeMalloc
 * 
 *   FCVREINIT  interfaces to CVReInit
 * 
 *   FCVDIAG    interfaces to CVDiag
 * 
 *   FCVDENSE   interfaces to CVDense
 * 
 *   FCVBAND    interfaces to CVBand
 * 
 *   FCVSPGMR, FCVSPGMRREINIT interface to CVSpgmr
 * 
 *   FCVODE     interfaces to CVode and CVodeGet*
 * 
 *   FCVDKY     interfaces to CVodeGetDky
 * 
 *   FCVFREE    interfaces to CVodeFree
 * 
 *   FNVFREES and FNVFREEP finalize serial and parallel vector
 *              computations, respectively
 * 
 * The user-supplied functions, each listed with the corresponding interface
 * function which calls it (and its type within CVODE), are as follows:
 *   FCVFUN    is called by the interface function FCVf of type RhsFn
 *   FCVDJAC   is called by the interface function FCVDenseJac of type CVDenseJacFn
 *   FCVBJAC   is called by the interface function FCVBandJac of type CVBandJacFn
 *   FCVPSOL   is called by the interface function FCVPSol of type CVSpgmrPrecSolveFn
 *   FCVPSET   is called by the interface function FCVPSet of type CVSpgmrPrecSetupFn
 *   FCVJTIMES is called by the interface function FCVJtimes of type CVSpgmrJacTimesVecFn
 * In contrast to the case of direct use of CVODE, and of most Fortran ODE
 * solvers, the names of all user-supplied routines here are fixed, in
 * order to maximize portability for the resulting mixed-language program.
 * 
 * Important note on portability.
 * In this package, the names of the interface functions, and the names of
 * the Fortran user routines called by them, appear as dummy names
 * which are mapped to actual values by a series of definitions in the
 * header file fcvode.h. These definition depend in turn on variables
 * SUNDIALS_UNDERSCORE_NONE and SUNDIALS_UNDERSCORE_TWO which can be 
 * set at the configuration stage.
 * 
 * ===============================================================================
 * 
 *                  Usage of the FCVODE Interface Package
 * 
 * The usage of FCVODE requires calls to six or seven interface
 * functions, depending on the method options selected, and one or more
 * user-supplied routines which define the problem to be solved.  These
 * function calls and user routines are summarized separately below.
 * 
 * Some details are omitted, and the user is referred to the user documents
 * on CVODE for more complete documentation.  Information on the
 * arguments of any given user-callable interface routine, or of a given
 * user-supplied function called by an interface function, can be found in
 * the documentation on the corresponding function in the CVODE package.
 * 
 * The number labels on the instructions below end with s for instructions
 * that apply to the serial version of CVODE only, and end with p for
 * those that apply to the parallel version only.
 * 
 * 
 * (1) User-supplied right-hand side routine: FCVFUN
 * The user must in all cases supply the following Fortran routine
 *       SUBROUTINE FCVFUN (T, Y, YDOT)
 *       DIMENSION Y(*), YDOT(*)
 * It must set the YDOT array to f(t,y), the right-hand side of the ODE 
 * system, as function of T = t and the array Y = y.  Here Y and YDOT
 * are distributed vectors.
 * 
 * (2) Optional user-supplied dense Jacobian approximation routine: FCVDJAC
 * As an option when using the DENSE linear solver, the user may supply a
 * routine that computes a dense approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FCVDJAC (NEQ, T, Y, FY, DJAC, EWT, H, WK1, WK2, WK3)
 *       DIMENSION Y(*), FY(*), EWT(*), DJAC(NEQ,*), WK1(*), WK2(*), WK3(*)
 * Typically this routine will use only NEQ, T, Y, and DJAC. It must compute
 * the Jacobian and store it columnwise in DJAC.
 * 
 * (3) Optional user-supplied band Jacobian approximation routine: FCVBJAC
 * As an option when using the BAND linear solver, the user may supply a
 * routine that computes a band approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FCVBJAC (NEQ, MU, ML, MDIM, T, Y, FY,
 *      1                    BJAC, EWT, H, WK1, WK2, WK3)
 *       DIMENSION Y(*), FY(*), EWT(*), BJAC(MDIM,*), WK1(*), WK2(*), WK3(*)
 * Typically this routine will use only NEQ, MU, ML, T, Y, and BJAC. 
 * It must load the MDIM by N array BJAC with the Jacobian matrix at the
 * current (t,y) in band form.  Store in BJAC(k,j) the Jacobian element J(i,j)
 * with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.
 * 
 * (4) Optional user-supplied Jacobian-vector product routine: FCVJTIMES
 * As an option when using the SPGMR linear solver, the user may supply a 
 * routine that computes the product of the system Jacobian J = df/dy and 
 * a given vector v.  If supplied, it must have the following form:
 *       SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, EWT, H, WORK, IER)
 *       DIMENSION V(*), FJV(*), Y(*), FY(*), EWT(*), WORK(*)
 * Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
 * compute the product vector Jv, where the vector v is stored in V, and store
 * the product in FJV.  On return, set IER = 0 if FCVJTIMES was successful,
 * and nonzero otherwise.
 * 
 * (5) Initialization:  FNVINITS / FNVINITP , FCVMALLOC, FCVREINIT
 * 
 * (5.1s) To initialize the serial machine environment, the user must make
 * the following call:
 *        CALL FNVINITS (NEQ, IER)
 * The arguments are:
 * NEQ     = size of vectors
 * IER     = return completion flag. Values are 0 = success, -1 = failure.
 * 
 * (5.1p) To initialize the parallel machine environment, the user must make 
 * the following call:
 *        CALL FNVINITP (NLOCAL, NGLOBAL, IER)
 * The arguments are:
 * NLOCAL  = local size of vectors on this processor
 * NGLOBAL = the system size, and the global size of vectors (the sum 
 *           of all values of NLOCAL)
 * IER     = return completion flag. Values are 0 = success, -1 = failure.
 * Note: If MPI was initialized by the user, the communicator must be
 * set to MPI_COMM_WORLD.  If not, this routine initializes MPI and sets
 * the communicator equal to MPI_COMM_WORLD.
 * 
 * (5.2) To set various problem and solution parameters and allocate
 * internal memory, make the following call:
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
 *          The optional inputs are MAXORD, MXSTEP, MXHNIL, SLDET, 
 *          MAXERRTESTFAILS, MAXNONLINITERS, MAXCONVFAILS, H0, HMAX,
 *          HMIN, TSTOP, NONLINCONVCOEF, stored in IOPT(1), IOPT(2), 
 *          IOPT(3), IOPT(14), IOPT(22), IOPT(23), IOPT(24), ROPT(1), ROPT(2), 
 *          ROPT(3), ROPT(8), ROPT(9), respectively.  If any of these optional 
 *          inputs are used, set the others to zero to indicate default values.
 *          The optional outputs are NST, NFE, NSETUPS, NNI, NCFN, NETF, QU, 
 *          QCUR, LENRW, LENIW, NOR, HU,HCUR, TCUR, TOLSF, stored in 
 *          IOPT(4) .. IOPT(13), IOPT(15), ROPT(4) .. ROPT(7), resp.  
 *          See the CVODE manual for details. 
 * IER    = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
 *          See printed message for details in case of failure.
 * 
 * (5.3) To re-initialize the CVODE solver for the solution of a new problem
 * of the same size as one already solved, make the following call:
 *       CALL FCVREINIT(T0, Y0, IATOL, RTOL, ATOL, INOPT, IOPT, ROPT, IER)
 * The arguments have the same names and meanings as those of FCVMALLOC,
 * except that NEQ, METH, and ITMETH  have been omitted from the argument list 
 * (being unchanged for the new problem).  
 * FCVREINIT performs the same initializations as FCVMALLOC, but does no memory 
 * allocation, using instead the existing internal memory created by the previous 
 * FCVMALLOC call.  The call to specify the linear system solution method may or 
 * may not be needed; see paragraph (6) below.
 * 
 * (6) Specification of linear system solution method.
 * In the case of a stiff system, the implicit BDF method involves the solution
 * of linear systems related to the Jacobian J = df/dy of the ODE system.
 * CVODE presently includes four choices for the treatment of these systems,
 * and the user of FCVODE must call a routine with a specific name to make the
 * desired choice.
 * 
 * (6.1) Diagonal approximate Jacobian.
 * This choice is appropriate when the Jacobian can be well approximated by
 * a diagonal matrix.  The user must make the call:
 *       CALL FCVDIAG(IER)
 * IER is an error return flag: 0 = success, -1 = memory failure.
 * There is no additional user-supplied routine.  Optional outputs specific
 * to the approximate diagonal Jacobian case are LRW and LIW, stored in
 * IOPT(16) and IOPT(17), respectively.  (See the CVODE manual for descriptions.)
 * 
 * (6.2) DENSE treatment of the linear system.
 * The user must make the call
 *       CALL FCVDENSE(NEQ, IER)
 * The argument is:
 * IER = error return flag: 0 = success , -1 = memory allocation failure,
 *                         -2 = illegal input. 
 * 
 * If the user program includes the FCVDJAC routine for the evaluation of the 
 * dense approximation to the Jacobian, the following call must be made
 *       CALL FCVDENSESETJAC(FLAG, IER)
 * The argument FLAG=0 specifies using the internal finite differences
 * approximation to the Jacobian, while FLAG=1 specifies that CVDJAC is
 * provided.
 * 
 *      Optional outputs specific to the DENSE case are LRW, LIW, and NJE
 * stored in IOPT(16), IOPT(17), and IOPT(18), respectively.  (See the CVODE
 * manual for descriptions.)
 * 
 * (6.3) BAND treatment of the linear system
 * The user must make the call
 *       CALL FCVBAND(NEQ, MU, ML, IER)
 * The arguments are:
 * MU  = upper bandwidth
 * ML  = lower bandwidth
 * IER = error return flag: 0 = success , -1 = memory allocation failure,
 *                         -2 = illegal input.     
 * 
 * If the user program includes the FCVBJAC routine for the evaluation of the 
 * band approximation to the Jacobian, the following call must be made
 *       CALL FCVBANDSETJAC(FLAG, IER)
 * The argument FLAG=0 specifies using the internal finite differences
 * approximation to the band Jacobian, while FLAG=1 specifies that CVBJAC is
 * provided.
 * 
 *      Optional outputs specific to the BAND case are LRW, LIW, and NJE
 * stored in IOPT(16), IOPT(17), and IOPT(18), respectively.  (See the CVODE
 * manual for descriptions.)
 * 
 * (6.4) SPGMR treatment of the linear systems.
 * For the Scaled Preconditioned GMRES solution of the linear systems,
 * the user must make the following call:
 *       CALL FCVSPGMR(IGSTYPE, MAXL, DELT, IER)              
 * 
 * The arguments are:
 * IPRETYPE = preconditioner type: 
 *              0 = none 
 *              1 = left only
 *              2 = right only
 *              3 = both sides
 * IGSTYPE  = Gram-schmidt process type: 
 *              1 = modified G-S
 *              2 = classical G-S.
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * DELT     = linear convergence tolerance factor; 0.0 indicates default.
 * IER      = error return flag: 
 *              0 = success
 *             -1 = memory allocation failure
 *             -2 = illegal input
 * 
 * If the user program includes the CVJTIMES routine for the evaluation of the 
 * Jacobian vector product, the following call must be made
 *       CALL FCVSPGMRSETJAC(FLAG, IER)
 * The argument FLAG=0 specifies using the internal finite differences
 * approximation to the Jacobian vector product, while FLAG=1 specifies that 
 * CVJTIMES is provided.
 * 
 * Usage of the user-supplied routine CVPSOL for solution of the preconditioner 
 * linear system is specified by calling
 *       CALL FCVSPGMRSETPSOL(FLAG, IER)
 * where FLAG=0 indicates no FCVPSOL (default) and FLAG=1 specifies using FCVPSOL.
 * The user-supplied routine FCVPSOL must be of the form:
 *       SUBROUTINE FCVPSOL (T, Y,FY, VT, GAMMA, EWT, H, DELTA, NFE, R, LR, Z, IER)
 *       DIMENSION Y(*), FY(*), VT(*), EWT(*), R(*), Z(*),
 * Typically this routine will use only NEQ, T, Y, GAMMA, R, LR, and Z.  It
 * must solve the preconditioner linear system Pz = r, where r = R is input, 
 * and store the solution z in Z.  Here P is the left preconditioner if LR = 1
 * and the right preconditioner if LR = 2.  The preconditioner (or the product
 * of the left and right preconditioners if both are nontrivial) should be an 
 * approximation to the matrix I - GAMMA*J (I = identity, J = Jacobian).
 * 
 * Usage of the user-supplied routine FCVPSET for construction of the preconditioner 
 * is specified by calling
 *       CALL FCVSPGMRSETPPSET(FLAG, IER)
 * where FLAG=0 indicates no FCVPSET (default) and FLAG=1 specifies using FCVPSET.
 * The user-supplied routine FCVPSET must be of the form:
 *       SUBROUTINE FCVPSET(T, Y, FY, JOK, JCUR, GAMMA, EWT, H, UROUND, 
 *      1                   NFE, V1, V2, V3, IER)
 *       DIMENSION Y(*), FY(*), EWT(*), V1(*), V2(*), V3(*) 
 * Typically this routine will use only NEQ, T, Y, JOK, and GAMMA. It must
 * perform any evaluation of Jacobian-related data and preprocessing needed
 * for the solution of the preconditioner linear systems by FCVPSOL.
 * The JOK argument allows for Jacobian data to be saved and reused:  If 
 * JOK = 0, this data should be recomputed from scratch.  If JOK = 1, a saved
 * copy of it may be reused, and the preconditioner constructed from it.
 * On return, set JCUR = 1 if Jacobian data was computed, and 0 otherwise.
 * Also on return, set IER = 0 if FCVPSET was successful, set IER positive if a 
 * recoverable error occurred, and set IER negative if a non-recoverable error
 * occurred.
 * 
 *      Optional outputs specific to the SPGMR case are LRW, LIW, NPE, NLI, NPS, 
 * and NCFL, stored in IOPT(16) ... IOPT(21), respectively.  (See the CVODE
 * manual for descriptions.)
 * 
 *      If a sequence of problems of the same size is being solved using the SPGMR
 * linear solver, then following the call to FCVREINIT, a call to the FCVSPGMRREINIT
 * routine may or may not be needed.  
 * 
 * (7) The integrator: FCVODE
 * Carrying out the integration is accomplished by making calls as follows:
 *       CALL FCVODE (TOUT, T, Y, ITASK, IER)
 * The arguments are:
 * TOUT  = next value of t at which a solution is desired (input)
 * T     = value of t reached by the solver on output
 * Y     = array containing the computed solution on output
 * ITASK = task indicator: 1 = normal mode (overshoot TOUT and interpolate)
 *         2 = one-step mode (return after each internal step taken)
 *         3 = normal tstop mode (like 0, but integration never proceeds past 
 *             TSTOP, which must be specified through the user input ROPT(8))
 *         4 = one step tstop (like 1, but integration never goes past TSTOP)
 * IER   = completion flag: 0 = success, values -1 ... -8 are various
 *         failure modes (see CVODE manual).
 * The current values of the optional outputs are available in IOPT and ROPT.
 * 
 * (8) Computing solution derivatives: FCVDKY
 * To obtain a derivative of the solution, of order up to the current method
 * order, make the following call:
 *       CALL FCVDKY (T, K, DKY, IER)
 * The arguments are:
 * T   = value of t at which solution derivative is desired
 * K   = derivative order (0 .le. K .le. QU)
 * DKY = array containing computed K-th derivative of y on return
 * IER = return flag: = 0 for success, < 0 for illegal argument.
 * 
 * (9) Memory freeing: FCVFREE and FNVFREES / FNVFREEP
 * To the free the internal memory created by the calls to FCVMALLOC and
 * FNVINITS or FNVINITP, depending on the version (serial/parallel), make
 * the following calls, in this order:
 *       CALL FCVFREE
 *       CALL FNVFREES or CALL FNVFREEP  
 * 
 * ===============================================================================
 */

/* Definitions of interface function names */

#if defined(SUNDIALS_UNDERSCORE_NONE)

#define FCV_MALLOC       fcvmalloc
#define FCV_REINIT       fcvreinit
#define FCV_DIAG         fcvdiag
#define FCV_DENSE        fcvdense
#define FCV_DENSESETJAC  fcvdensesetjac
#define FCV_BAND         fcvband
#define FCV_BANDSETJAC   fcvbandsetjac
#define FCV_SPGMR        fcvspgmr
#define FCV_SPGMRREINIT  fcvspgmrreinit
#define FCV_SPGMRSETJAC  fcvspgmrsetjac
#define FCV_SPGMRSETPSOL fcvspgmrsetpsol
#define FCV_SPGMRSETPSET fcvspgmrsetpset
#define FCV_CVODE        fcvode
#define FCV_DKY          fcvdky
#define FCV_FREE         fcvfree
#define FCV_FUN          fcvfun
#define FCV_DJAC         fcvdjac
#define FCV_BJAC         fcvbjac
#define FCV_PSOL         fcvpsol
#define FCV_PSET         fcvpset
#define FCV_JTIMES       fcvjtimes

#elif defined(SUNDIALS_UNDERSCORE_TWO)

#define FCV_MALLOC       fcvmalloc__
#define FCV_REINIT       fcvreinit__
#define FCV_DIAG         fcvdiag__
#define FCV_DENSE        fcvdense__
#define FCV_DENSESETJAC  fcvdensesetjac__
#define FCV_BAND         fcvband__
#define FCV_BANDSETJAC   fcvbandsetjac__
#define FCV_SPGMR        fcvspgmr__
#define FCV_SPGMRREINIT  fcvspgmrreinit__
#define FCV_SPGMRSETJAC  fcvspgmrsetjac__
#define FCV_SPGMRSETPSOL fcvspgmrsetpsol__
#define FCV_SPGMRSETPSET fcvspgmrsetpset__
#define FCV_CVODE        fcvode__
#define FCV_DKY          fcvdky__
#define FCV_FREE         fcvfree__
#define FCV_FUN          fcvfun__
#define FCV_DJAC         fcvdjac__
#define FCV_BJAC         fcvbjac__
#define FCV_PSOL         fcvpsol__
#define FCV_PSET         fcvpset__
#define FCV_JTIMES       fcvjtimes__

#else

#define FCV_MALLOC       fcvmalloc_
#define FCV_REINIT       fcvreinit_
#define FCV_DIAG         fcvdiag_
#define FCV_DENSE        fcvdense_
#define FCV_DENSESETJAC  fcvdensesetjac_
#define FCV_BAND         fcvband_
#define FCV_BANDSETJAC   fcvbandsetjac_
#define FCV_SPGMR        fcvspgmr_
#define FCV_SPGMRREINIT  fcvspgmrreinit_
#define FCV_SPGMRSETJAC  fcvspgmrsetjac_
#define FCV_SPGMRSETPSOL fcvspgmrsetpsol_
#define FCV_SPGMRSETPSET fcvspgmrsetpset_
#define FCV_CVODE        fcvode_
#define FCV_DKY          fcvdky_
#define FCV_FREE         fcvfree_
#define FCV_FUN          fcvfun_
#define FCV_DJAC         fcvdjac_
#define FCV_BJAC         fcvbjac_
#define FCV_PSOL         fcvpsol_
#define FCV_PSET         fcvpset_
#define FCV_JTIMES       fcvjtimes_

#endif


/* CVODE header files  */

#include "sundialstypes.h" /* definitions of type realtype         */
#include "cvode.h"         /* definition of type RHSFn             */
#include "nvector.h"       /* definition of type N_Vector, N_VSpec */
#include "dense.h"         /* definition of DenseMat               */
#include "band.h"          /* definition of BandMat                */

/* Prototypes: Functions Called by the CVODE Solver */
  
void FCVf(realtype t, N_Vector y, N_Vector ydot, void *f_data);

void FCVDenseJac(long int N, DenseMat J, realtype t, 
                 N_Vector y, N_Vector fy, void *jac_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

void FCVBandJac(long int N, long int mupper, long int mlower,
                BandMat J, realtype t, N_Vector y, N_Vector fy,
                void *jac_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FCVPSet(realtype tn, N_Vector y,N_Vector fy, booleantype jok,
            booleantype *jcurPtr, realtype gamma, void *P_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FCVPSol(realtype tn, N_Vector y, N_Vector fy, 
            N_Vector r, N_Vector z,
            realtype gamma, realtype delta,
            int lr, void *P_data, N_Vector vtemp);

int FCVJtimes(N_Vector v, N_Vector Jv, realtype t, 
              N_Vector y, N_Vector fy,
              void *jac_data, N_Vector work);


/* Declarations for global variables, shared among various routines */

extern N_Vector F2C_vec;

void *CV_cvodemem;
N_Vector CV_yvec, CV_atolvec;
booleantype CV_optin;
long int *CV_iopt;
realtype *CV_ropt;
int CV_nrtfn;
int CV_ls;    /* 1 = DENSE, 2 = BAND, 3 = DIAG, 4 = SPGMR */

#endif
