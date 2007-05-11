%IDAPrecSetupFn - type for user provided preconditioner setup function.
%
%   The user-supplied preconditioner setup function PSETFUN and
%   the user-supplied preconditioner solve function PSOLFUN
%   together must define a preconditoner matrix P which is an 
%   approximation to the Newton matrix M = J_yy - cj*J_yp.  
%   Here J_yy = df/dyy, J_yp = df/dyp, and cj is a scalar proportional 
%   to the integration step size h.  The solution of systems P z = r,
%   is to be carried out by the PrecSolve function, and PSETFUN
%   is to do any necessary setup operations.
%
%   The user-supplied preconditioner setup function PSETFUN
%   is to evaluate and preprocess any Jacobian-related data
%   needed by the preconditioner solve function PSOLFUN.
%   This might include forming a crude approximate Jacobian,
%   and performing an LU factorization on the resulting
%   approximation to M.  This function will not be called in
%   advance of every call to PSOLFUN, but instead will be called
%   only as often as necessary to achieve convergence within the
%   Newton iteration.  If the PSOLFUN function needs no
%   preparation, the PSETFUN function need not be provided.
%
%   For greater efficiency, the PSETFUN function may save
%   Jacobian-related data and reuse it, rather than generating it
%   from scratch.  In this case, it should use the input flag JOK
%   to decide whether to recompute the data, and set the output
%   flag JCUR accordingly.
%
%   Each call to the PSETFUN function is preceded by a call to
%   DAEFUN with the same (t,yy,yp) arguments.  Thus the PSETFUN
%   function can use any auxiliary data that is computed and
%   saved by the DAEFUN function and made accessible to PSETFUN.
%
%IVP Problem
%
%   The function PSETFUN must be defined as 
%        FUNCTION FLAG = PSETFUN(T,YY,YP,RR,CJ)
%   If successful, it must return FLAG=0. For a recoverable error (in    
%   which case the setup will be retried) it must set FLAG to a positive
%   integer value. If an unrecoverable error occurs, it must set FLAG
%   to a negative value, in which case the integration will be halted.
%   The input argument RR contains the current value of f(t,yy,yp).
%
%   If a user data structure DATA was specified in IDAMalloc, then
%   PSETFUN must be defined as
%        FUNCTION [FLAG,NEW_DATA] = PSETFUN(T,YY,YP,RR,CJ,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the flags JCUR
%   and FLAG, the PSETFUN function must also set NEW_DATA. Otherwise, it 
%   should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead
%   to unnecessary copying).
%
%Adjoint Problem
%
%   The function PSETFUNB must be defined either as
%        FUNCTION FLAG = PSETFUNB(T,YY,YP,YYB,YPB,RRB,CJB)
%   or as
%        FUNCTION [FLAG,NEW_DATA] = PSETFUNB(T,YY,YP,YYB,YPB,RRB,CJB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAMalloc.
%
%   See also IDAPrecSolveFn, IDASetOptions
%
%   NOTE: PSETFUN and PSETFUNB are specified through the property
%   PrecSetupFn to IDASetOptions and are used only if the property
%   LinearSolver was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $
