%IDAPrecSetupFn - type for preconditioner setup function.
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
%   Each call to the PSETFUN function is preceded by a call to
%   DAEFUN with the same (t,yy,yp) arguments.  Thus the PSETFUN
%   function can use any auxiliary data that is computed and
%   saved by the DAEFUN function and made accessible to PSETFUN.
%
%   The function PSETFUN must be defined as 
%        FUNCTION FLAG = PSETFUN(T,YY,YP,RR,CJ)
%   If successful, it must return FLAG=0. For a recoverable error (in    
%   which case the setup will be retried) it must set FLAG to a positive
%   integer value. If an unrecoverable error occurs, it must set FLAG
%   to a negative value, in which case the integration will be halted.
%   The input argument RR contains the current value of f(t,yy,yp).
%
%   If a user data structure DATA was specified in IDASetUserData, then
%   PSETFUN must be defined as
%        FUNCTION [FLAG,NEW_DATA] = PSETFUN(T,YY,YP,RR,CJ,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the flag
%   FLAG, the PSETFUN function must also set NEW_DATA. Otherwise, it 
%   should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead
%   to unnecessary copying).
%
%   See also IDAPrecSolveFn, IDASetOptions
%
%   NOTE: PSETFUN and PSETFUNB are specified through the property
%   PrecSetupFn to IDASetOptions and are used only if the property
%   LinearSolver was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:38:44 $
