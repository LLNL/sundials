%CVPrecSetupFn - type for user provided preconditioner setup function.
%
%   The user-supplied preconditioner setup function PSETFUN and
%   the user-supplied preconditioner solve function PSOLFUN
%   together must define left and right preconditoner matrices
%   P1 and P2 (either of which may be trivial), such that the
%   product P1*P2 is an approximation to the Newton matrix
%   M = I - gamma*J.  Here J is the system Jacobian J = df/dy,
%   and gamma is a scalar proportional to the integration step
%   size h.  The solution of systems P z = r, with P = P1 or P2,
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
%   ODEFUN with the same (t,y) arguments.  Thus the PSETFUN
%   function can use any auxiliary data that is computed and
%   saved by the ODEFUN function and made accessible to PSETFUN.
%
%   The function PSETFUN must be defined as 
%        FUNCTION [JCUR, FLAG] = PSETFUN(T,Y,FY,JOK,GAMMA)
%   and must return a logical flag JCUR (true if Jacobian information
%   was recomputed and false if saved data was reused). If PSETFUN
%   was successful, it must return FLAG=0. For a recoverable error (in    
%   which case the setup will be retried) it must set FLAG to a positive
%   integer value. If an unrecoverable error occurs, it must set FLAG
%   to a negative value, in which case the integration will be halted.
%   The input argument FY contains the current value of f(t,y).
%   If the input logical flag JOK is false, it means that
%   Jacobian-related data must be recomputed from scratch. If it is true,
%   it means that Jacobian data, if saved from the previous PSETFUN call
%   can be reused (with the current value of GAMMA).
%
%   If a user data structure DATA was specified in CVodeInit, then
%   PSETFUN must be defined as
%        FUNCTION [JCUR, FLAG, NEW_DATA] = PSETFUN(T,Y,FY,JOK,GAMMA,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the flags JCUR
%   and FLAG, the PSETFUN function must also set NEW_DATA. Otherwise, it 
%   should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead
%   to unnecessary copying).
%
%   See also CVPrecSolveFn, CVodeSetOptions
%
%   NOTE: PSETFUN is specified through the property PrecSetupFn to
%   CVodeSetOptions and is used only if the property LinearSolver was 
%   set to 'GMRES', 'BiCGStab', or 'TFQMR' and if the property PrecType
%   is not 'None'.

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
% $Revision$Date: 2007/05/11 18:51:33 $
