%CVPrecSetupFnB - type for user provided preconditioner setup function for backward problems.
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
%
%   The function PSETFUNB must be defined either as
%        FUNCTION [JCURB, FLAG] = PSETFUNB(T,Y,YB,FYB,JOK,GAMMAB)
%   or as
%        FUNCTION [JCURB, FLAG, NEW_DATA] = PSETFUNB(T,Y,YB,FYB,JOK,GAMMAB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the flags JCURB and FLAG.
%
%   See also CVPrecSolveFnB, CVodeSetOptions
%
%   NOTE: PSETFUNB is specified through the property PrecSetupFn to
%   CVodeSetOptions and is used only if the property LinearSolver was 
%   set to 'GMRES', 'BiCGStab', or 'TFQMR' and if the property PrecType
%   is not 'None'.

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2011/06/01 20:44:05 $
