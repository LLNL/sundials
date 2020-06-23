%KINPrecSolveFn - type for user provided preconditioner solve function.
%
%   The user-supplied preconditioner solve function PSOLFN
%   is to solve a linear system P z = r in which the matrix P is
%   the preconditioner matrix (possibly set implicitely by PSETFUN)
%
%   The function PSOLFUN must be defined as 
%        FUNCTION [Z, FLAG] = PSOLFUN(Y, YSCALE, FY, FSCALE, R)
%   and must return a vector Z containing the solution of Pz=r.
%   The input argument FY contains the current value of f(y), while YSCALE
%   and FSCALE are the scaling vectors for solution and system function,
%   respectively (as passed to KINSol)
%
%   If a user data structure DATA was specified in KINInit, then
%   PSOLFUN must be defined as
%        FUNCTION [Z, FLAG, NEW_DATA] = PSOLFUN(Y,YSCALE,FY,FSCALE,R,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector Z and
%   the flag FLAG, the PSOLFUN function must also set NEW_DATA. Otherwise,
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would
%   lead to unnecessary copying).
%
%   If successful, PSOLFUN must return FLAG=0. For a recoverable error it 
%   must set FLAG to a positive value (in which case the solver will attempt 
%   to correct). If an unrecoverable error occurs, it must set FLAG
%   to a negative value, in which case the solver will halt.
%
%   See also KINPrecSetupFn, KINSetOptions
%
%   NOTE: PSOLFUN is specified through the property PrecSolveFn to KINSetOptions 
%   and is used only if the property LinearSolver was set to 'GMRES' or 'BiCGStab'.

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
% $Revision$Date: 2007/05/11 18:48:46 $
