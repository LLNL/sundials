%KINPrecSolveFn - type for user provided preconditioner solve function.
%
%   The user-supplied preconditioner solve function PSOLFN
%   is to solve a linear system P z = r in which the matrix P is
%   the preconditioner matrix (possibly set implicitely by PSETFUN)
%
%   The function PSOLFUN must be defined as 
%        FUNCTION [Z,IER] = PSOLFUN(Y,YSCALE,FY,FSCALE,R)
%   and must return a vector Z containing the solution of Pz=r.
%   If successful, PSOLFUN must return IER=0. If an error occurs, then IER
%   must be set to a non-zero value.
%   The input argument FY contains the current value of f(y), while YSCALE
%   and FSCALE are the scaling vectors for solution and system function,
%   respectively (as passed to KINSol)
%
%   If a user data structure DATA was specified in KINMalloc, then
%   PSOLFUN must be defined as
%        FUNCTION [Z, IER, NEW_DATA] = PSOLFUN(Y,YSCALE,FY,FSCALE,R,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector Z and
%   the flag IER, the PSOLFUN function must also set NEW_DATA. Otherwise,
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would
%   lead to unnecessary copying).
%
%   See also KINPrecSetupFn, KINSetOptions
%
%   NOTE: PSOLFUN is specified through the property PrecSolveFn to KINSetOptions 
%   and is used only if the property LinearSolver was set to 'GMRES' or 'BiCGStab'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
