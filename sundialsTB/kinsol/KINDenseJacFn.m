%KINDenseJacFn - type for user provided dense Jacobian function.
%
%   The function DJACFUN must be defined as 
%        FUNCTION [J, IER] = DJACFUN(Y,FY)
%   and must return a matrix J corresponding to the Jacobian of f(y).
%   The input argument FY contains the current value of f(y).
%   If successful, IER should be set to 0. If an error occurs, IER should
%   be set to a nonzero value.
%   If a user data structure DATA was specified in KINMalloc, then
%   DJACFUN must be defined as
%        FUNCTION [J, IER, NEW_DATA] = DJACFUN(Y,FY,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the matrix J and
%   the flag IER, the DJACFUN function must also set NEW_DATA. Otherwise, 
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead
%   to unnecessary copying).
%
%   See also KINSetOptions
%
%   NOTE: DJACFUN is specified through the property JacobianFn to KINSetOptions 
%   and is used only if the property LinearSolver was set to 'Dense'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
