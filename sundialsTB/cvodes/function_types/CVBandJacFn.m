%CVBandJacFn - type for user provided banded Jacobian function.
%
%   The function BJACFUN must be defined as 
%        FUNCTION [J, FLAG] = BJACFUN(T, Y, FY)
%   and must return a matrix J corresponding to the banded Jacobian of f(t,y).
%   The input argument FY contains the current value of f(t,y).
%   If a user data structure DATA was specified in CVodeMalloc, then
%   BJACFUN must be defined as
%        FUNCTION [J, FLAG, NEW_DATA] = BJACFUN(T, Y, FY, DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the matrix J,
%   the BJACFUN function must also set NEW_DATA. Otherwise, it should 
%   set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to 
%   unnecessary copying).
%
%   The function BJACFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeSetOptions
%
%   See the CVODES user guide for more informaiton on the structure of
%   a banded Jacobian.
%
%   NOTE: BJACFUN is specified through the property JacobianFn to 
%   CVodeSetOptions and is used only if the property LinearSolver
%   was set to 'Band'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/15 19:31:24 $
