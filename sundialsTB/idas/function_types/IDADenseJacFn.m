%IDADenseJacFn - type for dense Jacobian function.
%
%   The function DJACFUN must be defined as 
%        FUNCTION [J, FLAG] = DJACFUN(T, YY, YP, RR, CJ)
%   and must return a matrix J corresponding to the Jacobian 
%   (df/dyy + cj*df/dyp).
%   The input argument RR contains the current value of f(t,yy,yp).
%   If a user data structure DATA was specified in IDAInit, then
%   DJACFUN must be defined as
%        FUNCTION [J, FLAG, NEW_DATA] = DJACFUN(T, YY, YP, RR, CJ, DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the matrix J,
%   the DJACFUN function must also set NEW_DATA. Otherwise, it should
%   set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to 
%   unnecessary copying).
%
%   The function DJACFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDASetOptions
%
%   NOTE: DJACFUN is specified through the property JacobianFn to 
%   IDASetOptions and is used only if the property LinearSolver 
%   was set to 'Dense'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/08/21 17:38:44 $
