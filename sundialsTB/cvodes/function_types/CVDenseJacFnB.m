%CVDenseJacFnB - type for user provided dense Jacobian function for backward problems.
%
%   The function DJACFUNB must be defined either as
%        FUNCTION [JB, FLAG] = DJACFUNB(T, Y, YB, FYB)
%   or as
%        FUNCTION [JB, FLAG, NEW_DATA] = DJACFUNB(T, Y, YB, FYB, DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the matrix JB, the
%   Jacobian of fB(t,y,yB), with respect to yB. The input argument
%   FYB contains the current value of f(t,y,yB).
%
%   The function DJACFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeSetOptions
%
%   NOTE: DJACFUNB is specified through the property JacobianFn to
%   CVodeSetOptions and is used only if the property LinearSolver
%   was set to 'Dense'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/15 19:31:24 $
