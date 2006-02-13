%CVBandJacFn - type for user provided banded Jacobian function.
%
%IVP Problem
%
%   The function BJACFUN must be defined as 
%        FUNCTION [J, FLAG] = BJACFUN(T,Y,FY)
%   and must return a matrix J corresponding to the banded Jacobian of f(t,y).
%   The input argument FY contains the current value of f(t,y).
%   If a user data structure DATA was specified in CVodeMalloc, then
%   BJACFUN must be defined as
%        FUNCTION [J, FLAG, NEW_DATA] = BJACFUN(T,Y,FY,DATA)
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
%Adjoint Problem
%
%   The function BJACFUNB must be defined either as
%        FUNCTION [JB, FLAG] = BJACFUNB(T,Y,YB,FYB)
%   or as
%        FUNCTION [JB, FLAG, NEW_DATA] = BJACFUNB(T,Y,YB,FYB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the matrix JB, the
%   Jacobian of fB(t,y,yB), with respect to yB. The input argument
%   FYB contains the current value of f(t,y,yB).
%
%   The function BJACFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeSetOptions
%
%   See the CVODES user guide for more informaiton on the structure of
%   a banded Jacobian.
%
%   NOTE: BJACFUN and BJACFUNB are specified through the property
%   JacobianFn to CVodeSetOptions and are used only if the property
%   LinearSolver was set to 'Band'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:41 $
