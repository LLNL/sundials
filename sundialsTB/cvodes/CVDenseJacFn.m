%CVDenseJacFn - type for user provided dense Jacobian function.
%
%IVP Problem
%
%   The function DJACFUN must be defined as 
%        FUNCTION J = DJACFUN(T,Y,FY)
%   and must return a matrix J corresponding to the Jacobian of f(t,y).
%   The input argument FY contains the current value of f(t,y).
%   If a user data structure DATA was specified in CVodeMalloc, then
%   DJACFUN must be defined as
%        FUNCTION [J, NEW_DATA] = DJACFUN(T,Y,FY,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the matrix J,
%   the DJACFUN function must also set NEW_DATA. Otherwise, it should
%   set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to 
%   unnecessary copying).
%
%Adjoint Problem
%
%   The function DJACFUNB must be defined either as
%        FUNCTION JB = DJACFUNB(T,Y,YB,FYB)
%   or as
%        FUNCTION [JB, NEW_DATA] = DJACFUNB(T,Y,YB,FYB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the matrix JB, the
%   Jacobian of fB(t,y,yB), with respect to yB. The input argument
%   FYB contains the current value of f(t,y,yB).
%
%   See also CVodeSetOptions
%
%   NOTE: DJACFUN and DJACFUNB are specified through the property
%   JacobianFn to CVodeSetOptions and are used only if the property
%   LinearSolver was set to 'Dense'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
