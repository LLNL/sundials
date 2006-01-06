%CVJacTimesVecFn - type for user provided Jacobian times vector function.
%
%IVP Problem
%
%   The function JTVFUN must be defined as 
%        FUNCTION JV = JTVFUN(T,Y,FY,V)
%   and must return a vector JV corresponding to the product of the 
%   Jacobian of f(t,y) with the vector v.
%   The input argument FY contains the current value of f(t,y).
%   If a user data structure DATA was specified in CVodeMalloc, then
%   JTVFUN must be defined as
%        FUNCTION [JV, NEW_DATA] = JTVFUN(T,Y,FY,V,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector JV,
%   the JTVFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%Adjoint Problem
%
%   The function JTVFUNB must be defined either as
%        FUNCTION JVB = JTVFUNB(T,Y,YB,FYB,VB)
%   or as
%        FUNCTION [JVB, NEW_DATA] = JTVFUNB(T,Y,YB,FYB,VB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the vector JVB, the
%   product of the Jacobian of fB(t,y,yB) with respect to yB and a vector
%   vB. The input argument FYB contains the current value of f(t,y,yB).
%
%   See also CVodeSetOptions
%
%   NOTE: JTVFUN and JTVFUNB are specified through the property
%   JacobianFn to CVodeSetOptions and are used only if the property
%   LinearSolver was set to 'GMRES' or 'BiCGStab'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
