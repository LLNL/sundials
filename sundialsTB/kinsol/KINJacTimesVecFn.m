%KINJacTimesVecFn - type for user provided Jacobian times vector function.
%
%   The function JTVFUN must be defined as 
%        FUNCTION [JV, FLAG, IER] = JTVFUN(Y,V,FLAG)
%   and must return a vector JV corresponding to the product of the 
%   Jacobian of f(y) with the vector v. On input, FLAG indicates if
%   the iterate has been updated in the interim. JV must be update
%   or reevaluated, if appropriate, unless FLAG=false. This flag must
%   be reset by the user.
%   If successful, IER should be set to 0. If an error occurs, IER should
%   be set to a nonzero value.
%   If a user data structure DATA was specified in KINMalloc, then
%   JTVFUN must be defined as
%        FUNCTION [JV, FLAG, IER, NEW_DATA] = JTVFUN(Y,V,FLAG,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector JV, and
%   flags FLAG and IER, the JTVFUN function must also set NEW_DATA. Otherwise, 
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   See also KINSetOptions
%
%   NOTE: JTVFUN is specified through the property JacobianFn to KINSetOptions 
%   and is used only if the property LinearSolver was set to 'GMRES' or 'BiCGStab'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
