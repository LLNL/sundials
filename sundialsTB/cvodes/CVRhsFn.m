%CVRhsFn - type for user provided RHS type
%
%IVP Problem
%
%   The function ODEFUN must be defined as 
%        FUNCTION YD = ODEFUN(T,Y)
%   and must return a vector YD corresponding to f(t,y).
%   If a user data structure DATA was specified in CVodeMalloc, then
%   ODEFUN must be defined as
%        FUNCTION [YD, NEW_DATA] = ODEFUN(T,Y,DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector YD,
%   the ODEFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%Adjoint Problem
%
%   The function ODEFUNB must be defined either as
%        FUNCTION YBD = ODEFUNB(T,Y,YB)
%   or as
%        FUNCTION [YBD, NEW_DATA] = ODEFUNB(T,Y,YB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the vector YBD
%   corresponding to fB(t,y,yB).
%
%   See also CVodeMalloc, CVodeMallocB
%
%   NOTE: ODEFUN and ODEFUNB are specified through the CVodeMalloc and
%   CVodeMallocB functions, respectively.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
