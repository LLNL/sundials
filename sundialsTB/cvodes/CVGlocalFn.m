%CVGlocalFn - type for user provided RHS approximation function (BBDPre).
%
%IVP Problem
%
%   The function GLOCFUN must be defined as 
%        FUNCTION G = GLOCFUN(T,Y)
%   and must return a vector G corresponding to an approximation to f(t,y)
%   which will be used in the BBDPRE preconditioner module. The case where
%   G is mathematically identical to F is allowed.
%   If a user data structure DATA was specified in CVodeMalloc, then
%   GLOCFUN must be defined as
%        FUNCTION [G, NEW_DATA] = GLOCFUN(T,Y,DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector G,
%   the GLOCFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%Adjoint Problem
%
%   The function GLOCFUNB must be defined either as
%        FUNCTION GB = GLOCFUNB(T,Y,YB)
%   or as
%        FUNCTION [GB, NEW_DATA] = GLOCFUNB(T,Y,YB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the vector GB
%   corresponding to an approximation to fB(t,y,yB).
%
%   See also CVGcommFn, CVodeSetOptions
%
%   NOTE: GLOCFUN and GLOCFUNB are specified through the GlocalFn property
%   in CVodeSetOptions and are used only if the property PrecModule
%   is set to 'BBDPre'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
