%KINGlocalFn - type for user provided RHS approximation function (BBDPre).
%
%   The function GLOCFUN must be defined as 
%        FUNCTION G = GLOCFUN(Y)
%   and must return a vector G corresponding to an approximation to f(y)
%   which will be used in the BBDPRE preconditioner module. The case where
%   G is mathematically identical to F is allowed.
%   If a user data structure DATA was specified in KINMalloc, then
%   GLOCFUN must be defined as
%        FUNCTION [G, NEW_DATA] = GLOCFUN(Y,DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector G,
%   the GLOCFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   See also KINGcommFn, KINSetOptions
%
%   NOTE: GLOCFUN is specified through the GlocalFn property in KINSetOptions 
%   and is used only if the property PrecModule is set to 'BBDPre'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
