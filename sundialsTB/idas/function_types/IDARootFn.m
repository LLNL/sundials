%IDARootFn - type for user provided root-finding function.
%
%   The function ROOTFUN must be defined as 
%        FUNCTION [G, FLAG] = ROOTFUN(T,YY,YP)
%   and must return a vector G corresponding to g(t,yy,yp).
%   If a user data structure DATA was specified in IDAMalloc, then
%   ROOTFUN must be defined as
%        FUNCTION [G, FLAG, NEW_DATA] = ROOTFUN(T,YY,YP,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector G,
%   the ROOTFUN function must also set NEW_DATA. Otherwise, it should 
%   set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to 
%   unnecessary copying).
%
%   The function ROOTFUN must set FLAG=0 if successful, or FLAG~=0 if
%   a failure occurred.
%
%   See also IDASetOptions
%
%   NOTE: ROOTFUN is specified through the RootsFn property in 
%   IDASetOptions and is used only if the property NumRoots is a
%   positive integer.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $
