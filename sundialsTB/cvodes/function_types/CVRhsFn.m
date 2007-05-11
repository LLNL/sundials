%CVRhsFn - type for user provided RHS function
%
%   The function ODEFUN must be defined as 
%        FUNCTION [YD, FLAG] = ODEFUN(T,Y)
%   and must return a vector YD corresponding to f(t,y).
%   If a user data structure DATA was specified in CVodeMalloc, then
%   ODEFUN must be defined as
%        FUNCTION [YD, FLAG, NEW_DATA] = ODEFUN(T,Y,DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector YD,
%   the ODEFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function ODEFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeInit
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/02/13 23:01:20 $
