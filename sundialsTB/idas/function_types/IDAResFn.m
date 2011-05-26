%IDAResFn - type for residual function
%
%   The function DAEFUN must be defined as 
%        FUNCTION [R, FLAG] = DAEFUN(T, YY, YP)
%   and must return a vector R corresponding to f(t,yy,yp).
%   If a user data structure DATA was specified in IDAInit, then
%   DAEFUN must be defined as
%        FUNCTION [R, FLAG, NEW_DATA] = DAEFUN(T, YY, YP, DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector YD,
%   the DAEFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function DAEFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/08/21 17:38:44 $
