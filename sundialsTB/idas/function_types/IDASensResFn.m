%IDASensRhsFn - type for user provided sensitivity RHS function.
%
%   The function DAESFUN must be defined as 
%        FUNCTION [RS, FLAG] = DAESFUN(T,YY,YP,YYS,YPS)
%   and must return a matrix RS corresponding to fS(t,yy,yp,yyS,ypS).
%   If a user data structure DATA was specified in IDAMalloc, then
%   DAESFUN must be defined as
%        FUNCTION [RS, FLAG, NEW_DATA] = DAESFUN(T,YY,YP,YYS,YPS,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the matrix YSD,
%   the ODESFUN function must also set NEW_DATA. Otherwise, it should
%   set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to 
%   unnecessary copying).
%
%   The function DAESFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDASetFSAOptions
%
%   NOTE: DAESFUN is specified through the property FSAResFn to 
%         IDASetFSAOptions.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $
