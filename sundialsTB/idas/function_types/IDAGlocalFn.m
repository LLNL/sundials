%IDAGlocalFn - type for RES approximation function (BBDPre).
%
%   The function GLOCFUN must be defined as 
%        FUNCTION [GLOC, FLAG] = GLOCFUN(T,YY,YP)
%   and must return a vector GLOC corresponding to an approximation to f(t,yy,yp)
%   which will be used in the BBDPRE preconditioner module. The case where
%   G is mathematically identical to F is allowed.
%   If a user data structure DATA was specified in IDAInit, then
%   GLOCFUN must be defined as
%        FUNCTION [GLOC, FLAG, NEW_DATA] = GLOCFUN(T,YY,YP,DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector G,
%   the GLOCFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function GLOCFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAGcommFn, IDASetOptions
%
%   NOTE: GLOCFUN and GLOCFUNB are specified through the GlocalFn property
%   in IDASetOptions and are used only if the property PrecModule
%   is set to 'BBDPre'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/08/21 17:38:44 $
