%CVQuadRhsFnB - type for user provided quadrature RHS function for backward problems
%
%   The function ODEQFUNB must be defined either as
%        FUNCTION [YQBD, FLAG] = ODEQFUNB(T,Y,YB)
%   or as
%        FUNCTION [YQBD, FLAG, NEW_DATA] = ODEQFUNB(T,Y,YB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the vector YQBD
%   corresponding to fQB(t,y,yB), the integrand for the integral to be 
%   evaluated on the backward phase.
%
%   The function ODEQFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeQuadInitB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/02/13 23:01:20 $
