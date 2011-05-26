%IDAQuadRhsFnB - type for quadrature RHS function for backward problems
%
%   The function QFUNB must be defined either as
%        FUNCTION [YQBD, FLAG] = QFUNB(T, YY, YP, YYB, YPB)
%   or as
%        FUNCTION [YQBD, FLAG, NEW_DATA] = QFUNB(T, YY, YP, YYB, YPB, DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. In either case, it must return the vector YQBD
%   corresponding to fQB(t,yy,yp,yyB,ypB), the integrand for the integral to be 
%   evaluated on the backward phase.
%
%   The function QFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAQuadInitB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/08/21 17:38:44 $
