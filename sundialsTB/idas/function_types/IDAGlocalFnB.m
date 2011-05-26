%IDAGlocalFnB - type for RES approximation function (BBDPre) for backward problems.
%
%   The function GLOCFUNB must be defined either as
%        FUNCTION [GLOCB, FLAG] = GLOCFUNB(T,YY,YP,YYB,YPB)
%   or as
%        FUNCTION [GLOCB, FLAG, NEW_DATA] = GLOCFUNB(T,YY,YP,YYB,YPB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. In either case, it must return the vector GLOCB
%   corresponding to an approximation to fB(t,yy,yp,yyB,ypB).
%
%   The function GLOCFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAGcommFnB, IDAGlocalFn, IDASetOptions
%
%   NOTE: GLOCFUN and GLOCFUNB are specified through the GlocalFn property
%   in IDASetOptions and are used only if the property PrecModule
%   is set to 'BBDPre'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/08/21 17:38:44 $
