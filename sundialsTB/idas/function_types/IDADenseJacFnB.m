%IDADenseJacFnb - type for dense Jacobian function for backward problems.
%
%   The function DJACFUNB must be defined either as
%        FUNCTION [JB, FLAG] = DJACFUNB(T, YY, YP, YYB, YPB, RRB, CJB)
%   or as
%        FUNCTION [JB,FLAG,NEW_DATA] = DJACFUNB(T,YY,YP,YYB,YPB,RRB,CJB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. In either case, it must return the matrix JB, the
%   Jacobian (dfB/dyyB + cjb*dfB/dypB). The input argument RRB contains 
%   the current value of f(t,yy,yp,yyB,ypB).
%
%   The function DJACFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDADenseJacFn, IDASetOptions
%
%   NOTE: DJACFUNB is specified through the property JacobianFn to 
%   IDASetOptions and is used only if the property LinearSolver was
%   set to 'Dense'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/08/21 17:38:44 $
