%IDAQuadRhsFn - type for user provided quadrature RHS function.
%
%IVP Problem
%
%   The function QFUN must be defined as 
%        FUNCTION [YQD, FLAG] = QFUN(T, YY, YP)
%   and must return a vector YQD corresponding to fQ(t,yy,yp), the 
%   integrand for the integral to be evaluated.
%   If a user data structure DATA was specified in IDAMalloc, then
%   QFUN must be defined as
%        FUNCTION [YQD, FLAG, NEW_DATA] = QFUN(T, YY, YP, DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector YQD,
%   the QFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function QFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%Adjoint Problem
%
%   The function QFUNB must be defined either as
%        FUNCTION [YQBD, FLAG] = QFUNB(T, YY, YP, YYB, YPB)
%   or as
%        FUNCTION [YQBD, FLAG, NEW_DATA] = QFUNB(T, YY, YP, YYB, YPB, DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAMalloc. In either case, it must return the vector YQBD
%   corresponding to fQB(t,yy,yp,yyB,ypB), the integrand for the integral to be 
%   evaluated on the backward phase.
%
%   The function QFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDASetOptions
%
%   NOTE: QFUN and QFUNB are specified through the property
%   QuadRhsFn to IDASetOptions and are used only if the property
%   Quadratures was set to 'on'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $
