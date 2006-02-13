%CVQuadRhsFn - type for user provided quadrature RHS function.
%
%IVP Problem
%
%   The function ODEQFUN must be defined as 
%        FUNCTION [YQD, FLAG] = ODEQFUN(T,Y)
%   and must return a vector YQD corresponding to fQ(t,y), the integrand
%   for the integral to be evaluated.
%   If a user data structure DATA was specified in CVodeMalloc, then
%   ODEQFUN must be defined as
%        FUNCTION [YQD, FLAG, NEW_DATA] = ODEQFUN(T,Y,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector YQD,
%   the ODEQFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function ODEQFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%Adjoint Problem
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
%   See also CVodeSetOptions
%
%   NOTE: ODEQFUN and ODEQFUNB are specified through the property
%   QuadRhsFn to CVodeSetOptions and are used only if the property
%   Quadratures was set to 'on'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:41 $
