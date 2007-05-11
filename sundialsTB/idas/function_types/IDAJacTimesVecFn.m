%IDAJacTimesVecFn - type for user provided Jacobian times vector function.
%
%IVP Problem
%
%   The function JTVFUN must be defined as 
%        FUNCTION [JV, FLAG] = JTVFUN(T,YY,YP,RR,V,CJ)
%   and must return a vector JV corresponding to the product of the 
%   Jacobian ( df/dyy + cj * df/dyp ) with the vector v.
%   The input argument RR contains the current value of f(t,yy,yp).
%   If a user data structure DATA was specified in IDAMalloc, then
%   JTVFUN must be defined as
%        FUNCTION [JV, FLAG, NEW_DATA] = JTVFUN(T,YY,YP,RR,V,CJ,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector JV,
%   the JTVFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function JTVFUN must set FLAG=0 if successful, or FLAG~=0 if
%   a failure occurred.
%
%Adjoint Problem
%
%   The function JTVFUNB must be defined either as
%        FUNCTION [JVB,FLAG] = JTVFUNB(T,YY,YP,YYB,YPB,RRB,VB,CJB)
%   or as
%        FUNCTION [JVB,FLAG,NEW_DATA] = JTVFUNB(T,YY,YP,YYB,YPB,RRB,VB,CJB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAMalloc. In either case, it must return the vector JVB, the
%   product of the Jacobian (dfB/dyyB + cj * dfB/dypB) and a vector
%   vB. The input argument RRB contains the current value of f(t,yy,yp,yyB,ypB).
%
%   The function JTVFUNB must set FLAG=0 if successful, or FLAG~=0 if
%   a failure occurred.
%
%   See also IDASetOptions
%
%   NOTE: JTVFUN and JTVFUNB are specified through the property
%   JacobianFn to IDASetOptions and are used only if the property
%   LinearSolver was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $
