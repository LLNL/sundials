%IDAResFn - type for user provided RHS type
%
%IVP Problem
%
%   The function DAEFUN must be defined as 
%        FUNCTION [R, FLAG] = DAEFUN(T, YY, YP)
%   and must return a vector R corresponding to f(t,yy,yp).
%   If a user data structure DATA was specified in IDAMalloc, then
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
%Adjoint Problem
%
%   The function DAEFUNB must be defined either as
%        FUNCTION [RB, FLAG] = DAEFUNB(T, YY, YP, YYB, YPB)
%   or as
%        FUNCTION [RB, FLAG, NEW_DATA] = DAEFUNB(T, YY, YP, YYB, YPB, DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAMalloc. In either case, it must return the vector RB
%   corresponding to fB(t,yy,yp,yyB,ypB).
%
%   The function DAEFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAMalloc, IDAMallocB
%
%   NOTE: DAEFUN and DAEFUNB are specified through the IDAMalloc and
%   IDAMallocB functions, respectively.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $
