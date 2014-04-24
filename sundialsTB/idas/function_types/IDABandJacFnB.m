%IDABandJacFnB - type for banded Jacobian function for backward problems.
%
%   The function BJACFUNB must be defined either as
%        FUNCTION [JB, FLAG] = BJACFUNB(T, YY, YP, YYB, YPB, RRB, CJB)
%   or as
%        FUNCTION [JB,FLAG,NEW_DATA] = BJACFUNB(T,YY,YP,YYB,YPB,RRB,CJB)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. In either case, it must return the matrix JB, the
%   Jacobian (dfB/dyyB + cjB*dfB/dypB)of fB(t,y,yB). The input argument
%   RRB contains the current value of f(t,yy,yp,yyB,ypB).
%
%   The function BJACFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDASetOptions
%
%   See the IDAS user guide for more information on the structure of
%   a banded Jacobian.
%
%   NOTE: BJACFUNB is specified through the property JacobianFn to 
%   IDASetOptions and is used only if the property LinearSolver 
%   was set to 'Band'.

% Radu Serban <radu@llnl.gov>
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
% $Revision$Date: 2007/08/21 17:38:44 $
