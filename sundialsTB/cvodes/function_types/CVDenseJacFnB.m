%CVDenseJacFnB - type for user provided dense Jacobian function for backward problems.
%
%   The function DJACFUNB must be defined either as
%        FUNCTION [JB, FLAG] = DJACFUNB(T, Y, YB, FYB)
%   or as
%        FUNCTION [JB, FLAG, NEW_DATA] = DJACFUNB(T, Y, YB, FYB, DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the matrix JB, the
%   Jacobian of fB(t,y,yB), with respect to yB. The input argument
%   FYB contains the current value of f(t,y,yB).
%
%   The function DJACFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeSetOptions
%
%   NOTE: DJACFUNB is specified through the property JacobianFn to
%   CVodeSetOptions and is used only if the property LinearSolver
%   was set to 'Dense'.

% Radu Serban <radu@llnl.gov>
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.2 $Date: 2007/05/11 18:51:33 $
