%CVJacTimesVecFnB - type for user provided Jacobian times vector function for backward problems.
%
%   The function JTVFUNB must be defined either as
%        FUNCTION [JVB, FLAG] = JTVFUNB(T,Y,YB,FYB,VB)
%   or as
%        FUNCTION [JVB, FLAG, NEW_DATA] = JTVFUNB(T,Y,YB,FYB,VB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the vector JVB, the
%   product of the Jacobian of fB(t,y,yB) with respect to yB and a vector
%   vB. The input argument FYB contains the current value of f(t,y,yB).
%
%   The function JTVFUNB must set FLAG=0 if successful, or FLAG~=0 if
%   a failure occurred.
%
%   See also CVodeSetOptions
%
%   NOTE: JTVFUNB is specified through the property JacobianFn to
%   CVodeSetOptions and is used only if the property LinearSolver
%   was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

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
