%CVPrecSolveFnB - type for user provided preconditioner solve function for backward problems.
%
%   The user-supplied preconditioner solve function PSOLFN
%   is to solve a linear system P z = r in which the matrix P is
%   one of the preconditioner matrices P1 or P2, depending on the
%   type of preconditioning chosen.
%
%   The function PSOLFUNB must be defined either as
%        FUNCTION [ZB, FLAG] = PSOLFUNB(T,Y,YB,FYB,RB)
%   or as
%        FUNCTION [ZB, FLAG, NEW_DATA] = PSOLFUNB(T,Y,YB,FYB,RB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the vector ZB and the
%   flag FLAG.
%
%   See also CVPrecSetupFnB, CVodeSetOptions
%
%   NOTE: PSOLFUNB is specified through the property PrecSolveFn to 
%   CVodeSetOptions and is used only if the property LinearSolver was
%   set to 'GMRES', 'BiCGStab', or 'TFQMR' and if the property PrecType
%   is not 'None'.

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
% $Revision$Date: 2007/05/11 18:51:33 $
