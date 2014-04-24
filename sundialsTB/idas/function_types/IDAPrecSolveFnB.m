%IDAPrecSolveFnB - type for preconditioner solve function.
%
%   The user-supplied preconditioner solve function PSOLFUNB
%   is to solve a linear system P z = r, where P is the
%   preconditioner matrix.
%
%   The function PSOLFUNB must be defined either as
%        FUNCTION [ZB,FLAG] = PSOLFUNB(T,YY,YP,YYB,YPB,RRB,RB)
%   or as
%        FUNCTION [ZB,FLAG,NEW_DATA] = PSOLFUNB(T,YY,YP,YYB,YPB,RRB,RB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. In either case, it must return the vector ZB and the
%   flag FLAG.
%
%   See also IDAPrecSetupFnB, IDAPrecSolveFn, IDASetOptions
%
%   NOTE: PSOLFUN and PSOLFUNB are specified through the property
%   PrecSolveFn to IDASetOptions and are used only if the property
%   LinearSolver was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

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
