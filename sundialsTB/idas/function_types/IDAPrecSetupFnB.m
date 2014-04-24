%IDAPrecSetupFnB - type for preconditioner setup function for backward problems.
%
%   The function PSETFUNB must be defined either as
%        FUNCTION FLAG = PSETFUNB(T,YY,YP,YYB,YPB,RRB,CJB)
%   or as
%        FUNCTION [FLAG,NEW_DATA] = PSETFUNB(T,YY,YP,YYB,YPB,RRB,CJB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDASetUserData.
%
%   See also IDAPrecSolveFnB, IDAPrecSetupFn, IDASetOptions
%
%   NOTE: PSETFUN and PSETFUNB are specified through the property
%   PrecSetupFn to IDASetOptions and are used only if the property
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
