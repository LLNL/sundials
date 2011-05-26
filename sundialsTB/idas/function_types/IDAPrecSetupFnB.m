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
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/08/21 17:38:44 $
