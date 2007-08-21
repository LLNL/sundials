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
%   IDAMalloc. In either case, it must return the vector ZB and the
%   flag FLAG.
%
%   See also IDAPrecSetupFnB, IDAPrecSolveFn, IDASetOptions
%
%   NOTE: PSOLFUN and PSOLFUNB are specified through the property
%   PrecSolveFn to IDASetOptions and are used only if the property
%   LinearSolver was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2007/05/11 18:48:45 $
