%CVPrecSolveFn - type for user provided preconditioner solve function.
%
%   The user-supplied preconditioner solve function PSOLFN
%   is to solve a linear system P z = r in which the matrix P is
%   one of the preconditioner matrices P1 or P2, depending on the
%   type of preconditioning chosen.
%
%IVP Problem
%
%   The function PSOLFUN must be defined as 
%        FUNCTION [Z, ERR] = PSOLFUN(T,Y,FY,R)
%   and must return a vector Z containing the solution of Pz=r.
%   If PSOLFUN was successful, it must return ERR=0. For a recoverable 
%   error (in which case the step will be retried) it must set ERR to a 
%   positive value. If an unrecoverable error occurs, it must set ERR
%   to a negative value, in which case the integration will be halted.
%   The input argument FY contains the current value of f(t,y).
%
%   If a user data structure DATA was specified in CVodeMalloc, then
%   PSOLFUN must be defined as
%        FUNCTION [Z, ERR, NEW_DATA] = PSOLFUN(T,Y,FY,R,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector Z and
%   the flag ERR, the PSOLFUN function must also set NEW_DATA. Otherwise,
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would
%   lead to unnecessary copying).
%
%Adjoint Problem
%
%   The function PSOLFUNB must be defined either as
%        FUNCTION [ZB, ERR] = PSOLFUNB(T,Y,YB,FYB,RB)
%   or as
%        FUNCTION [ZB, ERR, NEW_DATA] = PSOLFUNB(T,Y,YB,FYB,RB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeMalloc. In either case, it must return the vector ZB and the
%   flag ERR.
%
%   See also CVPrecSetupFn, CVodeSetOptions
%
%   NOTE: PSOLFUN and PSOLFUNB are specified through the property
%   PrecSolveFn to CVodeSetOptions and are used only if the property
%   LinearSolver was set to 'GMRES' or 'BiCGStab' and if the property 
%   PrecType is not 'None'.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$
