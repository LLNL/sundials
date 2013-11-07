%IDAPrecSolveFn - type for preconditioner solve function.
%
%   The user-supplied preconditioner solve function PSOLFUN
%   is to solve a linear system P z = r, where P is the
%   preconditioner matrix.
%
%   The function PSOLFUN must be defined as 
%        FUNCTION [Z, FLAG] = PSOLFUN(T,YY,YP,RR,R)
%   and must return a vector Z containing the solution of Pz=r.
%   If PSOLFUN was successful, it must return FLAG=0. For a recoverable 
%   error (in which case the step will be retried) it must set FLAG to a 
%   positive value. If an unrecoverable error occurs, it must set FLAG
%   to a negative value, in which case the integration will be halted.
%   The input argument RR contains the current value of f(t,yy,yp).
%
%   If a user data structure DATA was specified in IDAInit, then
%   PSOLFUN must be defined as
%        FUNCTION [Z, FLAG, NEW_DATA] = PSOLFUN(T,YY,YP,RR,R,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector Z and
%   the flag FLAG, the PSOLFUN function must also set NEW_DATA. Otherwise,
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would
%   lead to unnecessary copying).
%
%   See also IDAPrecSetupFn, IDASetOptions
%
%   NOTE: PSOLFUN and PSOLFUNB are specified through the property
%   PrecSolveFn to IDASetOptions and are used only if the property
%   LinearSolver was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

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
% $Revision: 1.3 $Date: 2007/08/21 17:38:44 $
