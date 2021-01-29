function [varargout] = IDASolve(tout,itask)
%IDASolve integrates the DAE.
%
%   Usage: [STATUS, T, Y] = IDASolve ( TOUT, ITASK ) 
%          [STATUS, T, Y, YQ] = IDASolve  (TOUT, ITASK )
%          [STATUS, T, Y, YS] = IDASolve ( TOUT, ITASK )
%          [STATUS, T, Y, YQ, YS] = IDASolve ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and returns 
%   Y(TOUT). If ITASK is 'OneStep', then the solver takes one internal time step 
%   and returns in Y the solution at the new internal time. In this case, TOUT 
%   is used only during the first call to IDASolve to determine the direction of 
%   integration and the rough scale of the problem. In either case, the time 
%   reached by the solver is returned in T.
%
%   If quadratures were computed (see IDAQuadInit), IDASolve will return their
%   values at T in the vector YQ.
%
%   If sensitivity calculations were enabled (see IDASensInit), IDASolve will 
%   return their values at T in the matrix YS. Each row in the matrix YS
%   represents the sensitivity vector with respect to one of the problem parameters.
%
%   In ITASK =' Normal' mode, to obtain solutions at specific times T0,T1,...,TFINAL
%   (all increasing or all decreasing) use TOUT = [T0 T1  ... TFINAL]. In this case
%   the output arguments Y and YQ are matrices, each column representing the solution
%   vector at the corresponding time returned in the vector T. If computed, the 
%   sensitivities are eturned in the 3-dimensional array YS, with YS(:,:,I) representing
%   the sensitivity vectors at the time T(I).
%
%   On return, STATUS is one of the following:
%     0: IDASolve succeeded and no roots were found.
%     1: IDASolve succeded and returned at tstop.
%     2: IDASolve succeeded, and found one or more roots. 
%    -1: An error occurred (see printed message).
%
%   See also IDASetOptions, IDAGetStats

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2011/05/26 00:05:36 $


mode = 20;

if nargin ~= 2
  error('Wrong number of input arguments');
end

if nargout < 3 || nargout > 5
  error('Wrong number of output arguments');
end

varargout = cell (nargout, 1);

[varargout{:}] = idm(mode,tout,itask);
