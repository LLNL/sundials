function [varargout] = CVodeB(tout,itask)
%CVodeB integrates all backwards ODEs currently defined.
%
%   Usage:  [STATUS, T, YB] = CVodeB ( TOUT, ITASK ) 
%           [STATUS, T, YB, YQB] = CVodeB ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and returns 
%   YB(TOUT). If ITASK is 'OneStep', then the solver takes one internal time step 
%   and returns in YB the solution at the new internal time. In this case, TOUT 
%   is used only during the first call to CVodeB to determine the direction of 
%   integration and the rough scale of the problem. In either case, the time 
%   reached by the solver is returned in T. 
%
%   If quadratures were computed (see CVodeQuadInitB), CVodeB will return their
%   values at T in the vector YQB.
%
%   In ITASK =' Normal' mode, to obtain solutions at specific times T0,T1,...,TFINAL
%   (all increasing or all decreasing) use TOUT = [T0 T1  ... TFINAL]. In this case
%   the output arguments YB and YQB are matrices, each column representing the solution
%   vector at the corresponding time returned in the vector T.
%
%   If more than one backward problem was defined, the return arguments are cell
%   arrays, with T{IDXB}, YB{IDXB}, and YQB{IDXB} corresponding to the backward
%   problem with index IDXB (as returned by CVodeInitB).
%
%   On return, STATUS is one of the following:
%     0: successful CVodeB return.
%     1: CVodeB succeded and return at a tstop value (internally set).
%    -1: an error occurred (see printed message).
%
%   See also CVodeSetOptions, CVodeGetStatsB

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:42:38 $

mode = 21;

if nargin ~= 2
  error('Wrong number of input arguments');
end

if nargout < 3 || nargout > 4
  error('Wrong number of output arguments');
end

varargout = cell (nargout, 1);
[varargout{:}] = cvm(mode,tout,itask);
