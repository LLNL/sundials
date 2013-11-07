function [varargout] = CVode(tout, itask)
%CVode integrates the ODE.
%
%   Usage: [STATUS, T, Y] = CVode ( TOUT, ITASK ) 
%          [STATUS, T, Y, YS] = CVode ( TOUT, ITASK )
%          [STATUS, T, Y, YQ] = CVode  (TOUT, ITASK )
%          [STATUS, T, Y, YQ, YS] = CVode ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and returns 
%   Y(TOUT). If ITASK is 'OneStep', then the solver takes one internal time step 
%   and returns in Y the solution at the new internal time. In this case, TOUT 
%   is used only during the first call to CVode to determine the direction of 
%   integration and the rough scale of the problem. In either case, the time 
%   reached by the solver is returned in T.
%
%   If quadratures were computed (see CVodeQuadInit), CVode will return their
%   values at T in the vector YQ.
%
%   If sensitivity calculations were enabled (see CVodeSensInit), CVode will 
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
%     0: successful CVode return.
%     1: CVode succeded and returned at tstop.
%     2: CVode succeeded and found one or more roots.
%    -1: an error occurred (see printed message).
%
%   See also CVodeSetOptions, CVodeGetStats

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
% $Revision: 1.7 $Date: 2007/05/16 17:12:56 $

mode = 20;

if nargin ~= 2
  error('Wrong number of input arguments');
end

if nargout < 3 || nargout > 5
  error('Wrong number of output arguments');
end

varargout = cell (nargout, 1);

[varargout{:}] = cvm(mode,tout,itask);
