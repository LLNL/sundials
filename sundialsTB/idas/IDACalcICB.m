function [status, varargout] = IDACalcICB(tout,icmeth)
%IDACalcICB computes consistent initial conditions for the backward phase.
%
%   Usage: STATUS = IDACalcICB ( TOUTB, ICMETHB )
%          [STATUS, YY0B, YP0B] = IDACalcIC ( TOUTB, ICMETHB )
%
%  See also: IDASetOptions, IDAInitB, IDAReInitB

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
% $Revision$Date: 2007/08/21 17:38:42 $

mode = 26;

if nargout == 1
  status = idm(mode, tout, icmeth);
elseif nargout == 3
  [status, yy, yp] = idm(mode, tout, icmeth);
  varargout(1) = {yy};
  varargout(2) = {yp};
else
  disp('IDACalcICB:: wrong number of output arguments');
end
