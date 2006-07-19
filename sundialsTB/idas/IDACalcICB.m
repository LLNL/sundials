function [status, varargout] = IDACalcICB(tout,icmeth)
%IDACalcICB computes consistent initial conditions for the backward phase.
%
%   Usage: STATUS = IDACalcICB ( TOUTB, ICMETHB )
%          [STATUS, YY0B, YP0B] = IDACalcIC ( TOUTB, ICMETHB )
%
%  See also: IDASetOptions, IDAMallocB, IDAReInitB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $

mode = 8;

if nargout == 1
  status = idm(mode, tout, icmeth);
elseif nargout == 3
  [status, yy, yp] = idm(mode, tout, icmeth);
  varargout(1) = {yy};
  varargout(2) = {yp};
else
  disp('IDACalcICB:: wrong number of output arguments');
end
