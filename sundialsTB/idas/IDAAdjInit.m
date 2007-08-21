function [] = IDAAdjInit(steps, interp)
%IDAAdjInit allocates and initializes memory for ASA with IDAS.
%
%   Usage: IDAAdjInit(STEPS, INTEPR) 
%
%   STEPS    specifies the (maximum) number of integration steps between two 
%            consecutive check points.
%   INTERP   Specifies the type of interpolation used for estimating the forward 
%            solution during the backward integration phase. INTERP should be
%            'Hermite', indicating cubic Hermite interpolation, or 'Polynomial',
%            indicating variable order polynomial interpolation.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $

mode = 4;

if nargin ~= 2
  error('Wrong number of input arguments');
end

idm(mode,steps,interp);
