function CVodeAdjInit(steps, interp)
%CVodeAdjInit allocates and initializes memory for ASA with CVODES.
%
%   Usage: CVodeAdjInit(STEPS, INTEPR) 
%
%   STEPS    specifies the (maximum) number of integration steps between two 
%            consecutive check points.
%   INTERP   Specifies the type of interpolation used for estimating the forward 
%            solution during the backward integration phase. INTERP should be
%            'Hermite', indicating cubic Hermite interpolation, or 'Polynomial',
%            indicating variable order polynomial interpolation.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/15 19:31:25 $

mode = 4;

if nargin ~= 2
  error('Wrong number of input arguments');
end

cvm(mode,steps,interp);
