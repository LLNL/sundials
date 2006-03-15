function [] = CVadjMalloc(steps, interp)
%CVadjMalloc allocates and initializes memory for ASA with CVODES.
%
%   Usage: CVadjMalloc(STEPS, INTEPR) 
%
%   STEPS    specifies the (maximum) number of integration steps between two 
%            consecutive check points.
%   INTERP   Specifies the type of interpolation used for estimating the forward 
%            solution during the backward integration phase. INTERP should be
%            'Hermite', indicating cubic Hermite interpolation, or 'Polynomial',
%            indicating variable order polynomial interpolation.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/03/07 01:19:50 $

mode = 3;

cvm(mode,steps,interp);
