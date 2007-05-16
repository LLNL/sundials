function CVodeInit(fct, t0, y0, options)
%CVodeInit allocates and initializes memory for CVODES.
%
%   Usage: CVodeInit ( ODEFUN, T0, Y0 [, OPTIONS ] ) 
%
%   ODEFUN   is a function defining the ODE right-hand side: y' = f(t,y).
%            This function must return a vector containing the current 
%            value of the righ-hand side.
%   T0       is the initial value of t.
%   Y0       is the initial condition vector y(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   See also: CVodeSetOptions, CVRhsFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/05/11 18:51:32 $

mode = 1;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  options = [];
end

cvm(mode, fct, t0, y0, options);
