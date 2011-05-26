function status = IDAInit(fct,t0,yy0,yp0,options)
%IDAInit allocates and initializes memory for IDAS.
%
%   Usage: IDAInit ( DAEFUN, T0, YY0, YP0 [, OPTIONS ] ) 
%
%   DAEFUN   is a function defining the DAE residual: f(t,yy,yp).
%            This function must return a vector containing the current 
%            value of the residual.
%   T0       is the initial value of t.
%   YY0      is the initial condition vector y(t0).  
%   YP0      is the initial condition vector y'(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%
%  See also: IDASetOptions, IDAResFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/12/05 21:58:18 $

mode = 1;

if nargin < 4
  error('Too few input arguments');
end

if nargin < 5
  options = [];
end

status = idm(mode, fct, t0, yy0, yp0, options);
