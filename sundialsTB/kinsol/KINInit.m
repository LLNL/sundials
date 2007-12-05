function status = KINInit(fct, n, options)
%KINInit allocates and initializes memory for KINSOL.
%
%   Usage:   KINInit ( SYSFUN, N [, OPTIONS ] );
%
%   SYSFUN   is a function defining the nonlinear problem f(y) = 0.
%            This function must return a column vector FY containing the
%            current value of the residual
%   N        is the (local) problem dimension.
%   OPTIONS  is an (optional) set of integration options, created with
%            the KINSetOptions function. 
%
%   See also: KINSetOptions, KINSysFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/01/06 19:00:02 $

mode = 1;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = kim(mode, fct, n, options);
