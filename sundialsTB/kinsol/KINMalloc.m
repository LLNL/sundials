function [] = KINMalloc(fct,n,varargin)
%KINMalloc allocates and initializes memory for KINSOL.
%
%   Usage:   KINMalloc ( SYSFUN, N [, OPTIONS [, DATA] ] );
%
%   SYSFUN   is a function defining the nonlinear problem f(y) = 0.
%            This function must return a column vector FY containing the
%            current value of the residual
%   N        is the (local) problem dimension.
%   OPTIONS  is an (optional) set of integration options, created with
%            the KINSetOptions function. 
%   DATA     is the (optional) problem data passed unmodified to all
%            user-provided functions when they are called. For example,
%            RES = SYSFUN(Y,DATA).
%
%   See also: KINSysFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

mode = 1;

if nargin < 2
  disp('KINMalloc:: too few parameters');
  return
end

options = [];
data =[];
if nargin > 2
  options = varargin{1};
end
if nargin > 3
  data = varargin{2};
end

kim(mode,fct,n,options,data);
