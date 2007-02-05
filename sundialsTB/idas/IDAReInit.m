function [] = IDAReInit(fct,t0,yy0,yp0,varargin)
%IDAReInit reinitializes memory for IDAS.
%   where a prior call to IDAMalloc has been made with the same
%   problem size N. IDAReInit performs the same input checking
%   and initializations that IDAMalloc does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage: IDAReInit ( DAEFUN, T0, YY0, YP0 [, OPTIONS [, DATA] ] ) 
%
%   DAEFUN   is a function defining the DAE residual: f(t,yy,yp).
%            This function must return a vector containing the current 
%            value of the residual.
%   T0       is the initial value of t.
%   YY0      is the initial condition vector y(t0).  
%   YP0      is the initial condition vector y'(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%   DATA     is (optional) problem data passed unmodified to all
%            user-provided functions when they are called. For example,
%            YD = DAEFUN(T,YY,YP,DATA).
%
%  See also: IDAMalloc, IDARhsFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date:$

mode = 11;

if nargin < 4
  disp('IDAReInit:: too few parameters');
  return
end

options = [];
data =[];
if nargin > 4
  options = varargin{1};
end
if nargin > 5
  data = varargin{2};
end

idm(mode,fct,t0,yy0,yp0,options,data);
