function [] = CVodeReInit(fct,t0,y0,varargin)
%CVodeReInit reinitializes memory for CVODES
%   where a prior call to CVodeMalloc has been made with the same
%   problem size N. CVodeReInit performs the same input checking
%   and initializations that CVodeMalloc does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage: CVodeReInit ( ODEFUN, T0, Y0 [, OPTIONS [, DATA] ] ) 
%
%   ODEFUN   is a function defining the ODE right-hand side: y' = f(t,y).
%            This function must return a vector containing the current 
%            value of the righ-hand side.
%   T0       is the initial value of t.
%   Y0       is the initial condition vector y(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%   DATA     is (optional) problem data passed unmodified to all
%            user-provided functions when they are called. For example,
%            YD = ODEFUN(T,Y,DATA).
%
%   See also: CVodeMalloc, CVRhsFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/07 19:08:40 $

mode = 11;

if nargin < 3
  disp('CVodeReInit:: too few parameters');
  return
end

options = [];
data =[];
if nargin > 3
  options = varargin{1};
end
if nargin > 4
  data = varargin{2};
end

cvm(mode,fct,t0,y0,options,data);
