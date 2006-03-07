function [] = CVodeMallocB(fctB,tB0,yB0,varargin)
%CVodeMallocB allocates and initializes backward memory for CVODES.
%
%   Usage:   CVodeMallocB ( FCTB, TB0, YB0 [, OPTIONSB] )
%
%   FCTB     is a function defining the adjoint ODE right-hand side.
%            This function must return a vector containing the current 
%            value of the adjoint ODE righ-hand side.
%   TB0      is the final value of t.
%   YB0      is the final condition vector yB(tB0).  
%   OPTIONSB is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   See also: CVRhsFn 
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:41 $

mode = 4;

if nargin < 3
  disp('CvodeMallocB:: too few parameters');
  return
end

options = [];
if nargin == 4
  optionsB = varargin{1};
end

cvm(mode,fctB,tB0,yB0,optionsB);
