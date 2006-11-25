function [] = CVodeReInitB(fctB,tB0,yB0,varargin)
%CVodeReInitB allocates and initializes backward memory for CVODES.
%   where a prior call to CVodeMallocB has been made with the same
%   problem size NB. CVodeReInitB performs the same input checking
%   and initializations that CVodeMallocB does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage:   CVodeReInitB ( FCTB, TB0, YB0 [, OPTIONSB] )
%
%   FCTB     is a function defining the adjoint ODE right-hand side.
%            This function must return a vector containing the current 
%            value of the adjoint ODE righ-hand side.
%   TB0      is the final value of t.
%   YB0      is the final condition vector yB(tB0).  
%   OPTIONSB is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   See also: CVodeMallocB, CVRhsFn 
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/07 19:08:40 $

mode = 14;

if nargin < 3
  disp('CVodeReInitB:: too few parameters');
  return
end

options = [];
if nargin == 4
  optionsB = varargin{1};
end

cvm(mode,fctB,tB0,yB0,optionsB);
