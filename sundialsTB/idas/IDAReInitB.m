function [] = IDAReInitB(fctB,tB0,yyB0,ypB0,varargin)
%IDAReInitB allocates and initializes backward memory for IDAS.
%   where a prior call to IDAMallocB has been made with the same
%   problem size NB. IDAReInitB performs the same input checking
%   and initializations that IDAMallocB does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage:   IDAReInitB ( FCTB, TB0, YYB0, YPB0 [, OPTIONSB] )
%
%   FCTB     is a function defining the adjoint ODE right-hand side.
%            This function must return a vector containing the current 
%            value of the adjoint ODE righ-hand side.
%   TB0      is the final value of t.
%   YYB0     is the final condition vector yB(tB0).  
%   YPB0     is the final condition vector yB'(tB0).
%   OPTIONSB is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%
%   See also: IDAMallocB, IDAResFn 
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/11/25 19:57:25 $

mode = 14;

if nargin < 4
  disp('IDAReInitB:: too few parameters');
  return
end

options = [];
if nargin == 5
  optionsB = varargin{1};
end

idm(mode,fctB,tB0,yyB0,ypB0,optionsB);
