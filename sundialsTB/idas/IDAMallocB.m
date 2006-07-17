function [] = IDAMallocB(fctB,tB0,yyB0,ypB0,varargin)
%IDAMallocB allocates and initializes backward memory for IDAS.
%
%   Usage:   IDAMallocB ( FCTB, TB0, YYB0, YPB0 [, OPTIONSB] )
%
%   FCTB     is a function defining the adjoint DAE residual.
%            This function must return a vector containing the current 
%            value of the adjoint DAE residual.
%   TB0      is the final value of t.
%   YYB0     is the final condition vector yB(tB0).  
%   YPB0     is the final condition vector yB'(tB0).  
%   OPTIONSB is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%
%   See also: IDAResFn 
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

mode = 4;

if nargin < 4
  disp('IDAMallocB:: too few parameters');
  return
end

options = [];
if nargin == 5
  optionsB = varargin{1};
end

idm(mode,fctB,tB0,yyB0,ypB0,optionsB);
