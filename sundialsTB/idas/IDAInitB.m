function idxB = IDAInitB(fctB, tB0, yyB0, ypB0, optionsB)
%IDAInitB allocates and initializes backward memory for CVODES.
%
%   Usage:   IDXB = IDAInitB ( DAEFUNB, TB0, YYB0, YPB0 [, OPTIONSB] )
%
%   DAEFUNB  is a function defining the adjoint DAE: F(t,y,y',yB,yB')=0
%            This function must return a vector containing the current 
%            value of the adjoint DAE residual.
%   TB0      is the final value of t.
%   YYB0     is the final condition vector yB(tB0).  
%   YPB0     is the final condition vector yB'(tB0).  
%   OPTIONSB is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%
%   IDAInitB returns the index IDXB associated with this backward
%   problem. This index must be passed as an argument to any subsequent
%   functions related to this backward problem.
%
%   See also: IDASetOptions, IDAResFnB
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2007/05/11 18:51:32 $

mode = 5;

if nargin < 4
  error('Too few input arguments');
end

if nargin < 5
  optionsB = [];
end

idxB = idm(mode, fctB, tB0, yyB0, ypB0, optionsB);
idxB = idxB+1;
