function [x, status, new_data] = kindiag_psol(y,yscale,fy,fscale,v,data)
%KINDIAG_PSOL - Preconditioner solve function for the KINDIAG example problem.
%
%   See also: kindiag, KINPrecSolveFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


P = data.P;

neq = length(y);

for i=1:neq
  x(i) = v(i) * P(i);
end

status = 0;
new_data = [];