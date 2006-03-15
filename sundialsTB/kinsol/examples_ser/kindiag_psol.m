function [x, flag, new_data] = kindiag_psol(y,yscale,fy,fscale,v,data)
%KINDIAG_PSOL - Preconditioner solve function for the KINDIAG example problem.
%
%   See also: kindiag, KINPrecSolveFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 19:00:06 $


P = data.P;

neq = length(y);

for i=1:neq
  x(i) = v(i) * P(i);
end

flag = 0;      % success
new_data = []; % data was not modified