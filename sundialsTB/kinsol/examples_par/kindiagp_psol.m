function [x, status, new_data] = kindiagp_psol(y,yscale,fy,fscale,v,data)
%KINDIAGP_PSOL - Preconditioner solve function for the KINDIAGP example.
%
%   See also: kindiagp, KINPrecSolveFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

nlocal = data.nlocal;
P = data.P;

for i=1:nlocal
  x(i) = v(i) * P(i);
end

status = 0;
new_data = [];