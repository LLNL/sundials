function [status, new_data] = kindiagp_pset(y,yscale,fy,fscale,data)
%KINDIAGP_PSET - Preconditioner setup function for the KINDIAGP example.
%
%   See also: kindiagp, KINPrecSetupFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

nlocal = data.nlocal;

for i = 1:nlocal
  P(i) = 0.5 / (y(i)+5.0);
end

new_data = data;
new_data.P = P;

status = 0;