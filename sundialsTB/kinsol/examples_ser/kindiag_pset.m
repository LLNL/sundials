function [flag, new_data] = kindiag_pset(y,yscale,fy,fscale,data)
%KINDIAG_PSET - Preconditioner setup function for the KINDIAG example problem.
%
%   See also: kindiag, KINPrecSetupFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 19:00:06 $


neq = length(y);

for i = 1:neq
  P(i) = 0.5 / (y(i)+5.0);
end

new_data.P = P;  % updated P in data structure

flag = 0;        % success