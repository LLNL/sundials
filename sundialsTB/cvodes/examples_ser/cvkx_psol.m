function [z, status, new_data] = cvkx_psol(t,y,fy,r,data)
%CVKX_PSOL - Preconditioner solve function for the CVKX example problem.
%
%   See also: cvkx, CVPrecSolveFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

P = data.P;
 
mx = data.mx;
my = data.my;

r = reshape(r,2,mx*my);

for i = 1:mx*my
  z(:,i) = P(:,:,i)^(-1)*r(:,i);
end

z = reshape(z,2*mx*my,1);

status = 0;

new_data = [];

