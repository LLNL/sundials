function [J, new_data] = cvbx_J(t, y, fy, data)
%CVBX_J - Jacobian functino for the CVBX example problem.
%
%   See also: cvbx, CVBandJacFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


mx = data.mx;
my = data.my;
hordc = data.hdcoef;
horac = data.hacoef;
verdc = data.vdcoef;

mu = my;
ml = my;
mband = mu + 1 + ml;

for i = 1:mx
  for j = 1:my
     k = j + (i-1)*my;
     J(mu+1,k) = -2.0 * (verdc + hordc);
     if  i ~= 1
       J(1,k) = hordc + horac;
     end
     if i ~= mx
       J(mband,k) = hordc - horac;
     end
     if j ~= 1
       J(mu,k) = verdc;
     end
     if j ~= my
       J(mu+2,k) = verdc;
     end
  end
end

new_data = [];
