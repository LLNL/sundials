function [ud, flag, new_data] = cvbx_f(t, u, data)
%CVBX_F - RHS function for the CVBX example problem
%
%   See also: cvbx, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/02/13 23:01:27 $

mx = data.mx;
my = data.my;
hordc = data.hdcoef;
horac = data.hacoef;
verdc = data.vdcoef;

for j = 1:my
  for i = 1:mx
    uij = u(j+(i-1)*my);
    if j == 1
      udn = 0.0;
    else
      udn = u(j-1+(i-1)*my);
    end
    if j == my
      uup = 0.0;
    else
      uup = u(j+1+(i-1)*my);
    end
    if i == 1
      ult = 0.0;
    else
      ult = u(j+(i-2)*my);
    end
    if i == mx
      urt = 0.0;
    else
      urt = u(j+i*my);
    end
    
    hdiff = hordc * (ult - 2*uij + urt);
    hadv = horac * (urt - ult);
    vdiff = verdc * (uup - 2*uij + udn);
    ud(j+(i-1)*my) = hdiff + hadv + vdiff;
  end
end

flag = 0;
new_data = [];


