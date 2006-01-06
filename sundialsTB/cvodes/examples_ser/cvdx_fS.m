function [ySd, new_data] = cvdx_fS(is,t,y,yd,yS,data)
%CVDX_FS - sensitivity RHS function for the CVFDX example problem.
%
%  See also: cvfdx, CVSensRhs1Fn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

ySd(1) = -r1*yS(1) + r2*y(3)*yS(2) + r2*y(2)*yS(3);
ySd(3) = 2*r3*y(2)*yS(2);
ySd(2) = -ySd(1)-ySd(3);

switch is
 case 1
  ySd(1) = ySd(1) - y(1);
  ySd(2) = ySd(2) + y(1);
 case 2
  ySd(1) = ySd(1) + y(2)*y(3);
  ySd(2) = ySd(2) - y(2)*y(3);
 case 3
  ySd(2) = ySd(2) - y(2)*y(2);
  ySd(3) = ySd(3) + y(2)*y(2);
end

new_data = [];
