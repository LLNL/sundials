function [ySd, flag, new_data] = cvdx_fS(t,y,yd,yS,data)
%CVDX_FS - sensitivity RHS function for the CVFDX example problem.
%   It returns the right-hand sides of the sensitivity equations
%   with respect to r1 and r2.
%   YS and YSD are NxNs = 3x2 matrices
%
%  See also: cvfdx, CVSensRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/02/13 23:01:27 $


r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

% r1

yS1 = yS(:,1);
yS1d = zeros(3,1);

yS1d(1) = -r1*yS1(1) + r2*y(3)*yS1(2) + r2*y(2)*yS1(3);
yS1d(3) = 2*r3*y(2)*yS1(2);
yS1d(2) = -yS1d(1)-yS1d(3);

yS1d(1) = yS1d(1) - y(1);
yS1d(2) = yS1d(2) + y(1);

% r2

yS2 = yS(:,2);
yS2d = zeros(3,1);

yS2d(1) = -r1*yS1(1) + r2*y(3)*yS2(2) + r2*y(2)*yS2(3);
yS2d(3) = 2*r3*y(2)*yS2(2);
yS2d(2) = -yS2d(1)-yS2d(3);

yS2d(1) = yS2d(1) + y(2)*y(3);
yS2d(2) = yS2d(2) - y(2)*y(3);

% Return values

ySd = [yS1d yS2d];
flag = 0;
new_data = [];
