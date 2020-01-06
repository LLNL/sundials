function mcvsOzone_FSA_dns()
%mcvsOzone_FSA_dns - CVODES example problem (serial, dense)
%        ozone depletion chemical mechanism (3 species)

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/10/26 16:30:47 $


% -------------------
% User data structure
% -------------------

data.p = [1.63e-16 ; 4.66e-16 ; 5.00e-11 ; 2.50e-4];

% ---------------------
% CVODES initialization
% ---------------------

options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-4,...
                          'AbsTol',1.e-8,...
                          'LinearSolver','Dense');

t0 = 0.0;
y0 = [1.0e6 ; 1.0e12 ; 3.7e16];



CVodeInit(@rhsfn, 'BDF', 'Newton', t0, y0, options);


% ------------------
% FSA initialization
% ------------------

Ns = 4;
yS0 = zeros(3,Ns);

FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
                                 'ErrControl', true,...
                                 'ParamField', 'p',...
                                 'ParamScales', data.p);
CVodeSensInit(Ns, [], yS0, FSAoptions);

% ----------------
% Problem solution
% ----------------

time(1,1) = t0;
sol(1,:) = y0';

t = t0;
tf = 0.25;
it = 1;
while t<tf
  it = it+1;
  [status, t, y, yS] = CVode(tf,'OneStep');
  time(it,1) = t;
  sol(it,:) = y';
end

si = CVodeGetStats

% -------------
% Plot solution
% -------------

figure
hold on
plot(time, sol(:,1),'r');
plot(time, sol(:,2),'g');
plot(time, sol(:,3),'b');
set(gca,'YScale','log');
set(gca,'XLim',[t0 tf]);
xlabel('time')

grid on
box on

% -----------
% Free memory
% -----------
CVodeFree;

return

% ===========================================================================

function [yd, flag, new_data] = rhsfn(t, y, data)
% Right-hand side function

p = data.p;

yd(1) = -p(1)*y(1)*y(2) - p(2)*y(1)*y(3) + 2*p(3)*y(2) + p(4)*y(3);
yd(2) = -p(1)*y(1)*y(2) + 2*p(2)*y(1)*y(3) + p(4)*y(3);
yd(3) =  p(1)*y(1)*y(2) - p(2)*y(1)*y(3) - p(4)*y(3);

flag = 0;
new_data = [];

return

% ===========================================================================

function [ySd, flag, new_data] = rhsSfn(t,y,yd,yS,data)
% Sensitivity right-hand side function

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

yS2d(1) = -r1*yS2(1) + r2*y(3)*yS2(2) + r2*y(2)*yS2(3);
yS2d(3) = 2*r3*y(2)*yS2(2);
yS2d(2) = -yS2d(1)-yS2d(3);

yS2d(1) = yS2d(1) + y(2)*y(3);
yS2d(2) = yS2d(2) - y(2)*y(3);

% Return values

ySd = [yS1d yS2d];
flag = 0;
new_data = [];

return
