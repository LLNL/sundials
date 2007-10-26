function mcvsDiscRHS_dns()
%mcvsDiscRHS_dns - CVODES example with RHS discontinuity
%  Trivial CVODES example to illustrate the proper 
%  way to integrate over a discontinuity in the RHS: 
%       y' = -y   ; y(0) = 1    ; t = [0,1]
%       z' = -5*z ; z(1) = y(1) ; t = [1,2]
%  The problem is solved twice, first by explicitly treating the 
%  discontinuity point and secondly by letting the integrator 
%  deal with the discontinuity.


% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $

t0 = 0.0;
t1 = 1.0;
t2 = 2.0;

y0 = 1.0;

% ---------------------------------------------------------------
% Discontinuity in RHS: Case 1 - let CVODES deal with it.
% ---------------------------------------------------------------

data.tdisc = t1;

% Initialize solver
options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-3,...
                          'AbsTol',1.e-4,...
                          'LinearSolver','Dense');
CVodeInit(@rhsfn1, 'BDF', 'Newton', t0, y0, options);

% Integrate over the point of discontinuity
t = t0;
i = 1;
tt1(1) = t0; yy1(1) = y0;
while t < t2
  [status, t, y] = CVode(t2,'OneStep');
  i = i+1;
  tt1(i) = t;
  yy1(i) = y;
end

% Free memory
CVodeFree;

% -------------------------------------------------------------
% Discontinuity in RHS: Case 1 - explicit treatment
% Note that, since we set tstop at the point of discontinuity,
% we could simply use the exact same RHS function as before. 
% However, we chose to use a flag set in the user data (to also
% illustrate the use of CVodeSet).
% -------------------------------------------------------------

% Initialize solver (data.flag = 1)

data.flag = 1;

options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-3,...
                          'AbsTol',1.e-4,...
                          'StopTime',t1,...
                          'LinearSolver','Dense');
CVodeInit(@rhsfn2, 'BDF', 'Newton', t0, y0, options);

% Integrate to the point of discontinuity
t = t0;
i = 1;
tt2(1) = t0; yy2(1) = y0;
while t < t2
  [status, t, y] = CVode(t2,'OneStep');
  i = i+1;
  tt2(i) = t;
  yy2(i) = y;
  if status == 1
    % Once tstop is reached, flag a change in RHS
    data.flag = 2;
    CVodeSet('UserData',data);
  end
end

% Free memory
CVodeFree;


% Plot the two solutions

figure

subplot(2,1,1)
hold on
plot(tt1,yy1);
plot(tt2,yy2,'r');
legend('Integrate over discontinuity','Stop at discontinuity');
title('Discontinuity in RHS');
xlabel('time');
ylabel('y');
box on

subplot(2,1,2)
hold on
plot(tt1,yy1,'b',tt1,yy1,'b.');
plot(tt2,yy2,'r',tt2,yy2,'r.');
set(gca,'XLim',[0.99 1.01],'YLim',[0.36 0.37]);
legend('Integrate over discontinuity','Stop at discontinuity');
title('Zoom on discontinuity');
xlabel('time');
ylabel('y');
grid on
box on

return

% ===========================================================================

function [yd, flag, new_data] = rhsfn1(t, y, data)
% Right-hand side function for case 1


if t <= data.tdisc
  yd(1) = -y(1);
else  
  yd(1) = -5*y(1);
end

flag = 0;
new_data = [];

return

% ===========================================================================


function [yd, flag, new_data] = rhsfn2(t, y, data)
% Right-hand side function for case 2

if data.flag == 1
  yd(1) = -y(1);
else  
  yd(1) = -5*y(1);
end

flag = 0;
new_data = [];

return


