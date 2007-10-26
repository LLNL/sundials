function mcvsDiscSOL_dns()
%mcvsDiscSOL_dns - CVODES example with solution discontinuity
%  Trivial CVODES example to illustrate the use of
%  CVodeReInit to integrate over a discontinuity in 
%  the solution:
%       y' = -y   ; y(0) = 1    ; t = [0,1]
%       y' = -y   ; y(1) = 1    ; t = [1,2]
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $

t0 = 0.0;
t1 = 1.0;
t2 = 2.0;

% Initialize solver
y = 1.0;
options = CVodeSetOptions('RelTol',1.e-3,...
                          'AbsTol',1.e-4,...
                          'StopTime',t1,...
                          'LinearSolver','Dense');
CVodeInit(@rhsfn,'BDF','Newton',t0,y,options);

% Integrate to the point of discontinuity
t = t0;
i = 1;
tt(i) = t0; yy(i) = y;
while t < t1
  [status, t, y] = CVode(t1,'OneStep');
  i = i+1;
  tt(i) = t;
  yy(i) = y;
end

% Add discontinuity and reinitialize solver
y = 1.0;
options = CVodeSetOptions(options,'StopTime',t2);
CVodeReInit(t1,y,options);

% Integrate from discontinuity to final time
t = t1;
i = i+1;
tt(i) = t1; yy(i) = y;
while t < t2
  [status, t, y] = CVode(t2,'OneStep');
  i = i+1;
  tt(i) = t;
  yy(i) = y;
end

% Plot solution
figure
plot(tt,yy)
title('Discontinuity in solution');
xlabel('time');
ylabel('y');

% Free memory
CVodeFree;

return

% ===========================================================================

function [yd, flag] = rhsfn(t, y)
% Right-hand side function

yd(1) = -y(1);
flag = 0;

return




