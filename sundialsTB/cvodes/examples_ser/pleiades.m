%PLEIADES - CVODES example problem (serial, nonstiff)
%
%   See also: pleiades_f, pleiades_J

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.4 $Date: 2006/03/15 19:31:26 $

iter = 'Newton';
reltol = 1.0e-7;
abstol = 1.0e-7;
itask = 'OneStep';

data = [];
neq = 28;

t0 = 0.0;
tout = 3.0;

y0 = zeros(neq,1);
y0(1)  =  3.0;
y0(2)  =  3.0;
y0(3)  = -1.0;
y0(4)  = -3.0;
y0(5)  =  2.0;
y0(6)  = -2.0;
y0(7)  =  2.0;
y0(8)  =  3.0;
y0(9)  = -3.0;
y0(10) =  2.0;
y0(13) = -4.0;
y0(14) =  4.0;
y0(20) =  1.75;
y0(21) = -1.5;
y0(25) = -1.25;
y0(26) =  1.0;

options = CVodeSetOptions('LMM', 'Adams',...
                          'NonlinearSolver', iter,...
                          'RelTol', reltol,...
                          'AbsTol', abstol,...
                          'StopTime',tout,...
                          'MaxNumSteps',2000);

%if strcmp(iter,'Newton')
%  options = CVodeSetOptions(options,'Jacobian',@pleiades_J);
%end

mondata.select = [1:14];
mondata.sol = true;
options = CVodeSetOptions(options,...
                          'MonitorFn', @CVodeMonitor,...
                          'MonitorData', mondata);

CVodeMalloc(@pleiades_f,t0,y0,options,data);

if strcmp(itask,'Normal')

  [status,t,y] = CVode(tout,'Normal');

else

  t = t0;
  i = 0;
  while t < tout
    i = i+1;
    [status,t,y] = CVode(tout,'OneStepTstop');
    time(i) = t;
    xx(:,i) = y(1:7);
    yy(:,i) = y(8:14);
  end

end
  
Stats = CVodeGetStats

LSInfo = Stats.LSInfo

CVodeFree;


colors = ['k','b','r','g','c','y','m'];
figure;
for i = 1:7
  plot(xx(i,:),yy(i,:),colors(i));
  hold on;
end
legend('Body 1','Body 2','Body 3','Body 4','Body 5','Body 6','Body 7');
title('Body Trajectories');
xlabel('x');
ylabel('y');
grid on;