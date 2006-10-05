%VDP - CVODES example problem (serial, dense)
%   This is the van der Pol Matlab example, solved with CVODES.
%
%   See also: vdp_f, vdp_J

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/03/07 01:19:54 $


data.mu = 100.0;

t0 = 0.0;
y0 = [2.0;0.0];

options = CVodeSetOptions('RelTol',1.e-3,...
                          'AbsTol',1e-6,...
                          'JacobianFn',@vdp_J);

mondata.mode = 'both';
mondata.sol = true;
options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,'MonitorData',mondata);

CVodeMalloc(@vdp_f,t0,y0,options,data);

% Using default options
%CVodeMalloc('vdp_f',t0,y0,[],data);

tout = 300.0;
i = 1;
while 1
  [status,t,y] = CVode(tout,'OneStep');
  if status < 0
    fprintf('At t = %f   status = %d\n',t,status);
  end
  tt(i) = t;
  yy(i) = y(1);
  i = i + 1;
  if t > tout
    break;
  end
end
  
CVodeFree;

figure;
plot(tt,yy);
