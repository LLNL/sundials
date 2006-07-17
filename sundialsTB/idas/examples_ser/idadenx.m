%IDADENX - IDAS example problem (serial, dense)
%
%   See also: idadenx_f, idadenx_g, idadenx_J

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/15 19:31:26 $

data.p = [0.04; 1.0e4; 3.0e7];

t0 = 0.0;
y0 = [1.0;0.0;0.0];
yp0 = [-0.04;0.04;0.0];

options = IDASetOptions('RelTol',1.e-4,...
                        'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                        'LinearSolver','Dense',...
                        'JacobianFn',@idadenx_J);

options = IDASetOptions(options,'RootsFn',@idadenx_g, 'NumRoots',2);

mondata.sol = true;
mondata.update = 100;
options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);

IDAMalloc(@idadenx_f,t0,y0,yp0,options,data);

t1 = 0.4;
tmult = 10.0;
nout = 12;

fprintf('-----------------------------------------------------------------------\n');
fprintf('  t             y1           y2           y3');
fprintf('      | nst  k      h\n');
fprintf('-----------------------------------------------------------------------\n');

iout = 0;
tout = t1;
while iout < nout

  [status,t,y,yp] = IDASolve(tout,'Normal');
  
% Extract statistics
  si = IDAGetStats;

% Print output
  if(status == 2)
    fprintf(' ... Root found  %d   %d\n',si.RootInfo.roots(1), si.RootInfo.roots(2));
  end
  fprintf('%10.4e %12.4e %12.4e %12.4e | %3d  %1d %12.4e\n',... 
          t, y(1), y(2), y(3), si.nst, si.qlast, si.hlast);

% Update output time
  if(status == 0)
    iout = iout+1;
    tout = tout*tmult;
  end
  
end

si = IDAGetStats;

IDAFree;

