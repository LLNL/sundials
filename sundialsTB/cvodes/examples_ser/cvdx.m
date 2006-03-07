%CVDX - CVODES example problem (serial, dense)
%   The following is a simple example problem, with the coding
%   needed for its solution by CVODES. The problem is from
%   chemical kinetics, and consists of the following three rate
%   equations:         
%      dy1/dt = -.04*y1 + 1.e4*y2*y3
%      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
%      dy3/dt = 3.e7*(y2)^2
%   on the interval from t = 0.0 to t = 4.e10, with initial
%   conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
%   While integrating the system, we also use the rootfinding
%   feature to find the points at which y1 = 1e-4 or at which
%   y3 = 0.01. This program solves the problem with the BDF method,
%   Newton iteration with the CVDENSE dense linear solver, and a
%   user-supplied Jacobian routine.
%   It uses a scalar relative tolerance and a vector absolute
%   tolerance. Output is printed in decades from t = .4 to t = 4.e10.
%   Run statistics (optional outputs) are printed at the end.
%
%   See also: cvdx_f, cvdx_g, cvdx_J

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:48 $

data.p = [0.04; 1.0e4; 3.0e7];

t0 = 0.0;
y0 = [1.0;0.0;0.0];

options = CVodeSetOptions('RelTol',1.e-8,...
                          'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                          'LinearSolver','Dense',...
                          'JacobianFn',@cvdx_J,...
                          'RootsFn',@cvdx_g, 'NumRoots',2);

mondata.sol = true;
mondata.mode = 'text';
mondata.skip = 10;
mondata.update = 100;
options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor','MonitorData',mondata);

CVodeMalloc(@cvdx_f,t0,y0,options,data);

t1 = 0.4;
tmult = 10.0;
nout = 10;

fprintf('---------------------------------------------------\n');
fprintf('     t        q       h\n');
fprintf(['          y   '...
         '          yd1'...
         '          yd2\n']);
fprintf('---------------------------------------------------\n');

iout = 0;
tout = t1;
while iout < nout

  [status,t,y] = CVode(tout,'Normal');
  
% Compute solution derivative in two ways:
  [yd1, new_data] = cvdx_f(t, y, data);
  yd2 = CVodeGet('DerivSolution', t, 1);

% Extract statistics
  si = CVodeGetStats;

% Print output
  fprintf('%0.4e    %1d  %0.4e',t, si.qlast, si.hlast);
  if(status == 2)
    fprintf(' ... Root found  %d   %d\n',si.RootInfo.roots(1), si.RootInfo.roots(2));
  else
    fprintf('\n');
  end
  for i = 1:3
    fprintf('  %14.6e  %14.6e  %14.6e\n', y(i), yd1(i), yd2(i));
  end
  fprintf('---------------------------------------------------\n');

% Update output time
  if(status == 0)
    iout = iout+1;
    tout = tout*tmult;
  end
  
end

si = CVodeGetStats;

CVodeFree;

