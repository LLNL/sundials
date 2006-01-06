%CVFDX - CVODES forward sensitivity example (serial, dense)
%   The following is a simple example problem, with the coding
%   needed for its solution by CVODES. The problem is from chemical
%   kinetics, and consists of the following three rate equations:
%      dy1/dt = -p1*y1 + p2*y2*y3
%      dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
%      dy3/dt =  p3*(y2)^2
%   on the interval from t = 0.0 to t = 4.e10, with initial
%   conditions y1 = 1.0, y2 = y3 = 0. The reaction rates are: p1=0.04,
%   p2=1e4, and p3=3e7. The problem is stiff.
%   This program solves the problem with the BDF method, Newton
%   iteration with the CVODES dense linear solver, and a
%   user-supplied Jacobian routine.
%   It uses a scalar relative tolerance and a vector absolute
%   tolerance.
%   Output is printed in decades from t = .4 to t = 4.e10.
%   Run statistics (optional outputs) are printed at the end.
%
%   Optionally, CVODES can compute sensitivities with respect to the
%   problem parameters p1, p2, and p3.
%   The sensitivity right hand side is given analytically through the
%   user routine fS (of type SensRhs1Fn).
%   Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
%   STAGGERED1) can be used and sensitivities may be included in the
%   error test or not (error control set on TRUE or FALSE,
%   respectively).
%
%   See also: cvdx_f, cvdx_J, cvdx_fS

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

% User data structure

data.p = [0.04; 1.0e4; 3.0e7];

% CVODES options

t0 = 0.0;
y0 = [1.0;0.0;0.0];

options = CVodeSetOptions('RelTol',1.e-4,...
                          'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                          'JacobianFn',@cvdx_J);
yS0 = zeros(3,3);

% Case 1: user-provided sensitivity RHS

options = CVodeSetOptions(options,...
                          'SensiAnalysis', 'FSA',...
                          'FSAInitCond', yS0,...
                          'FSAMethod', 'Staggered1',...
                          'FSAErrControl', 'on',...
                          'FSAParamScales', [0.04; 1.0e4; 3.0e7],...
                          'FSARhsType', 'One', 'FSARhsFn', @cvdx_fS);

% Case 2: inernal DQ approximation

%options = CVodeSetOptions(options,...
%                             'ForwardSensi', 'on',...
%                             'FSAInitCond', yS0,...
%                             'FSAMethod', 'Staggered1',...
%                             'FSAErrControl', 'on',...
%                             'ParamField', 'p',...
%                             'ParamList', [1 2 3],...
%                             'ParamScales', [0.04 1.0e4 3.0e7]);

options = CVodeSetOptions(options,'MonitorFn','CVodeMonitor');

CVodeMalloc(@cvdx_f,t0,y0,options,data);

% Problem solution

t1 = 0.4;
tmult = 10.0;
nout = 12;

iout = 0;
tout = t1;
while 1,
  [status,t,y, yS] = CVode(tout,'Normal');
  fprintf('At t = %0.4e  status = %d\n',t,status);
  fprintf('    y = %14.6e  %14.6e  %14.6e\n', y(1), y(2), y(3));
  fprintf('  yS1 = %14.6e  %14.6e  %14.6e\n', yS(1,1), yS(2,1), yS(3,1));
  fprintf('  yS2 = %14.6e  %14.6e  %14.6e\n', yS(1,2), yS(2,2), yS(3,2));
  fprintf('  yS3 = %14.6e  %14.6e  %14.6e\n', yS(1,3), yS(2,3), yS(3,3));
  if(status ==0)
    iout = iout+1;
    tout = tout*tmult;
  end
  
  if iout==nout
    break;
  end
end

si = CVodeGetStats

% Free memory

CVodeFree;



