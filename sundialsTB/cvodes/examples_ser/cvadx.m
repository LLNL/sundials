%CVADX - CVODES adjoint sensitivity example problem (serial, dense)
%   The following is a simple example problem, with the coding
%   needed for its solution by CVODES. The problem is from chemical
%   kinetics, and consists of the following three rate equations.
%      dy1/dt = -p1*y1 + p2*y2*y3
%      dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
%      dy3/dt =  p3*(y2)^2
%   on the interval from t = 0.0 to t = 4.e10, with initial
%   conditions: y1 = 1.0, y2 = y3 = 0. The reaction rates are:
%   p1=0.04, p2=1e4, and p3=3e7. The problem is stiff.
%   This program solves the problem with the BDF method, Newton
%   iteration with the CVODE dense linear solver, and a user-supplied
%   Jacobian routine.
%   It uses a scalar relative tolerance and a vector absolute
%   tolerance.
%   Output is printed in decades from t = .4 to t = 4.e10.
%   Run statistics (optional outputs) are printed at the end.
%   
%   Optionally, CVODES can compute sensitivities with respect to
%   the problem parameters p1, p2, and p3 of the following quantity:
%     G = int_t0^t1 g(t,p,y) dt
%   where
%     g(t,p,y) = y3
%          
%   The gradient dG/dp is obtained as:
%     dG/dp = int_t0^t1 (g_p - lambda^T f_p ) dt - lambda^T(t0)*y0_p
%           = - xi^T(t0) - lambda^T(t0)*y0_p
%   where lambda and xi are solutions of:
%     d(lambda)/dt = - (f_y)^T * lambda - (g_y)^T
%     lambda(t1) = 0
%   and
%     d(xi)/dt = - (f_p)^T * lambda + (g_p)^T
%     xi(t1) = 0
%   
%   During the backward integration, CVODES also evaluates G as
%     G = - phi(t0)
%   where
%     d(phi)/dt = g(t,y,p)
%     phi(t1) = 0
%
%   See also: cvdx_f, cvdx_q, cvdx_J, cvdx_fB, cvdx_qB, cvdx_JB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


% ----------------------------------------
% User data structure
% ----------------------------------------

data.p = [0.04; 1.0e4; 3.0e7];

% ----------------------------------------
% Forward CVODES options
% ----------------------------------------

t0 = 0.0;
y0 = [1.0;0.0;0.0];
q0 = 0.0;

options = CVodeSetOptions('RelTol',1.e-4,...
                          'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                          'LinearSolver','Dense',...
                          'JacobianFn',@cvdx_J);

options = CVodeSetOptions(options,...
                          'Quadratures','on',...
                          'QuadRhsFn',@cvdx_q,...
                          'QuadInitcond', q0,...
                          'QuadErrControl','on',...
                          'QuadRelTol',1.e-4,'QuadAbsTol',1.e-6);

options = CVodeSetOptions(options,...
                          'SensiAnalysis', 'ASA',...
                          'ASANumPoints', 150);

%options = CVodeSetOptions(options,...
%                          'MonitorFn','CVodeMonitor');
%
CVodeMalloc(@cvdx_f,t0,y0,options,data);

% ----------------------------------------
% Forward integration
% ----------------------------------------

tout = 4.e7;
[status,t,y,q] = CVode(tout,'Normal');
s = CVodeGetStats;
fprintf('G = %12.4e   (%d steps)\n',q, s.nst);

ck = CVodeGet('CheckPointsInfo');
fprintf('\n');
fprintf([' Address     Next'...
         '      t0         t1'...
         '     nstep  order'...
         '  step size\n']); 
for i = 1:length(ck)
  fprintf('%8x %8x  %8.3e  %8.3e  %4d     %1d   %10.5e\n',...
          ck(i).addr, ck(i).next, ck(i).t0, ck(i).t1, ck(i).nstep, ...
          ck(i).order, ck(i).step);
end
fprintf('\n');

% ----------------------------------------
% Backward CVODES options
% ----------------------------------------

tB1 = 4.e7;
yB1 = [0.0;0.0;0.0];
qB1 = [0.0;0.0;0.0];

optionsB = CVodeSetOptions('RelTol',1.e-6,...
                           'AbsTol',1.e-3,...
                           'LinearSolver','Dense',...
                           'JacobianFn',@cvdx_JB);

optionsB = CVodeSetOptions(optionsB,...
                           'Quadratures','on',...
                           'QuadRhsFn',@cvdx_qB,...
                           'QuadInitcond', qB1,...
                           'QuadErrControl','on',...
                           'QuadRelTol',1.e-6,'QuadAbsTol',1.e-3);
mondataB.dir = -1;
mondataB.cntr = false;
mondataB.sol = true;
optionsB = CVodeSetOptions(optionsB,...
                           'MonitorFn','CVodeMonitor',...
                           'MonitorData', mondataB);

CVodeMallocB(@cvdx_fB, tB1, yB1, optionsB);


fprintf('Start backward integration\n');

% ----------------------------------------
% Backward integration
% ----------------------------------------

[status,t,yB,qB] = CVodeB(t0,'Normal');

fprintf('tB1:        %12.4e\n',tB1);
fprintf('dG/dp:      %12.4e  %12.4e  %12.4e\n',...
        -qB(1),-qB(2),-qB(3));
fprintf('lambda(t0): %12.4e  %12.4e  %12.4e\n',...
        yB(1),yB(2),yB(3));



sB=CVodeGetStatsB;

% ----------------------------------------
% Free memory
% ----------------------------------------

CVodeFree;
