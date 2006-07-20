%IDADENX - IDAS example problem (serial, dense)
%
%   See also: idadenx_f, idadenx_g, idadenx_J

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/07/19 20:52:29 $

m = 40;
N = m^2;
data.m = m;
data.N = N;
data.dx = 1.0/(m-1);
data.c = 1.0/data.dx^2;


fp = figure;

[t0,yy0,yp0,id,cnstr] = idabanx_ic(data); 

figure(fp);
subplot(2,2,1);
hold on
hs1 = surf(reshape(yy0,m,m));
shading interp
set(hs1,'FaceAlpha',0.35);
box on
view(-30,30)
subplot(2,2,2);
hold on
hs2 = surf(reshape(yp0,m,m));
shading interp
set(hs2,'FaceAlpha',0.35);
box on
view(-30,30)

options = IDASetOptions('RelTol',0.0,...
                        'AbsTol',1.0e-3,...
                        'VariableTypes',id,...
                        'ConstraintTypes',cnstr,...
                        'LinearSolver','Band',...
                        'LowerBwidth',m,...
                        'UpperBwidth',m);

mondata.mode = 'text';
mondata.update = 100;
options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);

IDAMalloc(@idabanx_f,t0,yy0,yp0,options,data);

tout = 0.01;
[status, yy0_mod, yp0_mod] = IDACalcIC(tout, 'FindAlgebraic');

figure(fp);

subplot(2,2,1);
hs1 = surf(reshape(yy0_mod,m,m));
set(hs1,'FaceColor','none');
subplot(2,2,2);
hs2 = surf(reshape(yp0_mod,m,m));
set(hs2,'FaceColor','none');


subplot(2,2,3);
hold on
hs1 = surf(reshape(yy0_mod,m,m));
shading interp
view(-30,30)
zlim_yy = get(gca,'ZLim');
box on
subplot(2,2,4);
hold on
hs2 = surf(reshape(yp0_mod,m,m));
shading interp
view(-30,30)
zlim_yp = get(gca,'ZLim');
box on

fprintf('t = %.4f    [Press any key]\n',t0);
pause;

nout = 7;
tout = 0.01;

for iout = 1:nout
  [status,t,yy,yp] = IDASolve(tout,'Normal');
  tout = 2*tout;

  figure(fp);
  subplot(2,2,3);
  set(hs1,'FaceAlpha',0.15);
  hs1 = surf(reshape(yy,m,m));
  shading interp
  set(gca,'ZLim',zlim_yy);
  subplot(2,2,4);
  set(hs2,'FaceAlpha',0.15);
  hs2 = surf(reshape(yp,m,m));
  shading interp
  set(gca,'ZLim',zlim_yp);

  fprintf('t = %.4f    [Press any key]\n',t);
  pause;
  
end

IDAFree;

