%IDADENX - IDAS example problem (serial, dense)
%
%   See also: idadenx_f, idadenx_g, idadenx_J

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/15 19:31:26 $

m = 10;
N = m^2;
data.m = m;
data.N = N;
data.dx = 1.0/(m-1);
data.c = 1.0/data.dx^2;


[t0,yy0,yp0,id,cnstr] = idabanx_ic(data); 

options = IDASetOptions('RelTol',0.0,...
                        'AbsTol',1.0e-3,...
                        'ICcalculation','FindAlgebraic',...
                        'VariableTypes',id,...
                        'ConstraintTypes',cnstr,...
                        'LinearSolver','Band',...
                        'LowerBwidth',m,...
                        'UpperBwidth',m);

mondata.mode = 'text';
mondata.update = 100;
options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);

IDAMalloc(@idabanx_f,t0,yy0,yp0,options,data);


fprintf('\n   Output Summary (umax = max-norm of solution) \n\n');
fprintf('  time       ymax     k  nst  nni  nje   nre   nreLS    h      \n' );
fprintf(' ------------------------------------------------------------- \n');

nout = 11;
tout = 0.01;

for iout = 1:nout
  [status,t,yy,yp] = IDASolve(tout,'Normal');
  si = IDAGetStats;
  ymax = N_VMaxNorm(yy);
  fprintf(' %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e \n',...
          t, ymax, si.qlast, si.nst, si.nni, si.LSInfo.njeB, si.nre, si.LSInfo.nreB, si.hlast);

  tout = 2*tout;
end

IDAFree;

