function [] = pvnx(comm)
%PVNX - CVODES example problem (parallel, Adams, Functional)
%   This is a simple test for the CVODES solver. It solves a
%   set of decoupled ODEs.
%   
%   See also: mpirun pvnx_f

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

[status npes] = MPI_Comm_size(comm);
[status mype] = MPI_Comm_rank(comm);

nlocal = 20;
neq = npes * nlocal;


alpha = 10.0/neq;

data.alpha = alpha;
data.comm = comm;
data.nlocal = nlocal;
data.mype = mype;

t0 = 0.0;
for i = 1:nlocal
  y0(i,1) = 1.0;
end


rtol = 1.0e-5;
atol = 1.0e-10;

options = CVodeSetOptions('Adams','on',...
                          'NonlinearSolver','Functional',...
                          'Reltol',rtol,'AbsTol',atol);

mondata = [];
if mype ~= 0
  mondata.stats = false;
  mondata.cntr = false;
end
options = CVodeSetOptions(options,...
                          'MonitorFn','CVodeMonitor',...
                          'MonitorData',mondata);

CVodeMalloc(@pvnx_f,t0,y0,options,data);

if mype == 0
  fprintf('NEQ = %d  AlPHA = %f\n',neq, alpha);
  fprintf('RTOL = %e  ATOL = %e\n',rtol,atol);
  fprintf('NPE = %d\n',npes);
end


nout = 10;
dtout = 0.1;
tout = dtout;
for i = 1:nout
  [status,t,y] = CVode(tout,'Normal');
  if mype == 0
    si = CVodeGetStats;
    fprintf('t = %f  nst = %d  nfe = %d\n',t,si.nst,si.nfe);
  end
  tout = tout + dtout;
end


disp('DONE')

CVodeFree;
