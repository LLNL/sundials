function [] = pvnx(comm)
%PVNX - CVODES example problem (parallel, Adams, Functional)
%   This is a simple test for the CVODES solver. It solves a
%   set of decoupled ODEs.
%   
%   See also: mpirun pvnx_f

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:46 $

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

fprintf('\nPVNX example problem\n\n');
fprintf('  Processor %d/%d\n',mype,npes);
fprintf('  Global problem size: %d\n',neq);
fprintf('  Local problem size:  %d\n\n',nlocal);
if mype == 0
  fprintf('  alpha = %f\n',alpha);
  fprintf('  rtol = %e  atol = %e\n\n',rtol,atol);
end

options = CVodeSetOptions('LMM','Adams',...
                          'NonlinearSolver','Functional',...
                          'Reltol',rtol,'AbsTol',atol);

mondata = struct;

if mype == 0
  mondata.mode = 'both';
  mondata.sol = true;
else
%  mondata.post = false;
  mondata.sol = true;
  mondata.cntr = false;
  mondata.stats = false;
end
options = CVodeSetOptions(options,...
                          'MonitorFn','CVodeMonitor',...
                          'MonitorData',mondata);

CVodeMalloc(@pvnx_f,t0,y0,options,data);


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

CVodeFree;
