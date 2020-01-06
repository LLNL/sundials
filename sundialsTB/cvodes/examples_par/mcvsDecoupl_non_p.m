function [] = mcvsDecoupl_non_p(comm)
%mcvsDecoupl_non_p - CVODES example problem
%   (parallel, Adams, Functional)
%   This is a simple test for the CVODES solver. It solves a
%   set of decoupled ODEs.
%   
%   See also: mpirun

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 23:09:18 $

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

options = CVodeSetOptions('Reltol',rtol,'AbsTol',atol);

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

CVodeInit(@rhsfn,'Adams','Functional',t0,y0,options,data);

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

%
%---------------------------------------------------------
%

function [yd, flag, new_data] = rhsfn(t, y, data)

alpha  = data.alpha;
nlocal = data.nlocal;
mype   = data.mype;

for i = 1:nlocal
  yd(i) = -alpha * (mype*nlocal + i) * y(i);
end

flag = 0;
new_data = [];
