function [] = mcvsAdvDiff_FSA_non_p(comm)
%mcvsAdvDiff_FSA_non_p - CVODES forward sensitivity example
%                        (parallel, Adams, Functional)
%   Semi-discrete form of the advection-diffusion equation in 1-D:
%     du/dt = q1 * d^2 u / dx^2 + q2 * du/dx
%   on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
%   Homogeneous Dirichlet boundary conditions are posed, and the
%   initial condition is:
%     u(x,y,t=0) = x(2-x)exp(2x).
%   The PDE is discretized on a uniform grid of size MX+2 with
%   central differencing, and with boundary values eliminated,
%   leaving an ODE system of size NEQ = MX.
%
%   Optionally, sensitivities with respect to q1 and q2 are also computed. 
%
%   This program solves the problem with the option for nonstiff
%   systems: ADAMS method and functional iteration.
%   It uses scalar relative and absolute tolerances.
%   Output is printed at t = .5, 1.0, ..., 5.
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

sensi = true;

xmax = 2.0;
mx = 10;
dx = xmax/(mx+1);
neq = mx;

[status npes] = MPI_Comm_size(comm);
[status mype] = MPI_Comm_rank(comm);

nperpe = floor(neq/npes);
nrem = neq - npes*nperpe;

if mype < nrem
  nlocal = nperpe+1;
  my_base = mype * nlocal;
else
  nlocal = nperpe;
  my_base = mype * nperpe + nrem;
end

data.comm = comm;
data.nlocal = nlocal;
data.npes = npes;
data.mype = mype;
data.dx = dx;
data.p = [1.0 ; 0.5];

t0 = 0.0;
for i = 1:nlocal
  iglobal = my_base + i;
  x = iglobal * dx;
  u0(i,1) = x *(xmax-x)*exp(2.0*x);
end

rtol = 0.0;
atol = 1.0e-5;

options = CVodeSetOptions('Reltol',rtol,'AbsTol',atol);
CVodeInit(@rhsfn,'Adams','Functional',t0,u0,options,data);


if sensi

  Ns = 2;
  uS0 = zeros(nlocal,Ns);
  pbar = data.p;
  plist = [1;2];

  FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
                                   'ErrControl', 'on',...
                                   'ParamField', 'p',...
                                   'ParamList', plist,...
                                   'ParamScales', pbar);

  CVodeSensInit(Ns, [], uS0, FSAoptions);

end


if mype == 0
  fprintf('============================================================\n');
  fprintf('     T     Q       H      NST                    Max norm   \n');
  fprintf('============================================================\n');
end

nout = 10;
dtout = 0.5;
tout = dtout;
for i = 1:nout
  if sensi
    [status,t,u,uS] = CVode(tout,'Normal');
    PrintOutput(mype, comm, t, u, uS);
  else
    [status,t,u] = CVode(tout,'Normal');
    PrintOutput(mype, comm, t, u, []);
  end
  tout = tout + dtout;
end

CVodeFree;

%%
%%-------------------------------------------------------------------
%%

function [ud, flag, new_data] = rhsfn(t, u, data)

% Extract needed problem constants from data

dx = data.dx;
hordc = data.p(1) / dx^2;
horac = data.p(2) / (2.0*dx);

% Extract parameters for parallel computation

comm = data.comm;
npes = data.npes;
mype = data.mype;

nlocal = length(u);

% Compute related parameters

mype_m1 = mype-1;
mype_p1 = mype+1;
last_pe = npes-1;

% Local copy of state

y = [0.0 ; u ; 0.0];

% Pass needed data to processes before and after current one

if mype ~= 0
  MPI_Send(u(1), mype_m1, 0, comm);
end
if mype ~= last_pe
  MPI_Send(u(nlocal), mype_p1, 0, comm);
end

% Receive needed data from processes before and after current one

buf = 0.0;

if mype ~= 0
  MPI_Recv(buf, mype_m1, 0, comm);
  y(1) = buf;
else
  y(1) = 0.0; % zero BC
end

if mype ~= last_pe
  MPI_Recv(buf, mype_p1, 0, comm);
  y(nlocal+2) = buf;
else
  y(nlocal+2) = 0.0; % zero BC
end

for i = 2:nlocal+1
  ui = y(i);
  ul = y(i-1);
  ur = y(i+1);
  
  hdiff = hordc*(ul - 2.0*ui + ur);
  hadv  = horac * (ur-ul);
  
  ud(i-1) = hdiff + hadv;
  
end

flag = 0;
new_data = [];


%%
%%-------------------------------------------------------------------
%%

function [] = PrintOutput(mype, comm, t, u, uS)

umax = N_VMaxNorm(u,comm);

if ~isempty(uS)
  smax1 = N_VMaxNorm(uS(:,1),comm);
  smax2 = N_VMaxNorm(uS(:,2),comm);
end

if mype == 0
  si = CVodeGetStats;
  fprintf('%8.3e %2d  %8.3e %5ld\n', t,si.qlast,si.hlast,si.nst);
  fprintf('                                Solution       ');
  fprintf('%12.4e \n', umax);
  if ~isempty(uS)
    fprintf('                                Sensitivity 1  ');
    fprintf('%12.4e \n', smax1);
    fprintf('                                Sensitivity 2  ');
    fprintf('%12.4e \n', smax2);
  end
end