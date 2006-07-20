function [] = pvfnx(comm)
%PVFNX - CVODES forward sensitivity example (parallel, Adams, Functional)
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
%   See also: mpirun pvfnx_f

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/03/07 01:19:52 $

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

options = CVodeSetOptions('LMM','Adams',...
                          'NonlinearSolver','Functional',...
                          'Reltol',rtol,'AbsTol',atol);
CVodeMalloc(@pvfnx_f,t0,u0,options,data);


if sensi

  Ns = 2;
  uS0 = zeros(nlocal,Ns);
  pbar = data.p;
  plist = [1;2];

  FSAoptions = CVodeSetFSAOptions('SensErrControl', 'on',...
                                  'ParamField', 'p',...
                                  'ParamList', plist,...
                                  'ParamScales', pbar,...
                                  'SensDQtype', 'Centered',...
                                  'SensDQparam', 0.0);

  CVodeSensMalloc(Ns, 'Simultaneous', uS0, FSAoptions);

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