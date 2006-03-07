function [] = pvkx(comm)
%PVKX - CVODES example problem (parallel, BDF, Newton, GMRES, BBdPre)
%   This example solves a 3D advection-diffusion PDE with a distributed
%   source to simulate atmospheric dispersion.
%   
%   PVKX uses the BBDPre preconditioner module in CVODES.
%
%   See also: mpirun pvkx_f pvkx_fl

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:46 $

%---------------------------------
% Domain definition
%  xmin - left boundary
%  xmax - right boundary
%  m    - number of intervals
%  np   - number of processes
%---------------------------------

xmin(1) = 0.0; xmax(1) = 20.0; m(1) = 20; np(1) = 2;
xmin(2) = 0.0; xmax(2) = 20.0; m(2) = 40; np(2) = 2;
xmin(3) = 0.0; xmax(3) = 20.0; m(3) = 20; np(3) = 1;

%---------------------------------
% Get MPI id and no. of processes
%---------------------------------

[status npes] = MPI_Comm_size(comm);
[status myId] = MPI_Comm_rank(comm);

if npes ~= prod(np)
  error('Wrong number of processes');
end

%---------------------------------
% Set-up problem data
%---------------------------------

data = SetData(comm, npes, myId, xmin, xmax, m, np);

%--------------------------------
% Problem dimensions
%--------------------------------

nlocal = prod(data.ml);
neq = prod(data.m);

fprintf('\nPVKX example problem\n\n');
fprintf('  Processor %d/%d\n',myId,npes);
fprintf('  Global problem size: %d\n',neq);
fprintf('  Local problem size:  %d\n\n',nlocal);


%--------------------------------
% Initial conditions
%--------------------------------

% Initial time
t0 = 0.0;

% Initial states (concentrations)
y0 = zeros(nlocal,1);

% TEST communication pattern
fprintf('Local data structure\n');
fprintf('  myId = %d\n',data.myId);
fprintf('  xmin = %g %g %g\n',data.xmin);
fprintf('  xmax = %g %g %g\n',data.xmax);
fprintf('  dx   = %g %g %g\n',data.dx);
fprintf('  start  = %3d %3d %3d\n',data.start);
fprintf('  m      = %3d %3d %3d\n',data.m);
fprintf('  ml     = %3d %3d %3d\n',data.ml);
fprintf('  |yext| = %3d %3d %3d\n\n',size(data.yext));

%--------------------------------
% CVODES setup
%--------------------------------

% Tolerances
options = CVodeSetOptions('RelTol',1.e-8,...
                          'AbsTol',1.e-6);

% Linear solver
mldq = data.ml(1)+1;
mudq = data.ml(1)+1;
mlkeep = 2;
mukeep = 2;

options = CVodeSetOptions(options,...
                          'LinearSolver','GMRES',...
                          'PrecType','Left',...
                          'PrecModule','BBDPre',...
                          'GlocalFn',@pvkx_fl,...
                          'LowerBwidthDQ',mldq,...
                          'UpperBwidthDQ',mudq,...
                          'LowerBwidth',mlkeep,...
                          'UpperBwidth',mukeep);

% Monitoring
mondata = struct;
if myId ==0
  mondata.mode = 'text';
else
  mondata.post = false;
end
options = CVodeSetOptions(options,...
                          'MonitorFn','CVodeMonitor',...
                          'MonitorData',mondata);

% Memory allocation and initialization
CVodeMalloc(@pvkx_f,t0,y0,options,data);

%--------------------------------
% CVODES solution
%--------------------------------

tf = 0.01;

[status,t,y] = CVode(tf,'Normal');

if myId == 0
  si = CVodeGetStats
end


CVodeFree;

%
% ===========================================================
%

function d = SetData(comm, npes, myId, xmin, xmax, m, np)

%---------------------------------
% MPI stuff
%---------------------------------

d.comm = comm;
d.myId = myId;
d.np = np;

%---------------------------------
% Domain boundaries
%---------------------------------

d.xmin = xmin;
d.xmax = xmax;

%--------------------------------------
% Diffusion coefficient
%--------------------------------------

d.Dc = 1.0;

%--------------------------------------
% Velocity parameters: Poiseuille flow
% across y direction, max. velocity=1
%    v(y) = Vc*(L-y)*(L+y)
%--------------------------------------

d.L = 0.5 * ( xmax(2) - xmin(2) );
d.Vc = 1.0/d.L^2;


%--------------------------------------
% Grid spacing and differential volume
%   d.m -> number of internal points
%--------------------------------------

d.dx = (d.xmax - d.xmin) ./ m;
d.m = m - [1 1 1];
d.dOmega = prod(d.dx);

%------------------------------------------------
% Partitioning
%   d.left  -> left neighbours
%   d.right -> right neighbours
%   d.start -> left border in global index space
%   d.ml    -> length of subdomain
%-----------------------------------------------

npd = floor(d.m ./ np);

% in x direction

test = mod( myId , np(1) );

d.left(1) = myId-1;
d.right(1) = myId+1;
d.start(1) = npd(1) * test + 1;
d.ml(1) = npd(1);

if test == 0
  d.left(1) = myId;
end

if test == np(1)-1
  d.right(1) = myId;
  d.ml(1) = d.m(1) - d.start(1) + 1;
end

% in y direction

test = mod( floor(myId/np(1)) , np(2) );

d.left(2) = myId - np(1);
d.right(2) = myId + np(1);
d.start(2) = npd(2) * test + 1;
d.ml(2) = npd(2);

if test == 0
  d.left(2) = myId;
end

if test == np(2)-1
  d.right(2) = myId;
  d.ml(2) = d.m(2) - d.start(2) + 1;
end

% in z direction

test = mod( floor(myId/np(1)/np(2))  , np(3) );

d.left(3) = myId - np(1)*np(2);
d.right(3) = myId + np(1)*np(2);
d.start(3) = npd(3) * test + 1;
d.ml(3) = npd(3);

if test == 0
  d.left(3) = myId;
end

if test == np(3)-1
  d.right(3) = myId;
  d.ml(3) = d.m(3) - d.start(3) +1;
end

%--------------------------------------
% Space for extended local solution 
% 3D matrix
%--------------------------------------

d.yext = zeros([d.ml(1)+2 d.ml(2)+2 d.ml(3)+2]);

%--------------------------------------
% Source parameters: Gaussians with 
%   - A: amplitude
%   - S: sigma^2
%   - X: position
%--------------------------------------

d.A1 = 1.0;
d.S1 = 1.7^2;
d.X1 = 4.0; d.Y1 = 8.0; d.Z1 = 8.0;

d.A2 = 0.8;   
d.S2 = 3.0^2;
d.X2 = 16.0; d.Y2 = 12.0; d.Z2 = 12.0;

d.GMIN = 1.0e-5;

A1 = 1.0;   
S1 = 1.7^2;
X1 = [  4.0  8.0  8.0];

A2 = 0.8;   
S2 = 3.0^2;
X2 = [ 16.0 12.0 12.0];

GMIN = 1.0e-5;



d.s = zeros(d.ml);
for i = 1:d.ml(1)
  for j = 1:d.ml(2)
    for k = 1:d.ml(3)
      x = d.xmin + (d.start + [i-2 j-2 k-2] ) .* d.dx;
      s = A1 * prod( exp( -(X1-x).^2 / S1 )  ) + ...
          A2 * prod( exp( -(X2-x).^2 / S2 )  ) ;
      if s < GMIN
        s = 0.0;
      end
      d.s(i,j,k) = s;
    end
  end
end
