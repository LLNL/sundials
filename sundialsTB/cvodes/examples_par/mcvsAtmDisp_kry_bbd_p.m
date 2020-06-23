function [] = mcvsAtmDisp_kry_bbd_p(comm)
%mcvsAtmDisp_kry_bbd_p - CVODES example problem
%   (parallel, BDF, Newton, GMRES, BBDP)
%   This example solves a 3D advection-diffusion PDE with a
%   distributed source to simulate atmospheric dispersion.
%   
%   This example uses the BBDP preconditioner module in CVODES.
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
CVodeInit(@rhsfn,'BDF','Newton',t0,y0,options,data);

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

% ===========================================================

function [yd, flag, new_data] = rhsfn(t, y, data)

% Do all inter-process communication
% After this, data.yext contains all data needed for Finite Differences
data = rhsfn_comm(y, data);

% Compute right-hand side locally
[yd, flag, new_data] = rhsfn_local(t, y, data);

new_data = data;

%=================================================================

function [data] = rhsfn_comm(y, data)
%rhsfn_comm loads the local extended 3D solution matrix, by:
%   a) using the local solution in the interior of the local grid
%   b) communicating with neighbouring processes to obtain solution
%      on the internal boundaries
%   c) setting 0-flux B.C. on the external boundaries.

ml = data.ml;
comm = data.comm;
myId = data.myId;
left = data.left;
right = data.right;

% Reshape local solution into a cube
c = reshape(y,ml);

% Load local solution into extended matrix
data.yext(2:end-1,2:end-1,2:end-1) = c;

% Internal boundaries: loop over each dimension and exchange data.
%    The processor with lower ID always sends first.
% External boundaries: impose homogeneous Neumann B.C.

for dim = 1:3

  N = prod(ml)/ml(dim);     % dimension of communication buffers

% to the left

  nbr = left(dim);                                 % left neighbour

  bufS = reshape( get_slice(c,dim,1), N, 1);       % send buffer
  bufR = zeros(N, 1);                              % receive buffer

  if nbr == myId  % external BC
    data.yext = set_slice(data.yext,dim,1,bufS);
  else            % internal BC
    if myId < nbr
%      fprintf('  left send/recv  %d  N = %d\n',nbr,N);
      info = MPI_Send(bufS, nbr, 0, comm);
      if info ~= 0
        fprintf('Send to left: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
      [info, stat] = MPI_Recv(bufR, nbr, 0, comm);
      if info ~= 0
        fprintf('Receive from left: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
    else
%      fprintf('  left recv/send  %d  N = %d\n',nbr,N);
      [info, stat] = MPI_Recv(bufR, nbr, 0, comm);      
      if info ~= 0
        fprintf('Receive from left: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
      info = MPI_Send(bufS, nbr, 0, comm);
      if info ~= 0
        fprintf('Send to left: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
    end
    data.yext = set_slice(data.yext,dim,1,bufR);
  end
  
% to the right

  nbr = right(dim);                                % right neighbour

  bufS = reshape( get_slice(c,dim,ml(dim)), N, 1); % send buffer
  bufR = zeros(N, 1);                              % receive buffer

  if nbr == myId  % external BC
    data.yext = set_slice(data.yext,dim,ml(dim)+2,bufS);
  else            % internal BC
    if myId < nbr
%      fprintf('  right send/recv  %d  N = %d\n',nbr,N);
      info = MPI_Send(bufS, nbr, 0, comm);
      if info ~= 0
        fprintf('Send to right: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
      [info, stat] = MPI_Recv(bufR, nbr, 0, comm);
      if info ~= 0
        fprintf('Receive from right: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
    else
%      fprintf('  right recv/send  %d  N = %d\n',nbr,N);
      [info, stat] = MPI_Recv(bufR, nbr, 0, comm);      
      if info ~= 0
        fprintf('Receive from right: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
      info = MPI_Send(bufS, nbr, 0, comm);
      if info ~= 0
        fprintf('Send to right: myId %d  nbr %d  info %d',myId,nbr,info); 
        MPI_Abort(comm, 0);
      end
    end
    data.yext = set_slice(data.yext,dim,ml(dim)+2,bufR);
  end
  
  
end

function b = get_slice(a, dim, indx)
%get_slice extracts from the 3D matrix A, the 2D matrix slice B at 
%  index INDX in the dimension DIM

switch dim
 case 1
  b = a(indx,:,:);
 case 2
  b = a(:,indx,:);
 case 3
  b = a(:,:,indx);
end

function a = set_slice(a, dim, indx, b)
%set_slice loads the 2D matrix B at index INDX in the dimension DIM 
%  into the 3D matrix A. A has 2 more components than B in each 
%  dimension

[nr, nc] = size(b);   % number of rows and columns in B

switch dim
 case 1
  nr = size(a,2)-2;
  nc = size(a,3)-2;
  a(indx,2:end-1,2:end-1) = reshape(b,nr,nc);
 case 2
  nr = size(a,1)-2;
  nc = size(a,3)-2;
  a(2:end-1,indx,2:end-1) = reshape(b,nr,nc);
 case 3
  nr = size(a,1)-2;
  nc = size(a,2)-2;
  a(2:end-1,2:end-1,indx) = reshape(b,nr,nc);
end

% ===========================================================

function [yd, flag, new_data] = rhsfn_local(t, y, data)
%rhsfn_local - local RHS computation

xmin  = data.xmin;
ml    = data.ml;
start = data.start;
dx    = data.dx;
Dc    = data.Dc;
yext  = data.yext;

for i = 2:ml(1)+1
  for j = 2:ml(2)+1
    for k = 2:ml(3)+1

      x = xmin + (start + [i-2 j-2 k-2] ) .* dx;
      v = velocity(x, data);
      s = source(x,data);

      [c, cl, cr] = stencil(yext,i,j,k);

      adv = v .* (cr-cl) ./ (2.0*dx);
      dif = Dc * (cr - 2.0*c + cl) / dx.^2;
    
      yd(i-1,j-1,k-1) = s + sum(dif-adv);
      
    end
  end
end

yd = reshape(yd,prod(ml),1);

flag = 0;
new_data = [];


function [c,cl,cr] = stencil(yext,i,j,k)

c = yext(i,j,k) * ones(1,3);
cl(1) = yext(i-1,j,  k  ); cr(1) = yext(i+1,j,  k  );
cl(2) = yext(i,  j-1,k  ); cr(2) = yext(i,  j+1,k  );
cl(3) = yext(i,  j,  k-1); cr(3) = yext(i,  j,  k+1); 


function v = velocity(x, data)

L = data.L;
Vc = data.Vc;
xmin = data.xmin;

y = x(2) - xmin(2) - L;

v(1) = Vc * (L+y) * (L-y);
v(2) = 0.0;
v(3) = 0.0;


function s = source(x, data)

A1 = data.A1;  A2 = data.A2;
S1 = data.S1;  S2 = data.S2;
X1 = data.X1;  X2 = data.X2;
Y1 = data.Y1;  Y2 = data.Y2;
Z1 = data.Z1;  Z2 = data.Z2;

s1 = A1 * exp(-(X1-x(1))^2/S1) * exp(-(Y1-x(2))^2/S1) * exp(-(Z1-x(3))^2/S1);
s2 = A2 * exp(-(X2-x(1))^2/S2) * exp(-(Y2-x(2))^2/S2) * exp(-(Z2-x(3))^2/S2);

s = s1 + s2;

if s < data.GMIN
  s = 0.0;
end
