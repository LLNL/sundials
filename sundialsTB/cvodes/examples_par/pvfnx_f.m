function [ud, flag, new_data] = pvfnx_f(t, u, data)
%PVFNX_F - RHS function for the PVFNX example problem
%
%   See also: pvfnx, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:46 $

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
