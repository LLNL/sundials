function [yd, flag, new_data] = pvkx_f(t, y, data)
%PVKX_F - RHS function for the PVKX example problem.
%
%   see also: pvkx, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:46 $


% Do all inter-process communication
% After this, data.yext contains all data needed for Finite Differences
data = pvkx_comm(y, data);

% Compute right-hand side locally
[yd, flag, new_data] = pvkx_fl(t, y, data);

new_data = data;

%=================================================================

function [data] = pvkx_comm(y, data)
%PVKX_COMM loads the local extended 3D solution matrix, by:
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


%=================================================================

function b = get_slice(a, dim, indx)
%GET_SLICE extracts from the 3D matrix A, the 2D matrix slice B at 
%  index INDX in the dimension DIM

switch dim
 case 1
  b = a(indx,:,:);
 case 2
  b = a(:,indx,:);
 case 3
  b = a(:,:,indx);
end

%=================================================================

function a = set_slice(a, dim, indx, b)
%SET_SLICE loads the 2D matrix B at index INDX in the dimension DIM 
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
