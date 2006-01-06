function ret = N_VMin(x,comm)
%N_VMin returns the smallest element of x
%   Usage:  RET = N_VMin ( X [, COMM] )
%
%If COMM is not present, N_VMin returns the minimum value of 
%the local portion of X. Otherwise, it returns the global
%minimum value.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

if nargin == 1
  
  ret = min(x);
  
else
  
  lmin = min(x);
  gmin = 0.0;
  MPI_Allreduce(lmin,gmin,'MIN',comm);
  ret = gmin;
  
end