function ret = N_VMax(x,comm)
%N_VMax returns the largest element of x
%
%   Usage:  RET = N_VMax ( X [, COMM] )
%
%If COMM is not present, N_VMax returns the maximum value of 
%the local portion of X. Otherwise, it returns the global
%maximum value.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

if nargin == 1
  
  ret = max(x);
  
else
  
  lmax = max(x);
  gmax = 0.0;
  MPI_Allreduce(lmax,gmax,'MAX',comm);
  ret = gmax;
  
end