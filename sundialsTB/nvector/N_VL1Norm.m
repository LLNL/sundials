function ret = N_VL1Norm(x,comm)
%N_VL1Norm returns the L1 norm of x
%
%   Usage:  RET = N_VL1Norm ( X [, COMM] )
%
%If COMM is not present, N_VL1Norm returns the L1 norm of 
%the local portion of X. Otherwise, it returns the global
%L1 norm..

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

if nargin == 1
  
  ret = norm(x,1);
  
else
  
  lnrm = norm(x,1);
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'MAX',comm);
  ret = gnrm;
  
end