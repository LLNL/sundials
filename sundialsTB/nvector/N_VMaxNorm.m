function ret = N_VMaxNorm(x, comm)
%N_VMaxNorm returns the L-infinity norm of x
%
%   Usage:  RET = N_VMaxNorm ( X [, COMM] )
%
%If COMM is not present, N_VMaxNorm returns the L-infinity norm 
%of the local portion of X. Otherwise, it returns the global
%L-infinity norm..

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

if nargin == 1
  
  ret = norm(x,'inf');
  
else
  
  lnrm = norm(x,'inf');
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'MAX',comm);
  ret = gnrm;
  
end