function ret = N_VDotProd(x,y,comm)
%N_VDotProd returns the dot product of two vectors
%
%   Usage:  RET = N_VDotProd ( X, Y [, COMM] )
%
%If COMM is not present, N_VDotProd returns the dot product of the
%local portions of X and Y. Otherwise, it returns the global dot
%product.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


if nargin == 2
  
  ret = dot(x,y);
  
else
  
  ldot = dot(x,y);
  gdot = 0.0;
  MPI_Allreduce(ldot,gdot,'SUM',comm);
  ret = gdot;
  
end