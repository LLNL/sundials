function ret = N_VDotProd(x,y,comm)
%N_VDotProd returns the dot product of two vectors
%
%   Usage:  RET = N_VDotProd ( X, Y [, COMM] )
%
%If COMM is not present, N_VDotProd returns the dot product of the
%local portions of X and Y. Otherwise, it returns the global dot
%product.

% Radu Serban <radu@llnl.gov>
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.1 $Date$


if nargin == 2
  
  ret = dot(x,y);
  
else
  
  ldot = dot(x,y);
  gdot = 0.0;
  MPI_Allreduce(ldot,gdot,'SUM',comm);
  ret = gdot;
  
end