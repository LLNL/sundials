function [output, status] = CVodeGet(key, varargin)
%CVodeGet extracts data from the CVODES solver memory.
%
%   Usage: RET = CVodeGet ( KEY [, P1 [, P2] ... ]) 
%
%   CVodeGet returns internal CVODES information based on KEY. For some values
%   of KEY, additional arguments may be required and/or more than one output is
%   returned.
%
%   KEY is a string and should be one of:
%    o DerivSolution - Returns a vector containing the K-th order derivative
%       of the solution at time T. The time T and order K must be passed through 
%       the input arguments P1 and P2, respectively:
%       DKY = CVodeGet('DerivSolution', T, K)
%    o ErrorWeights - Returns a vector containing the current error weights.
%       EWT = CVodeGet('ErrorWeights')
%    o CheckPointsInfo - Returns an array of structures with check point information.
%       CK = CVodeGet('CheckPointInfo)

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.6 $Date: 2007/08/21 17:42:38 $

mode = 32;

if strcmp(key, 'DerivSolution')
  t = varargin{1};
  k = varargin{2};
  [output, status] = cvm(mode, 1, t, k);
elseif strcmp(key, 'ErrorWeights')
  [output, status] = cvm(mode, 2);
elseif strcmp(key, 'CheckPointsInfo')
  [output, status] = cvm(mode, 4);
else
  error('CVodeGet:: Unrecognized key');
end