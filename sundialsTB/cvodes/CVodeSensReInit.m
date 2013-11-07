function status = CVodeSensReInit(yS0, options)
%CVodeSensReInit reinitializes CVODES's FSA-related memory
%   assuming it has already been allocated in prior calls to CVodeInit 
%   and CVodeSensInit.
%   The number of sensitivities Ns is assumed to be unchanged since the 
%   previous call to CVodeSensInit.
%
%   Usage: CVodeSensReInit ( YS0 [, OPTIONS ] ) 
%
%   YS0      Initial conditions for sensitivity variables.
%            YS0 must be a matrix with N rows and Ns columns, where N is the problem
%            dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the CVodeSensSetOptions function. 
%
%   See also: CVodeSensSetOptions, CVodeReInit, CVodeSensInit

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
% $Revision: 1.4 $Date: 2007/08/21 17:42:38 $

mode = 13;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end

status = cvm(mode, yS0, options);
