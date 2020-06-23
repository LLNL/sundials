function status = CVodeSensInit(Ns,fctS,yS0,options)
%CVodeSensInit allocates and initializes memory for FSA with CVODES.
%
%   Usage: CVodeSensInit ( NS, SFUN, YS0 [, OPTIONS ] ) 
%
%   NS       is the number of parameters with respect to which sensitivities
%            are desired
%   SFUN     is a function defining the righ-hand sides of the sensitivity
%            ODEs yS' = fS(t,y,yS).
%   YS0      Initial conditions for sensitivity variables.
%            YS0 must be a matrix with N rows and Ns columns, where N is the problem
%            dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the CVodeSetFSAOptions function. 
%
%   See also CVodeSensSetOptions, CVodeInit, CVSensRhsFn
%

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:42:38 $

mode = 3;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  options = [];
end

status = cvm(mode, Ns, fctS, yS0, options);
