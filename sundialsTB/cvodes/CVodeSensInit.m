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
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision$Date: 2007/08/21 17:42:38 $

mode = 3;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  options = [];
end

status = cvm(mode, Ns, fctS, yS0, options);
