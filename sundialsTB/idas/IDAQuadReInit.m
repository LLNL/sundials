function status = IDAQuadReInit(yQ0, options)
%IDAQuadReInit reinitializes IDAS's quadrature-related memory
%   assuming it has already been allocated in prior calls to IDAInit 
%   and IDAQuadInit.
%
%   Usage: IDAQuadReInit ( YQ0 [, OPTIONS ] ) 
%
%   YQ0      Initial conditions for quadrature variables yQ(t0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the IDASetQuadOptions function. 
%
%   See also: IDASetQuadOptions, IDAQuadInit

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
% $Revision$Date: 2007/08/21 17:38:42 $

mode = 12;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end
  
status = idm(mode, yQ0, options);
