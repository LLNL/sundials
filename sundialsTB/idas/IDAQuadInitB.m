function status = IDAQuadInitB(idxB, fctQB, yQB0, optionsB)
%IDAQuadInitB allocates and initializes memory for backward quadrature integration.
%
%   Usage: IDAQuadInitB ( IDXB, QBFUN, YQB0 [, OPTIONS ] ) 
%
%   IDXB     is the index of the backward problem, returned by
%            IDAInitB.
%   QBFUN    is a function defining the righ-hand sides of the
%            backward ODEs yQB' = fQB(t,y,yB).
%   YQB0     is the final conditions vector yQB(tB0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the IDASetQuadOptions function. 
%
%   See also: IDAInitB, IDASetQuadOptions, IDAQuadRhsFnB 
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
% $Revision: 1.2 $Date: 2007/08/21 17:38:42 $

mode = 6;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  optionsB = [];
end

idxB = idxB-1;
status = idm(mode, idxB, fctQB, yQB0, optionsB);
