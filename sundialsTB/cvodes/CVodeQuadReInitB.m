function status = CVodeQuadReInitB(idxB, yQB0, optionsB)
%CVodeQuadReInitB reinitializes memory for backward quadrature integration.
%
%   Usage: CVodeQuadReInitB ( IDXB, YS0 [, OPTIONS ] ) 
%
%   IDXB     is the index of the backward problem, returned by
%            CVodeInitB.
%   YQB0     is the final conditions vector yQB(tB0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the CVodeSetQuadOptions function. 
%
%   See also: CVodeSetQuadOptions, CVodeReInitB, CVodeQuadInitB
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
% $Revision: 1.2 $Date: 2007/05/11 18:51:32 $

mode = 16;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  optionsB = [];
end
  
idxB = idxB-1;
status = cvm(mode, idxB, yQB0, optionsB);
