function status = IDASensInit(Ns,fctS,yyS0,ypS0,options)
%IDASensInit allocates and initializes memory for FSA with IDAS.
%
%   Usage: IDASensInit ( NS, SFUN, YYS0, YPS0 [, OPTIONS ] ) 
%
%   NS       is the number of parameters with respect to which sensitivities
%            are desired
%   SFUN     is a function defining the residual of the sensitivity DAEs
%            fS(t,y,yp,yS,ypS).
%   YYS0, YPS0   Initial conditions for sensitivity variables.
%            YYS0 and YPS0 must be matrices with N rows and Ns columns, where N is 
%            the problem dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the IDASetFSAOptions function. 
%
%   See also IDASensSetOptions, IDAInit, IDASensResFn

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:38:42 $

mode = 3;

if nargin < 4
  error('Too few input arguments');
end

if nargin < 5
  options = [];
end

status = idm(mode, Ns, fctS, yyS0, ypS0, options);
