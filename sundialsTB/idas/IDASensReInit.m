function status = IDASensReInit(yyS0,ypS0,options)
%IDASensReInit reinitializes IDAS's FSA-related memory
%   assuming it has already been allocated in prior calls to IDAInit
%   and IDASensInit.
%   The number of sensitivities Ns is assumed to be unchanged since the 
%   previous call to IDASensInit.
%
%   Usage: IDASensReInit ( YYS0, YPS0 [, OPTIONS ] ) 
%
%   YYS0, YPS0   Initial conditions for sensitivity variables.
%            YYS0 and YPS0 must be matrices with N rows and Ns columns, where N is 
%            the problem dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the IDASetFSAOptions function. 
%
%   See also: IDASensSetOptions, IDAReInit, IDASensInit

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
% $Revision$Date: 2007/08/21 17:38:43 $

mode = 13;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = idm(mode, yyS0, ypS0, options);
