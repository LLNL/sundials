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
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.3 $Date: 2007/08/21 17:38:43 $

mode = 13;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = idm(mode, yyS0, ypS0, options);
