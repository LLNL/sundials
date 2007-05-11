function CVodeSensReInit(yS0, options)
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
%            the CVodeSetFSAOptions function. 
%
%   See also: CVodeSetFSAOptions, CVodeReInit, CVodeSensInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/11/25 19:57:25 $

mode = 13;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end

cvm(mode, yS0, options);
