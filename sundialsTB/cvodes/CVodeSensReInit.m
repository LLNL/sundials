function [] = CVodeSensReInit(meth,yS0,varargin)
%CVodeSensReInit reinitializes CVODES's FSA-related memory
%   assuming it has already been allocated in prior calls to CVodeMalloc 
%   and CVodeSensMalloc.
%   The number of sensitivities Ns is assumed to be unchanged since the 
%   previous call to CVodeSensMalloc.
%
%   Usage: CVodeSensReInit ( METH, YS0 [, OPTIONS ] ) 
%
%   METHOD   FSA solution method [ 'Simultaneous' | {'Staggered'} ]
%            Specifies the FSA method for treating the nonlinear system solution for
%            sensitivity variables. In the simultaneous case, the nonlinear systems 
%            for states and all sensitivities are solved simultaneously. In the 
%            Staggered case, the nonlinear system for states is solved first and then
%            the nonlinear systems for all sensitivities are solved at the same time. 
%   YS0      Initial conditions for sensitivity variables.
%            YS0 must be a matrix with N rows and Ns columns, where N is the problem
%            dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the CVodeSetFSAOptions function. 
%
%   See also: CVodeSensMalloc

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

mode = 12;

if nargin < 2
  disp('CVodeSensReInit:: too few parameters');
  return
end

options = [];
if nargin > 2
  options = varargin{1};
end

cvm(mode,0,meth,yS0,options);
