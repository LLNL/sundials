function [] = IDASensMalloc(Ns,meth,yyS0,ypS0,varargin)
%IDASensMalloc allocates and initializes memory for FSA with IDAS.
%
%   Usage: IDASensMalloc ( NS, METH, YYS0, YPS0 [, OPTIONS ] ) 
%
%   NS       is the number of parameters with respect to which sensitivities
%            are desired
%   METHOD   FSA solution method [ 'Simultaneous' | {'Staggered'} ]
%            Specifies the FSA method for treating the nonlinear system solution for
%            sensitivity variables. In the simultaneous case, the nonlinear systems 
%            for states and all sensitivities are solved simultaneously. In the 
%            Staggered case, the nonlinear system for states is solved first and then
%            the nonlinear systems for all sensitivities are solved at the same time. 
%   YYS0, YPS0   Initial conditions for sensitivity variables.
%            YYS0 and YPS0 must be matrices with N rows and Ns columns, where N is 
%            the problem dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the IDASetFSAOptions function. 
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

mode = 2;

if nargin < 3
  disp('IDASensMalloc:: too few parameters');
  return
end

options = [];
if nargin > 3
  options = varargin{1};
end

idm(mode,Ns,meth,yyS0,ypS0,options);
