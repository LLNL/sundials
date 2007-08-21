function [] = IDASensInit(Ns,fctS,yyS0,ypS0,options)
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
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $

mode = 3;

if nargin < 4
  error('Too few input arguments');
end

if nargin < 5
  options = [];
end

idm(mode, Ns, fctS, yyS0, ypS0, options);
