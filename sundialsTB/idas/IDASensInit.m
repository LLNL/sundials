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

mode = 3;

if nargin < 4
  error('Too few input arguments');
end

if nargin < 5
  options = [];
end

status = idm(mode, Ns, fctS, yyS0, ypS0, options);
