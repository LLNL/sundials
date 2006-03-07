function [] = startup_STB(varargin)
% STARTUP_STB		PATH/environ script for sundialsTB
%			add to/replacement for startup.m

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/01/06 22:30:27 $

% Figure out the location of sundialsTB

if nargin > 0

% Get the location of sundialsTB from the input argument
  p = varargin{1};

else

% try first a local install
  p=[getenv('HOME') '/matlab/sundialsTB'];

% try next a system-wide install
  if ~exist(p, 'dir')
    p = [getenv('MATLAB') '/toolbox/sundialsTB'];
  end

end

if ~exist(p, 'dir')
  clear p
  warning('SUNDIALS Toolbox not found'); 
  return
end

% Add sundialsTB components to path

q = [p '/cvodes'];
if ~exist(q, 'dir')
  warning('SUNDIALS Toolbox not found (CVODES M files)');
end
addpath(q);

q = [p '/cvodes/cvm'];
if ~exist(q, 'dir')
  warning('SUNDIALS Toolbox not found (CVODES MEX files)');
end
addpath(q);

q = [p '/kinsol'];
if ~exist(q, 'dir')
  warning('SUNDIALS Toolbox not found (KINSOL M files)');
end
addpath(q);

q = [p '/kinsol/kim'];
if ~exist(q, 'dir')
  warning('SUNDIALS Toolbox not found (KINSOL MEX files)');
end
addpath(q);

q = [p '/nvector'];
if ~exist(q, 'dir')
  warning('SUNDIALS Toolbox not found (NVECTOR M files)');
end
addpath(q);

if ~isempty(getenv('MPITB_ROOT'))

  q = [p '/putils'];
  if ~exist(q, 'dir')
    warning('SUNDIALS Toolbox not found (PUTILS M files)');
  end
  addpath(q);
  
end

clear p q

