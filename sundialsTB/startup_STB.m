% STARTUP_STB		PATH/environ script for sundialsTB
%			add to/replacement for startup.m

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

p = getenv('STB_HOME');
if isempty(p)
  p=[getenv('HOME') '/matlab/sundialsTB'];
end

q = [p '/cvodes'];
if ~exist(q, 'dir'), clear p q, error('SUNDIALS Toolbox not found (CVODES M files)'), end
addpath(q);
q = [p '/cvodes/cvm'];
if ~exist(q, 'dir'), clear p q, error('SUNDIALS Toolbox not found (CVODES MEX files)'), end
addpath(q);

q = [p '/kinsol'];
if ~exist(q, 'dir'), clear p q, error('SUNDIALS Toolbox not found (KINSOL M files)'), end
addpath(q);
q = [p '/kinsol/kim'];
if ~exist(q, 'dir'), clear p q, error('SUNDIALS Toolbox not found (KINSOL MEX files)'), end
addpath(q);

q = [p '/nvector'];
if ~exist(q, 'dir'), clear p q, error('SUNDIALS Toolbox not found (NVECTOR M files)'), end
addpath(q);
q = [p '/nvector/nvm'];
if ~exist(q, 'dir'), clear p q, error('SUNDIALS Toolbox not found (NVECTOR MEX files)'), end
addpath(q);

if ~isempty(getenv('MPITB_ROOT'))

  q = [p '/putils'];
  if ~exist(q, 'dir'), clear p q, error('SUNDIALS Toolbox not found (PUTILS M files)'), end
  addpath(q);
  
end

clear p q

