function [] = mpirun(fct,npe,dbg)
%MPIRUN runs parallel examples.
%
%  Usage: MPIRUN ( FCT , NPE [, DBG] )
%
%  FCT - function to be executed on all MATLAB processes.
%  NPE - number of processes to be used (including the master).
%  DBG - flag for debugging [ true | {false} ]
%        If true, spawn MATLAB child processes with a visible xterm.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

ih = isa(fct,'function_handle');
is = isa(fct,'char');
if ih
  sh = functions(fct);
  fct_str = sh.function;
elseif is
  fct_str = fct;
else
  error('mpirun:: Unrecognized function');
end

if exist(fct_str) ~= 2
  err_msg = sprintf('mpirun:: Function %s not in search path.',fct_str);
  error(err_msg);
end
  
nslaves = npe-1;
mpistart(nslaves);

debug = false;
if (nargin > 2) & dbg
  debug = true;
end

cmd_slaves = sprintf('mpiruns(''%s'')',fct_str);

if debug
  cmd = 'xterm';
  args = {'-sb','-sl','5000','-e','matlab','-nosplash','-nojvm','-r',cmd_slaves};
else
  cmd = 'matlab';
  args = {'-nosplash','-nojvm','-r',cmd_slaves};
end

[info children errs] = MPI_Comm_spawn(cmd,args,nslaves,'NULL',0,'SELF');

[info NEWORLD] = MPI_Intercomm_merge(children,0);

nvm(1,NEWORLD);
feval(fct,NEWORLD);
nvm(2);

