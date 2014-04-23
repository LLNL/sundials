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
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision$Date: 2006/01/06 19:00:15 $

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

% Put the MPI communicator in the global workspace
global sundials_MPI_comm;
sundials_MPI_comm = NEWORLD;

% Get rank of current process and put it in the global workspace
[status mype] = MPI_Comm_rank(NEWORLD);
global sundials_MPI_rank;
sundials_MPI_rank = mype;

% Call the user main program
feval(fct,NEWORLD);

% Clear the global MPI communicator variable
clear sundials_MPI_comm

