function [] = mpiruns(fct)
%MPIRUNS runs the parallel example on a child MATLAB process.
%
%  Usage: MPIRUNS ( FCT )
%
%  This function should not be called directly. It is called
%  by mpirun on the spawned child processes.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

[dum hostname]=system('hostname');
fprintf('child MATLAB process on %s\n',hostname);

MPI_Init;

MPI_Errhandler_set('WORLD','RETURN');

[info parent] = MPI_Comm_get_parent;

fprintf('waiting for the parent to merge MPI intercommunicators ... ');
[info NEWORLD] = MPI_Intercomm_merge(parent,1);
fprintf('OK!\n');

MPI_Errhandler_set(NEWORLD,'RETURN');

nvm(1,NEWORLD);
feval(fct,NEWORLD);
nvm(2);

MPI_Finalize;