function [] = mpiruns(fct)
%MPIRUNS runs the parallel example on a child MATLAB process.
%
%  Usage: MPIRUNS ( FCT )
%
%  This function should not be called directly. It is called
%  by mpirun on the spawned child processes.

% Radu Serban <radu@llnl.gov>
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
% $Revision$Date: 2006/01/06 19:00:15 $

clc;

[dummy hostname]=system('hostname');
fprintf('mpiruns :: child MATLAB process on %s\n',hostname);

MPI_Init;

MPI_Errhandler_set('WORLD','RETURN');

[info parent] = MPI_Comm_get_parent;

fprintf('mpiruns :: waiting to merge MPI intercommunicators ... ');
[info NEWORLD] = MPI_Intercomm_merge(parent,1);
fprintf('OK!\n\n');

MPI_Errhandler_set(NEWORLD,'RETURN');

% Put the MPI communicator in the global workspace
global sundials_MPI_comm;
sundials_MPI_comm = NEWORLD;

% Get rank of current process and put it in the global workspace
[status mype] = MPI_Comm_rank(NEWORLD);
global sundials_MPI_rank;
sundials_MPI_rank = mype;

fprintf('mpiruns :: MPI rank: %d\n\n',mype);

fprintf('----------------------------------------------------------------\n\n');

% Call the user main program
feval(fct,NEWORLD);

% Clear the global MPI communicator variable
clear sundials_MPI_comm

% Finalize MPI on this slave
MPI_Finalize;