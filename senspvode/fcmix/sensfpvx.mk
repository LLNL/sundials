# File Makefile for SENSFPVODE examples.   Version of 25 August 2000.

# Makefile for examples using SensPVODE: diagnf, diagkf, diagkbf,
# sdiagnf, sdiagkf, sdiagkbf.

# Type 'make <example-name>' to generate executable file. 

# Notes on machine and environment dependencies:
#
# (1) The SensPVODE library file libpvode.a is assumed to be in the
# directory given by PVLIB below.
#
# (2) The library file for the SensPVODE Fortran/C interface package is
# assumed to be in the directory given by FPVLIB below.
#
# (3) Files for the MPI library are accessed as follows:
#     (a) The appropriate MPI library file is assumed to be in the directory
#         given by MPI_LIB below.
#     (b) The list of libraries required for MPI is given by MPI_LIBLIST below.
#
# (4) Files for Fortran 77 are accessed as follows:
#     (a) The Fortran 77 compiler and linker are given by FC and FLINKER below.
#     (b) Any needed Fortran libraries assumed to be in the directory given by 
#         F77_LIB below.
#     (c) The list of Fortran libraries required is given by F77_LIBLIST below.
#
# (5) If any additional libraries other than libm.a are required,
# their abbreviated names should be added to LIB_LIST below.
#
# (6) Environment-specific compiler flags must be included in CFLAGS below.
#
# Change these variables as needed for other environments.

# The command line to execute the example <ex> on N processors is as follows:
# On the Cray-T3D: <ex> -npes N
# On the CASC Sun cluster: mpirun -np N -machinefile <machine-list> <ex>
#    where <machine-list> is a file with the N machine names.


# PVODE library directory:
PVLIB = ../../lib

# Location of interface package library file:
FPVLIB = ../../lib


# MPI library and include file directories, and library lists:

# Version for the LLNL Cray-T3D with the EPCC version of MPI:
#MPI_LIB = /mpp/lib
#MPI_INC = /usr/include/mpp
#MPI_LIBLIST = -lmpi
# Version for the CASC Sun cluster:
MPI_LIB = /usr/local/mpi/mpich/lib/solaris/ch_p4
MPI_INC = /usr/local/mpi/mpich/include
MPI_LIBLIST = 

# Fortran compiler/linker, library directory, and library list
# Version for the LLNL Cray-T3D with the EPCC version of MPI:
#FC          = cf77
#FLINKER     = $(FC)
#F77_LIB =
#F77_LIBLIST = 
# Version for the CASC Sun cluster:
FC          = mpif77
FLINKER     = $(FC)
F77_LIB = 
F77_LIBLIST = 


# Library names and library citations:
LIB_LIST    = -lfpv -lpvode $(MPI_LIBLIST) $(F77_LIBLIST) -lm
LIBS = -L$(FPVLIB) -L$(PVLIB) -L$(MPI_LIB) -L$(F77_LIB) $(LIB_LIST)

# Compiler flags:
CFLAGS = -I$(MPI_INC)


diagnf: diagnf.o
	$(FLINKER) -o diagnf diagnf.o $(LIBS)

diagkf: diagkf.o
	$(FLINKER) -o diagkf diagkf.o $(LIBS)

diagkbf: diagkbf.o
	$(FLINKER) -o diagkbf diagkbf.o $(LIBS)

sdiagnf: sdiagnf.o
	$(FLINKER) -o sdiagnf sdiagnf.o $(LIBS)

sdiagkf: sdiagkf.o
	$(FLINKER) -o sdiagkf sdiagkf.o $(LIBS)

sdiagkbf: sdiagkbf.o
	$(FLINKER) -o sdiagkbf sdiagkbf.o $(LIBS)


%.o : %.f
	$(FC) $(CFLAGS) -c $<

clean:
	rm -f *.o
