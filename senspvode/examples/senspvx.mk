# File Makefile for SensPVODE examples.   Version of 25 August 2000.

# Makefile for the examples: pvnx.c, pvkx.c pvkxb.c, spvnx.c, spvkx.c,
#    and spvkxb.c.
# Type 'make pvnx', 'make pvkx', 'make pvkxb', 'make spvnx', 'make spvkx',
#    'make spvkxb' to generate executable file. 
#
# Notes on machine and environment dependencies:
#
# (1) Files for the MPI library are accessed as follows:
#     (a) The appropriate MPI library file is assumed to be in the directory
#         given by MPI_LIB below.
#     (b) The header file mpi.h is assumed to be in the directory given
#         by MPI_INC below.
#     (c) The list of libraries required for MPI is given by MPI_LIBLIST below.
#
# (2) Files from the SensPVODE solver are accessed as follows:
#     (a) The SensPVODE library files libpvode.a and libpvpre.a are assumed 
#         to be in the directory given by PVLIB below.
#     (b) All SensPVODE header files needed are assumed to be in the 
#         directory given by PVINC below.
#     (c) The header file for preconditioners is assumed to be
#         in the directory given by PVPRECON below.
#
# (3) The C compiler and linker are given by CC and CLINKER below.
#
# (4) If any additional libraries other than libm.a are required,
# their abbreviated names should be added to LIB_LIST below.
#
# (5) Environment-specific compiler flags must be included in CFLAGS below.
#
# Change these variables as needed for other environments.

# The command line to execute the example <ex> on N processors is as follows:
# On the Cray-T3D: <ex> -npes N
# On the CASC Sun cluster: mpirun -np N -machinefile <machine-list> <ex>
#    where <machine-list> is a file with the N machine names.


# MPI library and include file directories, and MPI library list:

# Version for the LLNL Cray-T3D with the EPCC version of MPI.
#MPI_LIB = /mpp/lib
#MPI_INC = /usr/include/mpp
#MPI_LIBLIST = -lmpi
# Version for the CASC Sun cluster.
MPI_LIB = /usr/local/mpi/mpich/lib/solaris/ch_p4
MPI_INC = /usr/local/mpi/mpich/include
MPI_LIBLIST = 

# SensPVODE library, include, and preconditioner directories:
PVLIB = ../lib
PVINC = ../include
PVPRECON = ../precon


# Compiler and linker:

# Version for the LLNL Cray-T3D:
#CC          = /mpp/bin/cc
#CLINKER     = $(CC)
# Version for CASC Sun cluster:
CC          = mpicc
CLINKER     = $(CC)

# Library names and library link directives:
LIB_LIST    = -lpvode -lpvpre $(MPI_LIBLIST) -lm
LIBS = -L$(PVLIB) -L$(MPI_LIB) $(LIB_LIST)

# Compiler flags:
CFLAGS =  -I$(MPI_INC) -I$(PVINC) -I$(PVPRECON)


all: pvkx pvnx pvkxb spvkx spvnx spvkxb

pvkx: pvkx.o $(PVLIB)/libpvode.a
	$(CLINKER) -o pvkx pvkx.o $(LIBS)

pvnx: pvnx.o $(PVLIB)/libpvode.a
	$(CLINKER) -o pvnx pvnx.o $(LIBS)

pvkxb: pvkxb.o $(PVLIB)/libpvode.a $(PVLIB)/libpvpre.a
	$(CLINKER) -o pvkxb pvkxb.o $(LIBS)

spvkx: spvkx.o $(PVLIB)/libpvode.a
	$(CLINKER) -o spvkx spvkx.o $(LIBS)

spvnx: spvnx.o $(PVLIB)/libpvode.a
	$(CLINKER) -o spvnx spvnx.o $(LIBS)

spvkxb: spvkxb.o $(PVLIB)/libpvode.a $(PVLIB)/libpvpre.a
	$(CLINKER) -o spvkxb spvkxb.o $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o
