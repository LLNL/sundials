#File Makefile.diagpf  Version of 6 July 2000

# Makefile for diagpf demo with MPI/KINSOL/FKINBBD/KINBBDPRE
#
# Notes on machine and environment dependencies:
#
# (1) Files from the KINSOL solver are accessed as follows:
#     (a) The KINSOL library file libkinsol.a is assumed to be
#	  given by KINLIB below.
#     (b) All KINSOL header files needed are assumed to be in the
#         directory given by KNC below.
#     (c) The FKINSOL files required for this Fortran interface example are
#	  assumed to be in KINLIB . This is done by running Makefile.fkinsolp .
# (2) Files for the MPI library are accessed as follows:
#     (a) The file libmpi.a is assumed to be in the directory given by
#	   MPI_LIB below.
#     (b) The header file mpi.h is assumed to be in the directory given
#         by MPI_INC below.
#     (c) The list of libraries required for MPI is given by MPI_LIBLIST below.
#
# (3) If any additional libraries other than libm.a are required,
# their abbreviated names should be added to LIB_LIST below.
#
# (4) The C compiler is given by CC below.
#
# Change these variables as needed for other environments.
# 
# The command line to execute this demo on N processors is as follows:
# On the Cray-T3D: diagscalep -npes N
# On the CASC Sun cluster:  mpirun -np N -machinefile <machine-list> diagpf
# where <machine-list> is a file with the N machine names.

VECTOR = parallel

# SUN Sparc version:
# This version uses the MPICH version of the MPI library.
MPI_LIB  = /usr/local/mpi/mpich/lib/solaris/ch_p4
MPI_INC = /usr/local/mpi/mpich/include
MPI_LIBLIST =
F77_LIBLIST =


# Version for the LLNL Cray-T3D with the EPCC version of MPI.
#MPI_LIB = /mpp/lib
#MPI_INC = /usr/include/mpp
#MPI_LIBLIST = -lmpi
#F77_LIBLIST = ????

#Version for the CASC Suns
CC          = mpicc
F77         = mpif77


#Version for the LLNL Cray-T3D.
#CC        = /mpp/bin/cc
#F77       = ????


#Version for the COMPASS DEC Alpha cluster.
#MPI_LIB = /usr/local/apps/mpi/lib/alpha/ch_shmem
#MPI_INC = /usr/local/apps/mpi/include
#MPI_LIBLIST = -lmpich -lrpc
#CC        = mpicc
#F77       = mpif77


LIB_LIST    = $(MPI_LIBLIST) $(F77_LIBLIST) -lm
KINLIB     = ../lib/libkinsolp.a
LIBS =  $(LIB_LIST)

KNC = ../include

diagpf: diagpf.o $(MPI_INC)/mpi.h $(MPI_INC)/mpif.h
	$(F77) -o diagpf diagpf.o $(KINLIB) $(LIBS)

diagpf.o: diagpf.f
	$(F77) -c -I$(MPI_INC) diagpf.f


