# Makefile for SENSFPVODE.  Version of 25 August 2000.

# When run in the directory SensPVODE/fcmix, this compiles the
# relevant source files there, for the Fortran/C interface packages for
# both SensPVODE and PVBBDPRE.  It leaves .o files in the same directory,
# and creates a lib file libfpv.a in the directory ../lib.

# Notes on machine and environment dependencies:
#
# (1) The C compiler is given by CC below.
#
# (2) The MPI header file mpi.h is assumed to be in the directory given by
# MPI_INC below.
#
# (3) The SensPVODE header files are assumed to be in the directory given
# by PVINC below.
#
# (4) The PVBBDPRE preconditioner module header and object files are assumed
# to be in the directory given by PVBBDPRE below.  (First run 'make' there
# if necessary.)
#
# (5) Other environment-specific compiler flags should appear in CFLAGS below.
#
# Change these variables as needed for other environments.


# Compiler:

# Version for the LLNL Cray-T3D:
#CC          = /mpp/bin/cc
# Version for the CASC Suns:
CC          = mpicc


# MPI include file directory:

# Version for the LLNL Cray-T3D with the EPCC version of MPI:
#MPI_INC = /usr/include/mpp
# Version for the CASC Sun cluster with MPI:
MPI_INC = /usr/local/mpi/mpich/include

# PVODE include directory:
PVINC = ../include

# Location of PVBBDPRE header and object files
PVBBDPRE = ../precon


# Compiler flags:
CFLAGS = -I$(MPI_INC) -I$(PVINC) -I$(PVBBDPRE)


# List of object files needed
OBJS = fpvode.o fpvspgmr1.o fpvspgmr2.o fpvbbd.o $(PVBBDPRE)/pvbbdpre.o $(PVBBDPRE)/band.o sensfpvode.o sensfpvspgmr1.o sensfpvspgmr2.o sensfpvbbd.o

mklib: $(OBJS)
	(ar rcv libfpv.a $(OBJS); mv libfpv.a ../lib)

fpvode.o: fpvode.c fpvode.h fcmixpar.h
	($(CC) $(CFLAGS) -c fpvode.c)

fpvspgmr1.o: fpvspgmr1.c fpvode.h fcmixpar.h
	($(CC) $(CFLAGS) -c fpvspgmr1.c)

fpvspgmr2.o: fpvspgmr2.c fpvode.h fcmixpar.h
	($(CC) $(CFLAGS) -c fpvspgmr2.c)

fpvbbd.o: fpvbbd.c fpvode.h fpvbbd.h fcmixpar.h
	($(CC) $(CFLAGS) -c fpvbbd.c)

sensfpvode.o: sensfpvode.c sensfpvode.h fcmixpar.h
	($(CC) $(CFLAGS) -c sensfpvode.c)

sensfpvspgmr1.o: sensfpvspgmr1.c 
	($(CC) $(CFLAGS) -c sensfpvspgmr1.c)

sensfpvspgmr2.o: sensfpvspgmr2.c 
	($(CC) $(CFLAGS) -c sensfpvspgmr2.c)

sensfpvbbd.o: sensfpvbbd.c sensfpvbbd.h fcmixpar.h
	($(CC) $(CFLAGS) -c sensfpvbbd.c)
