#File Makefile.senskinsolp  Version of 21 August 2000

# Makefile for parallel sens_kinsol (MPI/Sens_KINSOL)
#
# Notes on machine and environment dependencies:
#
# (1) Files from the KINSOL solver are accessed as follows:
#     (a) The KINSOL library file is assumed to be
#	  given by KINLIB below.
#     (b) All KINSOL header files needed are assumed to be in the
#         directory given by KININC below.
#
# (2) Files for the MPI library are accessed as follows:
#     (a) The header file mpi.h is assumed to be in the directory given
#         by MPI_INC below.
#
# (3) The C compiler and linker are given by CC and CLINKER below.
#
# 
# Change these variables as needed for other environments.
# 
# SUN Sparc version:
# This version uses the MPICH version of the MPI library through mpicc.
CC          = mpicc 

## Version for the LLNL Cray-T3D with the EPCC version of MPI.
#CC        = /mpp/bin/cc ii$(/usr/include/mpp)
#

#parallel case
KININC     = ../include
KINLIB  = ../lib/libsens_kinsolp.a

#*******************************************************
#definitions

# Code files

CODES = sens_kinsol.c sens_kinspgmr.c kinsol.c llnlmath.c kinspgmr.c \
	spgmr.c iterativ.c smalldense.c kinbbdpre.c band.c nvector.c

# Object files

OBJS = sens_kinsol.o sens_kinspgmr.o kinsol.o llnlmath.o kinspgmr.o \
	spgmr.o iterativ.o smalldense.o kinbbdpre.o band.o nvector.o

all:
	cp ./parallel/nvector.c nvector.c
	cp ../include/parallel/nvector.h ../include/nvector.h
	$(CC) -c $(CODES) -I$(KININC)
	ar rcv $(KINLIB) $(OBJS)
	rm -f *.o

