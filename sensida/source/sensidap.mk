# File SensMakefile.idap  Version of 26 July 2001

# Makefile for the parallel version of the SensIDA package

# Notes on machine and environment dependencies:
#
# (1) Files from the SensIDA solver are accessed as follows:
#     (a) All SensIDA header files needed are assumed to be in the
#         directory given by SensIDAINC below.
#     (b) The SensIDA library file is given by SensIDALIB below.
#
# (2) The C compiler is given by CC below.
#
# (3) Files for the MPI library are accessed as follows:
#     (a) The header file mpi.h, if needed, is assumed to be in the
#         directory given by MPI_INC below.
#
# Change these variables as needed for other environments.


SensIDAINC = ../include
SensIDALIB = ../lib/libidap.a


# SUN Sparc version:
# This version uses the MPICH version of the MPI library through mpicc.
CC          = mpicc

# Version for the LLNL Cray-T3D with the EPCC version of MPI.
#CC        = /mpp/bin/cc


# Source files:

CODES = ida.c llnlmath.c idaspgmr.c spgmr.c iterativ.c smalldense.c \
        band.c nvector.c idabbdpre.c sensida.c sensidaspgmr.c sensidaband.c \
	sensidadense.c

# Object files:

OBJS =  ida.o llnlmath.o idaspgmr.o spgmr.o iterativ.o smalldense.o \
        band.o nvector.o idabbdpre.o sensida.o sensidaspgmr.o sensidaband.o \
	sensidadense.o


# Command sequence:
# Note that the proper versions of the files nvector.c and nvector.h are
# copied into the proper directories for the compile, then removed.
# However, leave the appropriate version of nvector.h in /include .
# This is needed by the driver programs in /example .

all:
	cp ./parallel/sens_nvector.c nvector.c
	cp $(SensIDAINC)/parallel/sens_nvector.h $(SensIDAINC)/nvector.h
	cp $(SensIDAINC)/parallel/sens_nvector.h $(SensIDAINC)/parallel/nvector.h
	$(CC) -g -c $(CODES) -I$(SensIDAINC)
	ar rcv $(SensIDALIB) $(OBJS)
	rm -f *.o
#	rm $(SensIDAINC)/nvector.h
#	rm ./nvector.c
