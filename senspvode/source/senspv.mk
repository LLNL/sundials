# Makefile for SENS_PVODE.  Version of 25 August 2000.

# When run in the directory SensPVODE/source, this compiles
# all SensPVODE source files there, leaving .o files in the same directory,
# and creating a lib file libpvode.a in SensPVODE/lib.

# Notes on machine and environment dependencies:
#
# (1) The C compiler is given by CC below.
#
# (2) The MPI header file mpi.h is assumed to be in the directory given by
# MPI_INC below.
#
# (3) Other environment-specific compiler flags should appear in CFLAGS below.
#
# Change these variables as needed for other environments.


# Compiler:

# Version for the CASC Suns:
CC          = mpicc


# MPI include file directory:

# Version for the CASC Sun cluster with MPI:
MPI_INC = /usr/local/mpi/mpich/include

# PVODE include directory:
PVINC = ../include

# PVODE library file
LIBF = ../lib/libpvode.a

# Object files
OBJS = cvode.o nvector.o llnlmath.o cvspgmr.o spgmr.o iterativ.o smalldense.o \
       senscvspgmr.o  senspvode.o

# Compiler flags:
CFLAGS = -I$(PVINC) -g


HDRS = 	$(PVINC)/llnltyps.h $(PVINC)/cvode.h $(PVINC)/iterativ.h \
	$(PVINC)/cvspgmr.h $(PVINC)/sens_nvector.h \
	$(PVINC)/llnlmath.h $(MPI_INC)/mpi.h \
	$(PVINC)/senscvspgmr.h $(PVINC)/senspvode.h 

default:
	cp sens_nvector.c nvector.c
	cp ../include/sens_nvector.h ../include/nvector.h
	$(MAKE) -f Makefile mklib

mklib: $(OBJS) $(HDRS)
	ar rcv $(LIBF) $(OBJS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -f *.o
