# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# File Makefile.fsenskinsolp  Version of 07 September 2000 
#
# Makefile for fsenskinsolp 
#
# Makefile for senskinsolp (MPI/FSENSKINSOL)
#
# Notes on machine and environment dependencies:
#
# (1) Files from the SensKINSOL solver are accessed as follows:
#     (a) All SensKINSOL header files needed are assumed to be in the
#         directory given by KNC below.
#     (b) The FSENSKINSOL library file libsenskinsolp.a is assumed to be
#	  given by KINLIB below.
#
# (2) Files for the MPI library are accessed as follows:
#         The header file mpi.h is assumed to be in the directory given
#         by MPI_INC below.
#
# (3) The C compiler is given by CC below.
#
# (4) Environment-specific compiler flags must be included in CFLAGS below.
# 
# Change these variables as needed for other environments.
# 
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

all: makelib

# Version for the CASC Suns:
CC = mpicc
MPI_INC = /usr/local/mpi/mpich/include

## Version for the LLNL Cray-T3D.
#CC = /mpp/bin/cc
#MPI_INC = /usr/include/mpp


# parallel or serial version
VECTOR = parallel 


KNC = ../include
KINLIB = ../lib/libsens_kinsolp.a
CFLAGS =  -I$(MPI_INC) -I$(KNC) -I$(KNC)/$(VECTOR)

# List of object files needed

OBJS = fsenskinsolp.o fsenskinspgmr01.o fsenskinspgmr10.o \
       fsenskinspgmr11.o fsenskinspgmr20.o fsenskinspgmr21.o \
       fkinsolp.o fkinspgmr01.o fkinspgmr10.o fkinspgmr11.o  \
       fkinspgmr20.o fkinspgmr21.o fkinpreco.o fkinpsol.o \
       fkinuatimes.o

makelib: clearvect  $(OBJS)
	ar rcv $(KINLIB) $(OBJS)
	rm *.o


.c.o:
	$(CC) $(CFLAGS)  -c $*.c


clearvect:
	@if test -f $(KNC)/nvector.h ; then rm $(KNC)/nvector.h ; fi
