# File fpvode.mk.  Version of 28 December 2001
#--------------------------------------------------------------------------
# Makefile for FPVODE, the Fortran interface package for PVODE.
#
# This makefile compiles all the required sources, adds to the library
# file libpvode.a, and removes the .o files.
#
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------


#==========================================================================
# Machine-dependent variables
#==========================================================================

#----------------------
# MPI include directory
#----------------------
INC_MPI = /usr/local/mpi/mpich/include

#-----------
# C compiler
#-----------
# SUN cluster (SunOS)
CC = mpicc

# COMPASS cluster (OSF1)
#CC	= mpicc

# BLUE cluster (AIX)
#CC	= mpcc

# LLNL Cray-T3D with the EPCC version of MPI.
#CC	= /mpp/bin/cc


#==========================================================================
# Other variables
#==========================================================================

#-------------
# Library file
#-------------
LIBFILE = ../../lib/libpvode.a

#--------------
# Include files
#--------------

INC_SHARED = ../../shared/include

INC_SOLVER = ../include

#-------------
# Source files
#-------------
SRCS = fpvode.c fpvpreco.c fpvpsol.c fpvjtimes.c fpvspgmr01.c \
       fpvspgmr10.c fpvspgmr11.c fpvspgmr20.c fpvspgmr21.c \
       fpvbbd.c fpvbbdin1.c

#-------------
# Object files
#-------------
OBJ = fpvode.o fpvpreco.o fpvpsol.o fpvjtimes.o fpvspgmr01.o \
      fpvspgmr10.o fpvspgmr11.o fpvspgmr20.o fpvspgmr21.o \
      fpvbbd.o fpvbbdin1.o

#---------------
# Compiler flags
#---------------
CFLAGS = -I$(INC_MPI) -I$(INC_SOLVER) -I$(INC_SHARED)


#==========================================================================
# Make rules
#==========================================================================

all: 
	@(echo 'Compile sources...')
	@($(CC) $(CFLAGS) -c $(SRCS))
	@(echo 'Create archive...')
	@(ar rc $(LIBFILE) $(OBJ))
	@(echo 'Clean up...')
	@(rm -f *.o)
