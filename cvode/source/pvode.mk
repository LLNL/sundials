# File pvode.mk.  Version of 14 December 2001
#--------------------------------------------------------------------------
# Makefile for the PVODE package (parallel)
#
# This makefile compiles all sources required for the solver, creates a
# library file, and removes the .o files.
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
#CC = /mpp/bin/cc


#==========================================================================
# Other variables
#==========================================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

#-------------
# Library name
#-------------
LIBFILE    = ../../lib/libpvode.a

#--------------
# Include files
#--------------

INC_SHARED = $(SHARED_DIR)/include

INC_SOLVER = ../include

#-------------
# Source files
#-------------

SRC_SHARED  = $(SHARED_DIR)/source

SRCS_SHARED = $(SRC_SHARED)/llnlmath.c $(SRC_SHARED)/nvector.c \
              $(SRC_SHARED)/spgmr.c $(SRC_SHARED)/iterativ.c \
              $(SRC_SHARED)/band.c $(SRC_SHARED)/smalldense.c

SRCS_SOLVER = cvode.c cvdiag.c cvspgmr.c pvbbdpre.c

#-------------
# Object files
#-------------
OBJ         = cvode.o cvdiag.o cvspgmr.o pvbbdpre.o llnlmath.o \
              nvector.o spgmr.o iterativ.o band.o smalldense.o

#---------------
# Compiler flags
#---------------
CFLAGS      = -I$(INC_MPI) -I$(INC_SOLVER) -I$(INC_SHARED)


#==========================================================================
# Make rules
#==========================================================================

all: 
	@(echo 'Compile sources...')
	@($(CC) $(CFLAGS) -c $(SRCS_SOLVER) $(SRCS_SHARED))
	@(echo 'Create archive...')
	@(ar rc $(LIBFILE) $(OBJ))
	@(echo 'Clean up...')
	@(rm -f *.o)
