# File pvsmake  Version of 12 November 2001
#------------------------------------------------------
# Makefile for the parallel version of the PVODES package
#------------------------------------------------------

#======================================================
# Machine dependent variables
#======================================================

#----------------------
# MPI include directory
#----------------------
INC_MPI = 

#-----------
# C compiler
#-----------
# SUN cluster (SunOS)
CC = mpicc

# COMPASS cluster (OSF1)
#CC = mpicc

# BLUE cluster (AIX)
#CC = mpcc

# LLNL Cray-T3D with the EPCC version of MPI.
#CC = /mpp/bin/cc

#======================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

#-------------
# Library name
#-------------
PVODESLIB  = ../../lib/libpvodes.a

#--------------
# Include files
#--------------

INC_SHARED = $(SHARED_DIR)/include
INC_PVS    = ../include

#-------------
# Source files
#-------------

SRC_SHARED = $(SHARED_DIR)/source/band.c \
             $(SHARED_DIR)/source/smalldense.c \
             $(SHARED_DIR)/source/iterativ.c \
             $(SHARED_DIR)/source/llnlmath.c \
             $(SHARED_DIR)/source/nvector.c \
             $(SHARED_DIR)/source/spgmr.c

SRC_PVS = cvodes.c cvsband.c pvsbbdpre.c cvsdiag.c cvsspgmr.c

#-------------
# Object files
#-------------

OBJ     = cvodes.o cvsband.o pvsbbdpre.o cvsdiag.o cvsspgmr.o \
          band.o smalldense.o iterativ.o llnlmath.o nvector.o spgmr.o

#---------------
# Compiler flags
#---------------
CFLAGS = -I$(INC_PVS) -I$(INC_SHARED)

#======================================================

all: 
	@(echo 'Compile sources...')
	@($(CC) $(CFLAGS) -c $(SRC_PVS) $(SRC_SHARED))
	@(echo 'Create archive...')
	@(ar rc $(PVODESLIB) $(OBJ))
	@(echo 'Clean up...')
	@(rm -f *.o)

#======================================================
