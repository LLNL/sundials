# File fkinsols.mk.  Version of 13 December 2001
#--------------------------------------------------------------------------
# Makefile for the serial version of the FKINSOL package.
#
# This makefile compiles all sources required for the solver, adds to the
# existing library file, and removes the .o files.
#
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------


#==========================================================================
# Machine-dependent variables
#==========================================================================

#-----------
# C compiler
#-----------
CC = gcc


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
LIBFILE    = ../../lib/libkinsols.a

#--------------
# Include files
#--------------

INC_SHARED = $(SHARED_DIR)/include

INC_SOLVER = ../include

#-------------
# Source files
#-------------

SRCS_SOLVER = fkinsols.c fkinspgmr01.c fkinspgmr10.c fkinspgmr11.c \
              fkinspgmr20.c fkinspgmr21.c fkinpreco.c fkinpsol.c \
              fkinuatimes.c


#-------------
# Object files
#-------------
OBJ         = fkinsols.o fkinspgmr01.o fkinspgmr10.o fkinspgmr11.o \
              fkinspgmr20.o fkinspgmr21.o fkinpreco.o fkinpsol.o \
              fkinuatimes.o


#---------------
# Compiler flags
#---------------
CFLAGS      = -I$(INC_SOLVER) -I$(INC_SHARED)


#==========================================================================
# Make rules
#==========================================================================

all: 
	@(echo 'Compile sources...')
	@($(CC) $(CFLAGS) -c $(SRCS_SOLVER))
	@(echo 'Create archive...')
	@(ar rc $(LIBFILE) $(OBJ))
	@(echo 'Clean up...')
	@(rm -f *.o)


