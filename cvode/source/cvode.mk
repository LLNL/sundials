# File cvode.mk.  Version of 14 December 2001
#--------------------------------------------------------------------------
# Makefile for the CVODE package (serial)
#
# This makefile compiles all sources required for the solver, creates a
# library file, and removes the .o files.
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
LIBFILE    = ../../lib/libcvode.a

#--------------
# Include files
#--------------

INC_SHARED = $(SHARED_DIR)/include

INC_SOLVER = ../include

#-------------
# Source files
#-------------

SRC_SHARED = $(SHARED_DIR)/source

SRCS_SHARED = $(SRC_SHARED)/llnlmath.c \
              $(SRC_SHARED)/spgmr.c \
              $(SRC_SHARED)/iterativ.c \
              $(SRC_SHARED)/band.c \
              $(SRC_SHARED)/dense.c \
              $(SRC_SHARED)/smalldense.c \
              $(SRC_SHARED)/nvector.c

SRCS_SOLVER = cvode.c cvdense.c cvband.c cvdiag.c cvspgmr.c \
              cvbandpre.c

#-------------
# Object files
#-------------
OBJ         = cvode.o cvdense.o cvband.o cvdiag.o cvspgmr.o \
              spgmr.o iterativ.o band.o dense.o smalldense.o \
              cvbandpre.o llnlmath.o nvector.o

#---------------
# Compiler flags
#---------------
CFLAGS      = -I$(INC_SOLVER) -I$(INC_SHARED)


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
