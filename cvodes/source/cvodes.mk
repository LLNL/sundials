# File cvodes.mk  Version of 6 February 2002
#------------------------------------------------------
# Makefile for the serial version of the CVODES package
#------------------------------------------------------

#======================================================
# Machine dependent variables
#======================================================

#-----------
# C compiler
#-----------
CC = gcc

#======================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

#-------------
# Library name
#-------------
CVODESLIB = ../../lib/libcvodes.a

#--------------
# Include files
#--------------
INC_SHARED = $(SHARED_DIR)/include
INC_CVS    = ../include

#--------------
# Source files
#--------------
SRC_SHARED = $(SHARED_DIR)/source/band.c \
             $(SHARED_DIR)/source/dense.c \
             $(SHARED_DIR)/source/smalldense.c \
             $(SHARED_DIR)/source/iterativ.c \
             $(SHARED_DIR)/source/llnlmath.c \
             $(SHARED_DIR)/source/nvector.c \
             $(SHARED_DIR)/source/spgmr.c

SRC_CVS    = cvodes.c cvsband.c cvsbandpre.c \
             cvsdense.c cvsdiag.c cvsspgmr.c

#-------------
# Object files
#-------------
OBJ        = cvodes.o cvsband.o cvsbandpre.o \
             cvsdense.o cvsdiag.o cvsspgmr.o \
             band.o dense.o smalldense.o iterativ.o \
             llnlmath.o nvector.o spgmr.o

#---------------
# Compiler flags
#---------------
#CFLAGS = -O3 -ffloat-store -I$(INC_CVS) -I$(INC_SHARED)
CFLAGS = -ffloat-store -Wall -I$(INC_CVS) -I$(INC_SHARED)
#CFLAGS = -I$(INC_CVS) -I$(INC_SHARED)

#======================================================

all: 
	@(echo 'Compile sources...')
	@($(CC) $(CFLAGS) -c $(SRC_CVS) $(SRC_SHARED))
	@(echo 'Create archive...')
	@(ar rc $(CVODESLIB) $(OBJ))
	@(echo 'Clean up...')
	@(rm -f *.o)

#======================================================
