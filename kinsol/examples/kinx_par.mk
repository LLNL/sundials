# File kinx_par.mk.  Version of 7 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel KINSOL examples
#
# For each example requested, this makefile compiles the source and
# generates the executable
#
# Usage:  make -f kinx_par.mk          [to display list of examples]
#         make -f kinx_par.mk <ex>     [where <ex> is kinwebp or kinwebbbd]
#         make -f kinx_par.mk examples [to make all examples]
#         make -f kinx_par.mk purge    [to remove all example executable files]
# 
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------

SHELL = /bin/sh

#--------------------------
# Top of SUNDIALS directory
#--------------------------
SUNDIALS_DIR = ../..

#---------------------
# Path to header files
#---------------------
INC_DIR	= $(SUNDIALS_DIR)/include

#----------------------
# Path to library files
#----------------------
LIB_DIR = $(SUNDIALS_DIR)/lib

#-------------
# Architecture
#-------------
ARCH = `uname -s`.`uname -m`

#=====================================
# Machine-dependent variables
#=====================================

#----------------------
# MPI include directory
#----------------------
INC_MPI	=	/usr/local/mpi/mpich/include

#-----------
# C compiler
#-----------

# CASC Sun cluster.
CC     = mpicc
CFLAGS = -I$(INC_DIR) -L$(LIB_DIR)

# COMPASS cluster.
#CC     = mpicc
#CFLAGS = -I$(INC_DIR) -L$(LIB_DIR)

# BLUE cluster
#CC     = mpcc
#CFLAGS = -I$(INC_DIR) -L$(LIB_DIR)

#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of parallel KINSOL examples (using the parallel NVECTOR module):')
	@(echo '  kinwebp:  2-D food web system, block-diagonal preconditioner')
	@(echo '  kinwebbbd: kinwebp problem with BBD preconditioner')

kinwebp: kinwebp.c
	@echo '...Compile kinwebp...'
	@$(CC) $(CFLAGS) -o kinwebp kinwebp.c -lkinsol.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f kinwebp.o

kinwebbbd: kinwebbbd.c
	@echo '...Compile kinwebbbd...'
	@$(CC) $(CFLAGS) -o kinwebbbd kinwebbbd.c -lkinsol.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f kinwebbbd.o

examples: kinwebp kinwebbbd

purge:
	@(rm -f kinwebp)
	@(rm -f kinwebbbd)

#---End of kinx_par.mk---