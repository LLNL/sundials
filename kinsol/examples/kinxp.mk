# File kinxp.mk.  Version of 13 December 2001
#--------------------------------------------------------------------------
# Makefile for the parallel KINSOL examples
#
# For each example, this makefile compiles the source, links,
# and removes the .o files.
#
# Usage:  make -f kinxp.mk   [to display list of examples]
#         make -f kinxp.mk <ex>   [where <ex> is kinwebp or kinwebbbd]
#         make -f kinxp.mk examples   [to make all examples]
#         make -f kinxp.mk purge   [to remove all example executable files]
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

# CASC Sun cluster.
CC      = mpicc
CLINKER = $(CC)

# COMPASS cluster.
#CC      = mpicc
#CLINKER = $(CC)

# BLUE cluster
#CC      = mpcc
#CLINKER = $(CC)


#======================================================
# Other variables
#======================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

#--------------
# Include files
#--------------
INC_SHARED = $(SHARED_DIR)/include
INC_SOLVER = ../include

#---------------
# Compiler flags
#---------------
CFLAGS     = -I$(INC_MPI) -I$(INC_SOLVER) -I$(INC_SHARED)

#----------
# Libraries
#----------
LIB        = ../../lib/
LIB_LIST   = -lkinsolp -lm
LIBS       = -L$(LIB) $(LIB_LIST)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of parallel KINSOL examples:')
	@(echo '  kinwebp:  2-D food web system, block-diagonal preconditioner')
	@(echo '  kinwebbbd: kinwebp problem with BBD preconditioner')

kinwebp: 
	$(CC) $(CFLAGS) -c kinwebp.c
	$(CLINKER) -o kinwebp kinwebp.o $(LIBS)
	@(rm -f kinwebp.o)

kinwebbbd: 
	$(CC) $(CFLAGS) -c kinwebbbd.c
	$(CLINKER) -o kinwebbbd kinwebbbd.o $(LIBS)
	@(rm -f kinwebbbd.o)

examples: kinwebp kinwebbbd

purge:
	@(rm -f kinwebp)
	@(rm -f kinwebbbd)
