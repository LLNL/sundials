# File idaxp.mk.  Version of 17 December 2001
#--------------------------------------------------------------------------
# Makefile for the parallel IDA examples
#
# For each example, this makefile compiles the source, links,
# and removes the .o files.
#
# Usage:  make -f idaxp.mk   [to display list of examples]
#         make -f idaxp.mk <ex>   [where <ex> is iheatpk, iheatbbd, iwebpk,
#                                  or iwebbbd]
#         make -f idaxp.mk examples   [to make all examples]
#         make -f idaxp.mk purge   [to remove all example executable files]
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
#CC          = mpicc
#CLINKER     = $(CC)

# BLUE cluster
#CC          = mpcc
#CLINKER     = $(CC)


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
LIB	     = ../../lib/
LIB_LIST = -lidap -lm
LIBS     = -L$(LIB) $(LIB_LIST)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of parallel IDA examples:')
	@(echo '  iheatpk:  2-D heat equation, diagonal preconditioner')
	@(echo '  iheatbbd: 2-D heat equation, BBD preconditioner')
	@(echo '  iwebpk:   2-D food web, block-diagonal preconditioner')
	@(echo '  iwebbbd:  2-D food web, BBD preconditioner')

iheatpk:
	$(CC) $(CFLAGS) -c iheatpk.c
	$(CLINKER) -o iheatpk iheatpk.o $(LIBS)
	@(rm -f iheatpk.o)

iheatbbd:
	$(CC) $(CFLAGS) -c iheatbbd.c
	$(CLINKER) -o iheatbbd iheatbbd.o $(LIBS)
	@(rm -f iheatbbd.o)

iwebpk:
	$(CC) $(CFLAGS) -c iwebpk.c
	$(CLINKER) -o iwebpk iwebpk.o $(LIBS)
	@(rm -f iwebpk.o)

iwebbbd:
	$(CC) $(CFLAGS) -c iwebbbd.c
	$(CLINKER) -o iwebbbd iwebbbd.o $(LIBS)
	@(rm -f iwebbbd.o)

examples: iheatpk iheatbbd iwebpk iwebbbd

purge:
	@(rm -f iheatpk iheatbbd iwebpk iwebbbd)
