# File idax_par.mk.  Version of 8 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel IDA examples
#
# For each example requested, this makefile compiles the source and
# generates the executable
#
# Usage:  make -f idax_par.mk          [to display list of examples]
#         make -f idax_par.mk <ex>     [where <ex> is iheatpk, iheatbbd, iwebpk,
#                                       or iwebbbd]
#         make -f idax_par.mk examples [to make all examples]
#         make -f idax_par.mk purge    [to remove all example executable files]
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
	@(echo 'List of parallel IDA examples:')
	@(echo '  iheatpk  : 2-D heat equation, diagonal preconditioner')
	@(echo '  iheatbbd : 2-D heat equation, BBD preconditioner')
	@(echo '  iwebpk   : 2-D food web, block-diagonal preconditioner')
	@(echo '  iwebbbd  : 2-D food web, BBD preconditioner')

iheatpk: iheatpk.c
	@echo '...Compile iheatpk...'
	@$(CC) $(CFLAGS) -o iheatpk iheatpk.c -lida.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f iheatpk.o

iheatbbd: iheatbbd.c
	@echo '...Compile iheatbbd...'
	@$(CC) $(CFLAGS) -o iheatbbd iheatbbd.c -lida.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f iheatbbd.o

iwebpk: iwebpk.c
	@echo '...Compile iwebpk...'
	@$(CC) $(CFLAGS) -o iwebpk iwebpk.c -lida.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f iwebpk.o

iwebbbd: iwebbbd.c
	@echo '...Compile iwebbbd...'
	@$(CC) $(CFLAGS) -o iwebbbd iwebbbd.c -lida.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f iwebbbd.o

examples: iheatpk iheatbbd iwebpk iwebbbd

purge:
	@(rm -f iheatpk)
	@(rm -f iheatbbd)
	@(rm -f iwebpk)
	@(rm -f iwebbbd)

#---End of idax_par.mk---