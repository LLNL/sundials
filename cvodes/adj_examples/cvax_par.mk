# File cvax_par.mk.  Version of 22 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel adjoint CVODES examples
#
# For each example requested, this makefile compiles the source and
# generates the executable
#
# Usage:  make -f cvax_par.mk          [to display list of examples]
#         make -f cvax_par.mk <ex>     [make example <ex>]
#         make -f cvax_par.mk examples [to make all examples]
#         make -f cvax_par.mk purge    [to remove all example executable files]
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
	@(echo 'List of parallel CVODES examples (using the parallel NVECTOR module):')
	@(echo '  pvanx  : 1-D advection-diffusion; adjoint sensitivity')

examples: pvanx

pvanx:
	@echo '...Compile pvanx...'
	@$(CC) $(CFLAGS) -o pvanx pvanx.c -lcvodes.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvanx.o

purge:
	@(rm -f pvanx)

#---End of cvax_par.mk---