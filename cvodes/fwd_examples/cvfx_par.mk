# File cvfx_par.mk.  Version of 20 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel CVODES examples
#
# For each example requested, this makefile compiles the source and
# generates the executable
#
# Usage:  make -f cvfx_par.mk          [to display list of examples]
#         make -f cvfx_par.mk <ex>     [where <ex> is pvfnx or pvfkx]
#         make -f cvfx_par.mk examples [to make all examples]
#         make -f cvfx_par.mk purge    [to remove all example executable files]
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
	@(echo '  pvfnx: 1-D advection difusion PDE;')
	@(echo '         Adams with Functional iteration')
	@(echo '  pvfkx: 2-D 2-species diurnal advection-diffusion PDE;')
	@(echo '         BDF with Newton GMRES')

pvfnx: pvfnx.c
	@echo '...Compile pvfnx...'
	@$(CC) $(CFLAGS) -o pvfnx pvfnx.c -lcvodes.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvfnx.o

pvfkx: pvfkx.c
	@echo '...Compile pvfkx...'
	@$(CC) $(CFLAGS) -o pvfkx pvfkx.c -lcvodes.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvfkx.o

examples: pvfnx pvfkx

purge:
	@(rm -f pvfnx)
	@(rm -f pvfkx)

#---End of cvfx_par.mk---