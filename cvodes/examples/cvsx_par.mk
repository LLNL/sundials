# File cvsx_par.mk.  Version of 20 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel CVODES examples
#
# For each example requested, this makefile compiles the source and
# generates the executable
#
# Usage:  make -f cvsx_par.mk          [to display list of examples]
#         make -f cvsx_par.mk <ex>     [where <ex> is pvsnx or pvskx]
#         make -f cvsx_par.mk examples [to make all examples]
#         make -f cvsx_par.mk purge    [to remove all example executable files]
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
	@(echo '    pvsnx: 1-D advection difusion PDE;')
	@(echo '           Adams with Functional iteration')
	@(echo '    pvskx: 2-D 2-species diurnal advection-diffusion PDE;')
	@(echo '           BDF with Newton GMRES')

pvsnx: pvsnx.c
	@echo '...Compile pvsnx...'
	@$(CC) $(CFLAGS) -o pvsnx pvsnx.c -lcvodes.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvsnx.o

pvskx: pvskx.c
	@echo '...Compile pvskx...'
	@$(CC) $(CFLAGS) -o pvskx pvskx.c -lcvodes.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvskx.o

examples: pvsnx pvskx

purge:
	@(rm -f pvsnx)
	@(rm -f pvskx)

#---End of cvsx_par.mk---