# File cvx_par.mk.  Version of 5 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel CVODE examples
#
# For each example requested, this makefile compiles the source and
# generates the executable
#
# Usage:  make -f cvx_par.mk          [to display list of examples]
#         make -f cvx_par.mk <ex>     [where <ex> is pvnx, pvkx, or pvkxb]
#         make -f cvx_par.mk examples [to make all examples]
#         make -f cvx_par.mk purge    [to remove all example executable files]
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
	@(echo 'List of parallel CVODE examples (using the parallel NVECTOR module):')
	@(echo '  pvnx  :  1-D advection-diffusion; nonstiff')
	@(echo '  pvkx  :  2-D 2-species diurnal advection-diffusion')
	@(echo '  pvkxb :  pvkx problem with BBD preconditioner')

pvnx: 
	@echo '...Compile pvnx...'
	@$(CC) $(CFLAGS) -o pvnx pvnx.c -lcvode.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvnx.o

pvkx: 
	@echo '...Compile pvbx...'
	@$(CC) $(CFLAGS) -o pvkx pvkx.c -lcvode.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvkx.o

pvkxb: 
	@echo '...Compile pvkxb...'
	@$(CC) $(CFLAGS) -o pvkxb pvkxb.c -lcvode.$(ARCH) -lnvecparallel.$(ARCH) -lshared.$(ARCH) -lm
	@rm -f pvkxb.o

examples: pvnx pvkx pvkxb

purge:
	@(rm -f pvnx)
	@(rm -f pvkx)
	@(rm -f pvkxb)

#---End of cvx_par.mk---