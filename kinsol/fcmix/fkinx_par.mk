# File fkinx_par.mk.mk.  Version of 12 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel FKINSOL example.
#
# This makefile compiles the source and generates the executable.
#
# Usage:  make -f fkinx_par.mk          [to display list of examples]
#         make -f fkinx_par.mk <ex>     [to make example <ex>, where <ex> is kindiagpf]
#         make -f fkinx_par.mk examples [to make all examples]
#         make -f fkinx_par.mk purge    [to remove all example executable files]
#
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------

SHELL = /bin/sh

#--------------------------
# Top of SUNDIALS directory
#--------------------------
SUNDIALS_DIR = ../..

#----------------------
# Path to library files
#----------------------
LIB_DIR = $(SUNDIALS_DIR)/lib

#-------------
# Architecture
#-------------
ARCH = `uname -s`.`uname -m`

#==========================================================================
# Machine-dependent variables
#==========================================================================

#----------------------------------
# MPI library and include directory
#----------------------------------
MPI_LIB     = /usr/local/mpi/mpich/lib
MPI_INC     = /usr/local/mpi/mpich/include
MPI_LIBLIST = 

#------------------------
# Fortran compiler/linker
#------------------------
# Sun cluster
FC      = mpif77
FFLAGS  = -L$(LIB_DIR) -L$(MPI_LIB)


#==========================================================================
# Make rules
#==========================================================================

all:
	@(echo 'List of parallel FKINSOL examples (using the parallel NVECTOR module):')
	@(echo '  kindiagpf:  simple diagonal test with Fortran interface')

kindiagpf: 
	@echo '...Compile kindiagps...'
	@$(FC) $(FFLAGS) -o kindiagpf kindiagpf.f -lkinsol.$(ARCH) -lshared.$(ARCH) -lnvecparallel.$(ARCH)
	@rm -f kindiagpf.o

examples: kindiagpf

purge:
	@(rm -f kindiagpf)

#---end of fkinx_par.mk---