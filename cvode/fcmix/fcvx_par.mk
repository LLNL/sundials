# File fcvx_par.mk.mk.  Version of 27 March 2002
#--------------------------------------------------------------------------
# Makefile for the parallel FCVODE example.
#
# This makefile compiles the source and generates the executable.
#
# Usage:  make -f fcvx_par.mk          [to display list of examples]
#         make -f fcvx_par.mk <ex>     [to make example <ex>, where <ex>
#                                       pvdiagnf, pvdiagkf, or pvdiagkbf ]
#         make -f fcvx_par.mk examples [to make all examples]
#         make -f fcvx_par.mk purge    [to remove all example executable files]
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
	@(echo 'List of parallel FCVODE examples (using the parallel NVECTOR module):')
	@(echo '  pvdiagnf  : Diagonal ODE example. Non-stiff case (ADAMS/FUNCTIONAL)')
	@(echo '  pvdiagkf  : Diagonal ODE example. Stiff case (BDF/SPGMR)')
	@(echo '  pvdiagkbf : Diagonal ODE example. Stiff case (BDF/SPGMR with FCVBBD preconditioner)')

pvdiagnf: 
	@echo '...Compile pvdiagnf...'
	@$(FC) $(FFLAGS) -o pvdiagnf pvdiagnf.f -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecparallel.$(ARCH)
	@rm -f pvdiagnf.o

pvdiagkf: 
	@echo '...Compile pvdiagkf...'
	@$(FC) $(FFLAGS) -o pvdiagkf pvdiagkf.f -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecparallel.$(ARCH)
	@rm -f pvdiagkf.o

pvdiagkbf: 
	@echo '...Compile pvdiagkbf...'
	@$(FC) $(FFLAGS) -o pvdiagkbf pvdiagkbf.f -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecparallel.$(ARCH)
	@rm -f pvdiagkbf.o

examples: pvdiagnf pvdiagkf pvdiagkbf

purge:
	@(rm -f pvdiagnf)
	@(rm -f pvdiagkf)
	@(rm -f pvdiagkbf)

#---end of fcvx_par.mk---