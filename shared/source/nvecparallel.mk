#--------------------------------------------------------------------------------------
# File        : nvecparallel.mk
# Programmers : Radu Serban @ LLNL
# Version of  : 7 March 2002
#--------------------------------------------------------------------------------------
# Makefile for the parallel NVECTOR module.
#
#   generate the parallel NVECTOR library in sundials/shared/lib
#   add symbolic link to the parallel NVECTOR library in sundials/lib 
#   add symbolic link to nvector_parallel.h in sundials/include
#
#--------------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../include

# NVECPARALLEL library name and location
ARCH             = `uname -s`.`uname -m`
SHARED_LIB_DIR   = ../lib
NVECPAR_LIB_NAME = libnvecparallel.$(ARCH).a

#==========================================================================
# Machine-dependent variables
#==========================================================================

#----------------------
# MPI include directory
#----------------------
MPI_INC_DIR	= /usr/local/mpi/mpich/include

#-----------
# C compiler
#-----------

# CASC Sun cluster.
CC          = mpicc
CLINKER     = $(CC)

# COMPASS cluster.
#CC          = mpicc
#CLINKER     = $(CC)

# BLUE cluster
#CC          = mpcc
#CLINKER     = $(CC)

#---------------
# Compiler flags
#---------------
CFLAGS = -I$(MPI_INC_DIR) -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------


lib:
	@(echo '...Compile parallel NVECTOR object files...')
	@($(CC) $(CFLAGS) -c nvector_parallel.c)
	@(echo '...Create parallel NVECTOR library file...')
	@(ar rc $(SHARED_LIB_DIR)/$(NVECPAR_LIB_NAME) nvector_parallel.o)
	@(rm -f nvector_parallel.o)
	@(echo '...Create symbolic links to parallel NVECTOR...')
	@(cd ../../lib;     rm -f $(NVECPAR_LIB_NAME);  ln -fs ../shared/lib/$(NVECPAR_LIB_NAME) .)
	@(cd ../../include; rm -f nvector_parallel.h; ln -fs ../shared/include/nvector_parallel.h .)

purge:
	@(rm -f $(SHARED_LIB_DIR)/$(NVECPAR_LIB_NAME))
	@(cd ../../lib;     rm -f $(NVECPAR_LIB_NAME);)
	@(cd ../../include; rm -f nvector_parallel.h;)

#---End of nvecparallel.mk---