# File fkinsolpx.mk.  Version of 14 December 2001
#--------------------------------------------------------------------------
# Makefile for the parallel FKINSOL example.
#
# This makefile compiles the source, links, and removes the .o file.
#
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------


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
FLINKER = $(FC)


#==========================================================================
# Other variables
#==========================================================================

#----------
# Libraries
#----------
LIB      = ../../lib
LIB_LIST = -lkinsolp $(MPI_LIBLIST) -lm
LIBS     = -L$(LIB) -L$(MPI_LIB) $(LIB_LIST)


#---------------
# Compiler flags
#---------------
CFLAGS   = -I$(MPI_INC)


#==========================================================================
# Make rules
#==========================================================================

all: 
	$(FC) -c $(CFLAGS) kindiagpf.f
	$(FLINKER) -o kindiagpf kindiagpf.o $(LIBS)
	@(rm -f kindiagpf.o)

