# File fpvodex.mk.   Version of 14 December 2001.
#--------------------------------------------------------------------------
# Makefile for the FPVODE examples.
#
# This makefile makes all of the FPVODE examples.  For each one, it
# compiles the source, links, and removes the .o file.
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
# SUN cluster (SunOS)
FC       = mpif77
FLINKER  = $(FC)


#==========================================================================
# Other variables
#==========================================================================

#----------
# Libraries
#----------
LIB      = ../../lib
LIB_LIST = -lpvode $(MPI_LIBLIST) -lm
LIBS     = -L$(LIB) -L$(MPI_LIB) $(LIB_LIST)


#---------------
# Compiler flags
#---------------
CFLAGS   = -I$(MPI_INC)


#==========================================================================
# Make rules
#==========================================================================

all:
	$(FC) $(CFLAGS) -c pvdiagnf.f
	$(FLINKER) -o pvdiagnf pvdiagnf.o $(LIBS)
	@(rm -f pvdiagnf.o)
	$(FC) $(CFLAGS) -c pvdiagkf.f
	$(FLINKER) -o pvdiagkf pvdiagkf.o $(LIBS)
	@(rm -f pvdiagkf.o)
	$(FC) $(CFLAGS) -c pvdiagkbf.f
	$(FLINKER) -o pvdiagkbf pvdiagkbf.o $(LIBS)
	@(rm -f pvdiagkbf.o)
