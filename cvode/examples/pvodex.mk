# File pvodex.mk.  Version of 14 December 2001
#--------------------------------------------------------------------------
# Makefile for the PVODE examples
#
# For each example requested, this makefile compiles the source, links,
# and removes the .o files.
#
# Usage:  make -f pvodex.mk   [to display list of examples]
#         make -f pvodex.mk <ex>   [where <ex> is pvnx, pvkx, or pvkxb]
#         make -f pvodex.mk examples   [to make all examples]
#         make -f pvodex.mk purge   [to remove all example executable files]
# 
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------


#==========================================================================
# Machine-dependent variables
#==========================================================================

#----------------------
# MPI include directory
#----------------------
INC_MPI	=	/usr/local/mpi/mpich/include

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


#======================================================
# Other variables
#======================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR	=	../../shared

#--------------
# Include files
#--------------
INC_SHARED	=	$(SHARED_DIR)/include
INC_SOLVER	=	../include

#---------------
# Compiler flags
#---------------
CFLAGS		=	-I$(INC_MPI) -I$(INC_SOLVER) -I$(INC_SHARED)

#----------
# Libraries
#----------
LIB		=	../../lib/
LIB_LIST	=	-lpvode -lm
LIBS		=	-L$(LIB) $(LIB_LIST)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of PVODE examples:')
	@(echo '  pvnx:  1-D advection-diffusion; nonstiff')
	@(echo '  pvkx:  2-D 2-species diurnal advection-diffusion')
	@(echo '  pvkxb: pvkx problem with BBD preconditioner')

pvnx: 
	$(CC) $(CFLAGS) -c pvnx.c
	$(CLINKER) -o pvnx pvnx.o $(LIBS)
	@(rm -f pvnx.o)

pvkx: 
	$(CC) $(CFLAGS) -c pvkx.c
	$(CLINKER) -o pvkx pvkx.o $(LIBS)
	@(rm -f pvkx.o)

pvkxb: 
	$(CC) $(CFLAGS) -c pvkxb.c
	$(CLINKER) -o pvkxb pvkxb.o $(LIBS)
	@(rm -f pvkxb.o)

examples: pvnx pvkx pvkxb

purge:
	@(rm -f pvnx)
	@(rm -f pvkx)
	@(rm -f pvkxb)
