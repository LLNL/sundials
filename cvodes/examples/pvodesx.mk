# File pvxmake  Version of 12 November 2001
#----------------------------------------------------
# Makefile for the PVODES examples (parallel versions)
#----------------------------------------------------

#======================================================
# Machine dependent variables
#======================================================

#-----------
# C compiler
#-----------

# CASC Sun cluster.
CC      = mpicc
CLINKER = $(CC)

# COMPASS cluster.
#CC          = mpicc
#CLINKER     = $(CC)

# BLUE cluster
#CC          = mpcc
#CLINKER     = $(CC)

#======================================================
# Machine independent variables
#======================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

#--------------
# Include files
#--------------
INC_SHARED = $(SHARED_DIR)/include
INC_PVS    = ../include

#---------------
# Compiler flags
#---------------
CFLAGS     = -I$(INC_PVS) -I$(INC_SHARED)

#----------
# Libraries
#----------
PVSLIB     = ../../lib/
LIB_LIST   = -lpvodes -lm
LIBS       = -L$(PVSLIB) $(LIB_LIST)

#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo 'List of examples (parallel)')
	@(echo '---------------------------')
	@(echo '    pvsnx: 1-D advection difusion PDE;')
	@(echo '           Adams with Functional iteration')
	@(echo '    pvskx: 2-D 2-species diurnal advection-diffusion PDE;')
	@(echo '           BDF with Newton GMRES')
	@(echo)

examples: pvsnx pvskx

pvsnx: pvsnx.o
	$(CLINKER) -o pvsnx pvsnx.o $(LIBS)
	@(rm -f pvsnx.o)

pvskx: pvskx.o
	$(CLINKER) -o pvskx pvskx.o $(LIBS)
	@(rm -f pvskx.o)

#======================================================

purge:
	@(rm -f pvsnx)
	@(rm -f pvskx)

#======================================================

%.o : %.c
	$(CC) $(CFLAGS) -c $<

#======================================================
