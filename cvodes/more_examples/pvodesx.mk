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
	@(echo '    advd_p: 1-D advection-diffusion')
	@(echo '    diur_p: 2-D 2-species diurnal advection-diffusion')
	@(echo)

advd_p: advd_p.o
	$(CLINKER) -o advd_p advd_p.o $(LIBS)
	@(rm -f advd_p.o)

diur_p: diur_p.o
	$(CLINKER) -o diur_p diur_p.o $(LIBS)
	@(rm -f diur_p.o)

#======================================================

purge:
	@(rm -f advd_p)
	@(rm -f diur_p)

#======================================================

%.o : %.c
	$(CC) $(CFLAGS) -c $<

#======================================================
