# File cvsxmake  Version of 12 November 2001
#---------------------------------------------------
# Makefile for the CVODES examples (serial versions)
#---------------------------------------------------

#======================================================
# Machine dependent variables
#======================================================

#-----------
# C compiler
#-----------
CC      = gcc
CLINKER = $(CC)

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
INC_CVS    = ../include

#---------------
# Compiler flags
#---------------
CFLAGS     = -I$(INC_CVS) -I$(INC_SHARED)

#----------
# Libraries
#----------
CVSLIB     = ../../lib/
LIB_LIST   = -lcvodes -lm
LIBS       = -L$(CVSLIB) $(LIB_LIST)

#======================================================
# Make rules
#======================================================

all: 
	@(echo)
	@(echo 'List of examples (serial)')
	@(echo '-------------------------')
	@(echo '    advd_s: 1-D advection-diffusion')
	@(echo '    chem_s: Chemical kinetics')
	@(echo '    diur_s: 2-D 2-species diurnal advection-diffusion')
	@(echo '    heat_s: 2-D heat PDE')
	@(echo '    vpol_s: Van der Pol oscillator')
	@(echo)

advd_s: advd_s.o 
	$(CLINKER) -o advd_s advd_s.o $(LIBS)
	@(rm -f advd_s.o)

chem_s: chem_s.o
	$(CLINKER) -o chem_s chem_s.o $(LIBS)
	@(rm -f chem_s.o)

diur_s: diur_s.o
	$(CLINKER) -o diur_s diur_s.o $(LIBS)
	@(rm -f diur_s.o)

food_s: food_s.o
	$(CLINKER) -o food_s food_s.o $(LIBS)
	@(rm -f food_s.o)

heat_s: heat_s.o
	$(CLINKER) -o heat_s heat_s.o $(LIBS)
	@(rm -f heat_s.o)

vpol_s: vpol_s.o
	$(CLINKER) -o vpol_s vpol_s.o $(LIBS)
	@(rm -f vpol_s.o)

#======================================================

purge:
	@(rm -f advd_s)
	@(rm -f chem_s)
	@(rm -f diur_s)
	@(rm -f food_s)
	@(rm -f heat_s)
	@(rm -f vpol_s)

#======================================================

%.o : %.c
	$(CC) $(CFLAGS) -c $<

#======================================================
