# File cvxmake  Version of 12 November 2001
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
CFLAGS = -ffloat-store -Wall -g -I$(INC_CVS) -I$(INC_SHARED)
#CFLAGS = -I$(INC_CVS) -I$(INC_SHARED)

#----------
# Libraries
#----------
CVSLIB   = ../../lib/

#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo 'chemfwd - Chemical kinetics; forward sensitivity')
	@(echo 'chemadj - Chemical kinetics; adjoint sensitivity')
	@(echo)
	@(echo 'ad      - Advection-diffusion; simulation')
	@(echo 'adadj   - Advection-diffusion; adjoint sensitivity')
	@(echo)
	@(echo 'diuirnadj - Diurnal adv-dif; adjoint sensitivity')
	@(echo)

examples: chemadj adadj

chemfwd: chemfwd.o
	$(CLINKER) -o chemfwd chemfwd.o -L$(CVSLIB) -lcvodes -lm
	@(rm -f chemfwd.o)

chemadj: chemadj.o
	$(CLINKER) -o chemadj chemadj.o -L$(CVSLIB) -lcvodea -lm
	@(rm -f chemadj.o)

ad: ad.o
	$(CLINKER) -g -o ad ad.o -L$(CVSLIB) -lcvodes -lm
	@(rm -f ad.o)

adadj: adadj.o
	$(CLINKER) -g -o adadj adadj.o -L$(CVSLIB) -lcvodea -lm
	@(rm -f adadj.o)

diurnadj: diurnadj.o
	$(CLINKER) -g -o diurnadj diurnadj.o -L$(CVSLIB) -lcvodea -lm
	@(rm -f diurnadj.o)
#======================================================

purge:
	@(rm -f chemfwd)
	@(rm -f chemadj)
	@(rm -f ad)
	@(rm -f adadj)
	@(rm -f diurnadj)

#======================================================

%.o : %.c
	$(CC) $(CFLAGS) -c $<

#======================================================
