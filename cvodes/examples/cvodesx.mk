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
CFLAGS     = -g -I$(INC_CVS) -I$(INC_SHARED)

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
	@(echo '    cvsnx: 1-D advection difusion PDE;')
	@(echo '           Adams with Functional iteration')
	@(echo '    cvsdx: chemical kinetics ODEs;')
	@(echo '           BDF with Newton Dense')
	@(echo '    cvskx: 2-D 2-species diurnal advection-diffusion PDE;')
	@(echo '           BDF with Newton GMRES')
	@(echo)

examples: cvsnx cvsdx cvskx

cvsnx: cvsnx.o
	$(CLINKER) -o cvsnx cvsnx.o $(LIBS)
	@(rm -f cvsnx.o)

cvsdx: cvsdx.o
	$(CLINKER) -o cvsdx cvsdx.o $(LIBS)
	@(rm -f cvsdx.o)

cvskx: cvskx.o
	$(CLINKER) -o cvskx cvskx.o $(LIBS)
	@(rm -f cvskx.o)

#======================================================

purge:
	@(rm -f cvsnx)
	@(rm -f cvsdx)
	@(rm -f cvskx)

#======================================================

%.o : %.c
	$(CC) $(CFLAGS) -c $<

#======================================================
