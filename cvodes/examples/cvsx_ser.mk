# File cvsx_ser.mk.  Version of 20 March 2002
#--------------------------------------------------------------------------
# Makefile for the CVODES examples using the serial N_Vector implementation.
#
# For each example requested, this makefile compiles the source and
# generates the executable.
#
# Usage:  make -f cvsx_ser.mk          [to display list of examples]
#         make -f cvsx_ser.mk <ex>     [to make example <ex>, where <ex> is 
#                                       cvsdx, cvskx, cvsnx]
#         make -f cvsx_ser.mk examples [to make all examples]
#         make -f cvsx_ser.mk purge    [to remove all example executable files]
# 
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------

SHELL = /bin/sh

#--------------------------
# Top of SUNDIALS directory
#--------------------------
SUNDIALS_DIR = ../..

#---------------------
# Path to header files
#---------------------
INC_DIR	= $(SUNDIALS_DIR)/include

#----------------------
# Path to library files
#----------------------
LIB_DIR = $(SUNDIALS_DIR)/lib

#-------------
# Architecture
#-------------
ARCH = `uname -s`.`uname -m`

#======================================================
# Machine-dependent variables
#======================================================

#------------------------------
# C compiler and compiler flags
#------------------------------

CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(INC_DIR) -L$(LIB_DIR)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of serial CVODES examples (using the serial NVECTOR module):')
	@(echo '    cvsnx: 1-D advection difusion PDE;')
	@(echo '           Adams with Functional iteration')
	@(echo '    cvsdx: chemical kinetics ODEs;')
	@(echo '           BDF with Newton Dense')
	@(echo '    cvskx: 2-D 2-species diurnal advection-diffusion PDE;')
	@(echo '           BDF with Newton GMRES')

cvsnx: cvsnx.c
	@echo '...Compile cvsnx...'
	@$(CC) $(CFLAGS) -o cvsnx cvsnx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvsnx.o

cvsdx: cvsdx.c
	@echo '...Compile cvsdx...'
	@$(CC) $(CFLAGS) -o cvsdx cvsdx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvsdx.o

cvskx: cvskx.c
	@echo '...Compile cvsnx...'
	@$(CC) $(CFLAGS) -o cvskx cvskx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvskx.o

examples: cvsnx cvsdx cvskx

purge:
	@(rm -f cvsnx)
	@(rm -f cvsdx)
	@(rm -f cvskx)

#---End of cvsx_ser.mk---