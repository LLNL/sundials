# File cvfx_ser.mk.  Version of 20 March 2002
#--------------------------------------------------------------------------
# Makefile for the CVODES examples using the serial N_Vector implementation.
#
# For each example requested, this makefile compiles the source and
# generates the executable.
#
# Usage:  make -f cvfx_ser.mk          [to display list of examples]
#         make -f cvfx_ser.mk <ex>     [to make example <ex>, where <ex> is 
#                                       cvfdx, cvfkx, cvfnx]
#         make -f cvfx_ser.mk examples [to make all examples]
#         make -f cvfx_ser.mk purge    [to remove all example executable files]
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
	@(echo '    cvfnx: 1-D advection difusion PDE;')
	@(echo '           Adams with Functional iteration')
	@(echo '    cvfdx: chemical kinetics ODEs;')
	@(echo '           BDF with Newton Dense')
	@(echo '    cvfkx: 2-D 2-species diurnal advection-diffusion PDE;')
	@(echo '           BDF with Newton GMRES')

cvfnx: cvfnx.c
	@echo '...Compile cvfnx...'
	@$(CC) $(CFLAGS) -o cvfnx cvfnx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvfnx.o

cvfdx: cvfdx.c
	@echo '...Compile cvfdx...'
	@$(CC) $(CFLAGS) -o cvfdx cvfdx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvfdx.o

cvfkx: cvfkx.c
	@echo '...Compile cvfkx...'
	@$(CC) $(CFLAGS) -o cvfkx cvfkx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvfkx.o

examples: cvfnx cvfdx cvfkx

purge:
	@(rm -f cvfnx)
	@(rm -f cvfdx)
	@(rm -f cvfkx)

#---End of cvfx_ser.mk---