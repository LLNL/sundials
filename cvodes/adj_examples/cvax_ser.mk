# File cvax_ser.mk  Version of 22 March 2002
#------------------------------------------------------------------------------------
# Makefile for the adjoint CVODES examples  using the serial N_Vector implementation.
#
# For each example requested, this makefile compiles the source and
# generates the executable.
#
# Usage:  make -f cvax_ser.mk          [to display list of examples]
#         make -f cvax_ser.mk <ex>     [to make example <ex>]
#         make -f cvax_ser.mk examples [to make all examples]
#         make -f cvax_ser.mk purge    [to remove all example executable files]
# 
# Check, and modify if necessary, the machine-dependent variables below.
#------------------------------------------------------------------------------------

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
	@(echo '  cvadx  : Chemical kinetics; adjoint sensitivity')
	@(echo '  cvabx  : Advection-diffusion; adjoint sensitivity')
	@(echo '  cvakx  : Food web; adjoint sensitivity for G') 
	@(echo '  cvakxb : Food web; adjoint sensitivity for g') 

cvadx: cvadx.c
	@echo '...Compile cvadx...'
	@$(CC) $(CFLAGS) -o cvadx cvadx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvadx.o

cvabx: cvabx.c
	@echo '...Compile cvabx...'
	@$(CC) $(CFLAGS) -o cvabx cvabx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvabx.o

cvakx: cvakx.c
	@echo '...Compile cvakx...'
	@$(CC) $(CFLAGS) -o cvakx cvakx.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvakx.o

cvakxb: cvakxb.c
	@echo '...Compile cvakxb...'
	@$(CC) $(CFLAGS) -o cvakxb cvakxb.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvakxb.o

examples: cvadx cvabx cvakx cvakxb

purge:
	@(rm -f cvadx)
	@(rm -f cvabx)
	@(rm -f cvakx)
	@(rm -f cvakxb)

#---End of cvax_ser.mk---