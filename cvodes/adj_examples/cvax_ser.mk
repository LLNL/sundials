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
	@(echo 'List of serial adjoint CVODES examples (using the serial NVECTOR module):')
	@(echo 'chemadj  - Chemical kinetics; adjoint sensitivity')
	@(echo 'adadj    - Advection-diffusion; adjoint sensitivity')
	@(echo 'foodadj  - Food web; adjoint sensitivity for G') 
	@(echo 'foodadj1 - Food web; adjoint sensitivity for g') 

chemadj: chemadj.c
	@echo '...Compile chemadj...'
	@$(CC) $(CFLAGS) -o chemadj chemadj.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f chemadj.o

adadj: adadj.c
	@echo '...Compile adadj...'
	@$(CC) $(CFLAGS) -o adadj adadj.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f adadj.o

foodadj: foodadj.c
	@echo '...Compile foodadj...'
	@$(CC) $(CFLAGS) -o foodadj foodadj.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f foodadj.o

foodadj1: foodadj1.c
	@echo '...Compile foodadj1...'
	@$(CC) $(CFLAGS) -o foodadj1 foodadj1.c -lcvodes.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f foodadj1.o

examples: chemadj adadj foodadj foodadj1

purge:
	@(rm -f chemadj)
	@(rm -f adadj)
	@(rm -f foodadj)
	@(rm -f foodadj1)

#---End of cvax_ser.mk---