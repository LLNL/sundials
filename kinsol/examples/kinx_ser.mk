# File kinx_ser.mk.  Version of 7 March 2002
#--------------------------------------------------------------------------
# Makefile for the serial KINSOL examples.
#
# For the example requested, this makefile compiles the source and
# generates the executable.
#
# Usage:  make -f kinx_ser.mk          [to display list of examples]
#         make -f kinx_ser.mk <ex>     [to make example <ex>, where <ex> is kinwebs]
#         make -f kinx_ser.mk examples [to make all examples]
#         make -f kinx_ser.mk purge    [to remove all example executable files]
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

#==========================================================================
# Machine-dependent variables
#==========================================================================

#-----------
# C compiler
#-----------
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(INC_DIR) -L$(LIB_DIR)



#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of serial KINSOL examples (using the serial NVECTOR module):')
	@(echo '  kinwebs:  2-D food web system, block-diagonal preconditioner')

kinwebs: kinwebs.c
	@echo '...Compile kinwebs...'
	@$(CC) $(CFLAGS) -o kinwebs kinwebs.c -lkinsol.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f kinwebs.o

examples: kinwebs

purge:
	@(rm -f kinwebs)

#---End of kinx_ser.mk---