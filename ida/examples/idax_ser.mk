# File idax_ser.mk.  Version of 8 March 2002
#--------------------------------------------------------------------------
# Makefile for the serial IDA example using the serial N_Vector implementation.
#
# For each example requested, this makefile compiles the source and
# generates the executable.
#
# Usage:  make -f idax_ser.mk          [to display list of examples]
#         make -f idax_ser.mk <ex>     [where <ex> is irobx, iheatsb, iheatsk,
#                                       or iwebsb]
#         make -f idax_ser.mk examples [to make all examples]
#         make -f idax_ser.mk purge    [to remove all example executable files]
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
	@(echo 'List of serial IDA examples (using the serial NVECTOR module):')
	@(echo '  irobx    :  3-species Robertson kinetics system')
	@(echo '  iheatsb  :  2-D heat equation, banded Jacobian')
	@(echo '  iheatsk  :  2-D heat equation, diagonal preconditioner')
	@(echo '  iwebsb   :  2-D food web system, banded Jacobian')

irobx: irobx.c
	@echo '...Compile irobx...'
	@$(CC) $(CFLAGS) -o irobx irobx.c -lida.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f irobx.o

iheatsb: iheatsb.c
	@echo '...Compile iheatsb...'
	@$(CC) $(CFLAGS) -o iheatsb iheatsb.c -lida.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f iheatsb.o

iheatsk: iheatsk.c
	@echo '...Compile iheatsk...'
	@$(CC) $(CFLAGS) -o iheatsk iheatsk.c -lida.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f iheatsk.o

iwebsb: iwebsb.c
	@echo '...Compile iwebsb...'
	@$(CC) $(CFLAGS) -o iwebsb iwebsb.c -lida.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f iwebsb.o

examples: irobx iheatsb iheatsk iwebsb

purge:
	@(rm -f irobx)
	@(rm -f iheatsb)
	@(rm -f iheatsk)
	@(rm -f iwebsb)

#---End of idax_ser.mk---