# File fcvx_ser.mk.  Version of 27 March 2002
#--------------------------------------------------------------------------
# Makefile for the serial FCVODE example.
#
# This makefile compiles the source and generates the executable.
#
# Usage:  make -f fcvx_ser.mk          [to display list of examples]
#         make -f fcvx_ser.mk <ex>     [to make example <ex>]
#         make -f fcvx_ser.mk examples [to make all examples]
#         make -f fcvx_ser.mk purge    [to remove all example executable files]
#
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------

SHELL = /bin/sh

#--------------------------
# Top of SUNDIALS directory
#--------------------------
SUNDIALS_DIR = ../..

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

#------------------------
# Fortran compiler/linker
#------------------------

FC      = f77
FFLAGS  = -L$(LIB_DIR)

#==========================================================================
# Make rules
#==========================================================================

all:
	@(echo 'List of serial FCVODE examples (using the serial NVECTOR module):')
	@(echo '  NO examples at this time...')

#serial_ex:
#	@echo '...Compile serial_ex...'
#	@$(FC) $(FFLAGS) -o serial_ex serial_ex.f -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH)
#	@rm -f serial_ex.o

examples: serial_ex

purge:
	@(rm -f serial_ex)

#---end of fcvx_ser.mk---
