# File fkinx_ser.mk.  Version of 12 March 2002
#--------------------------------------------------------------------------
# Makefile for the serial FKINSOL example.
#
# This makefile compiles the source and generates the executable.
#
# Usage:  make -f fkinx_ser.mk          [to display list of examples]
#         make -f fkinx_ser.mk <ex>     [to make example <ex>, where <ex> is kindiagsf]
#         make -f fkinx_ser.mk examples [to make all examples]
#         make -f fkinx_ser.mk purge    [to remove all example executable files]
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
# Common libraries
#==========================================================================

LIBS = -lfkinsol.$(ARCH) -lkinsol.$(ARCH) \
       -lfnvecserial.$(ARCH) -lnvecserial.$(ARCH) \
       -lshared.$(ARCH) 

#==========================================================================
# Make rules
#==========================================================================

all:
	@(echo 'List of serial FKINSOL examples (using the serial NVECTOR module):')
	@(echo '  kindiagsf:  simple diagonal test with Fortran interface')

kindiagsf:
	@echo '...Compile kindiagsf...'
	@$(FC) $(FFLAGS) -o kindiagsf kindiagsf.f $(LIBS)
	@rm -f kindiagsf.o

examples: kindiagsf

purge:
	@(rm -f kindiagsf)

#---end of fkinx_ser.mk---
