# File fcvx_ser.mk.  Version of 22 July 2002
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------

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

LIBS = -lfcvode.$(ARCH) -lcvode.$(ARCH) \
       -lfnvecserial.$(ARCH) -lnvecserial.$(ARCH) \
       -lshared.$(ARCH) 

#==========================================================================
# Make rules
#==========================================================================

all:
	@(echo 'List of serial FCVODE examples (using the serial NVECTOR module):')
	@(echo '  cvdensef  : chemical kinetics example (BDF/DENSE)')
	@(echo '  cvbandf   : advection-diffusion example (BDF/BAND)')
	@(echo '  cvkryf    : kinetics-transport example (BDF/SPGMR)')

cvdensef:
	@echo '...Compile cvdensef...'
	@$(FC) $(FFLAGS) -o cvdensef cvdensef.f $(LIBS)
	@rm -f cvdensef.o

cvbandf:
	@echo '...Compile cvbandf...'
	@$(FC) $(FFLAGS) -o cvbandf cvbandf.f $(LIBS)
	@rm -f cvbandf.o

cvkryf:
	@echo '...Compile cvkryf...'
	@$(FC) $(FFLAGS) -o cvkryf cvkryf.f $(LIBS)
	@rm -f cvkryf.o

examples: cvdensef cvbandf cvkryf

purge:
	@(rm -f cvdensef cvbandf cvkryf)

#---end of fcvx_ser.mk---
