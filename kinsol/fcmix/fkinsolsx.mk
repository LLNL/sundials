# File fkinsolsx.mk.  Version of 14 December 2001
#--------------------------------------------------------------------------
# Makefile for the serial FKINSOL example.
#
# This makefile compiles the source, links, and removes the .o file.
#
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------


#==========================================================================
# Machine-dependent variables
#==========================================================================

#------------------------
# Fortran compiler/linker
#------------------------
FC      = f77
FLINKER = $(FC)


#==========================================================================
# Other variables
#==========================================================================

#----------
# Libraries
#----------
LIB  = ../../lib/libkinsols.a
LIBS = $(LIB) -lm

#---------------
# Compiler flags
#---------------
CFLAGS =


#==========================================================================
# Make rules
#==========================================================================

all: 
	$(FC) $(CFLAGS) -c kindiagsf.f
	$(FLINKER) -o kindiagsf kindiagsf.o $(LIBS)
	@(rm -f kindiagsf.o)

