# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for FKINSOL, the Fortran interface package for KINSOL and the
# associated examples.
#
# This makefile invokes fkinsols.mk, fkinsolp.mk, fkinsolsx.mk, or
# fkinsolpx.mk.  In the first two cases, it adds to the library
# file libkinsols.a or libkinsolp.a; it first copies the appropriate
# version of the NVECTOR module files into place, then removes these files.
# In the second two cases, it generates the example requested.
#
# Type 'make' to see usage instructions.
#--------------------------------------------------------------------------


#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared
INC_SHARED = $(SHARED_DIR)/include
SRC_SHARED = $(SHARED_DIR)/source


#==========================================================================
# Make rules
#==========================================================================

all:
	@(echo)
	@(echo 'Usage: make fkinsols - generate the serial FKINSOL')
	@(echo '                       (libkinsols.a library)')
	@(echo '       make fkinsolp - generate the parallel FKINSOL')
	@(echo '                       (libkinsolp.a library)')
	@(echo '       make fkinsolsx - generate the serial FKINSOL example')
	@(echo '       make fkinsolpx - generate the parallel FKINSOL example')
	@(echo '       make purge     - remove example executable files')
	@(echo)

fkinsols:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_serial.c $(SRC_SHARED)/nvector.c)
	@(make -f fkinsols.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

fkinsolp:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_parallel.c $(SRC_SHARED)/nvector.c)
	@(make -f fkinsolp.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

fkinsolsx:
	@(make -f fkinsolsx.mk)

fkinsolpx:
	@(make -f fkinsolpx.mk)

purge:
	@(rm -f kindiagsf)
	@(rm -f kindiagpf)
