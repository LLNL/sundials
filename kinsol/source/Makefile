# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for the serial and/or parallel versions of KINSOL package.
#
# This makefile invokes kinsols.mk or kinsolp.mk and generates a library
# file libkinsols.a or libkinsolp.a.  It first copies the appropriate
# version of the NVECTOR module files into place, then removes these files.
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
	@(echo 'Usage: make kinsols - generate the serial KINSOL')
	@(echo '                      (libkinsols.a library)')
	@(echo '       make kinsolp - generate the parallel KINSOL')
	@(echo '                      (libkinsolp.a library)')
	@(echo)

kinsols:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_serial.c $(SRC_SHARED)/nvector.c)
	@(make -f kinsols.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

kinsolp:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_parallel.c $(SRC_SHARED)/nvector.c)
	@(make -f kinsolp.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

