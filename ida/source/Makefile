# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for the serial and/or parallel versions of IDA package.
#
# This makefile invokes idas.mk or idap.mk and generates a library
# file libidas.a or libidap.a.  It first copies the appropriate
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
	@(echo 'Usage: make idas - generate the serial IDA')
	@(echo '                      (libidas.a library)')
	@(echo '       make idap - generate the parallel IDA')
	@(echo '                      (libidap.a library)')
	@(echo)

idas:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_serial.c $(SRC_SHARED)/nvector.c)
	@(make -f idas.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

idap:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_parallel.c $(SRC_SHARED)/nvector.c)
	@(make -f idap.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

