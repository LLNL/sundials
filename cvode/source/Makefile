# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for the serial and/or parallel packages CVODE/PVODE.
#
# This makefile invokes cvode.mk or pvode.mk and generates a library file
# libcvode.a or libpvode.a.  It first copies the appropriate version of
# the NVECTOR module files into place, and then removes these files.
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
	@(echo 'Usage: make cvode - generate the serial package CVODE')
	@(echo '                    (libcvode.a library)')
	@(echo '       make pvode - generate the parallel package PVODE')
	@(echo '                    (libpvode.a library)')
	@(echo)

cvode:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_serial.c $(SRC_SHARED)/nvector.c)
	@(make -f cvode.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

pvode:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_parallel.c $(SRC_SHARED)/nvector.c)
	@(make -f pvode.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

