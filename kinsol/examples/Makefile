# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for the KINSOL examples (serial and parallel).
#
# This makefile invokes kinxs.mk or kinxp.mk to make the executable files
# for the examples requested.  It first copies the appropriate version of
# the NVECTOR module header file into place, and then removes this file.
#
# Type 'make' to see usage instructions and list of examples.
#--------------------------------------------------------------------------

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared
INC_SHARED = $(SHARED_DIR)/include


#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo 'Usage: make kinxs   - make all serial KINSOL examples')
	@(echo '       make kinxp   - make all parallel KINSOL examples')
	@(echo '       make <ex>     - make example <ex>')
	@(echo '       make examples - make all KINSOL examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(make -f kinxs.mk)
	@(make -f kinxp.mk)

examples: kinxs kinxp

kinxs:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(make -f kinxs.mk examples)
	@(rm $(INC_SHARED)/nvector.h)

kinxp:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(make -f kinxp.mk examples)
	@(rm $(INC_SHARED)/nvector.h)

kinwebs:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(make -f kinxs.mk $@)
	@(rm $(INC_SHARED)/nvector.h)

kinwebp kinwebbbd:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(make -f kinxp.mk $@)
	@(rm $(INC_SHARED)/nvector.h)

purge:
	@(make -f kinxs.mk purge)
	@(make -f kinxp.mk purge)

