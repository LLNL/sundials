# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for the IDA examples (serial and parallel).
#
# This makefile invokes idaxs.mk or idaxp.mk to make the executable files
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
	@(echo 'Usage: make idaxs   - make all serial IDA examples')
	@(echo '       make idaxp   - make all parallel IDA examples')
	@(echo '       make <ex>     - make example <ex>')
	@(echo '       make examples - make all IDA examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(make -f idaxs.mk)
	@(make -f idaxp.mk)

examples: idaxs idaxp

idaxs:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(make -f idaxs.mk examples)
	@(rm $(INC_SHARED)/nvector.h)

idaxp:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(make -f idaxp.mk examples)
	@(rm $(INC_SHARED)/nvector.h)

irobx iheatsb iheatsk iwebsb:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(make -f idaxs.mk $@)
	@(rm $(INC_SHARED)/nvector.h)

iheatpk iheatbbd iwebpk iwebbbd:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(make -f idaxp.mk $@)
	@(rm $(INC_SHARED)/nvector.h)

purge:
	@(make -f idaxs.mk purge)
	@(make -f idaxp.mk purge)

