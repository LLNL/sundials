# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for the CVODE (serial) and PVODE (parallel) examples.
#
# This makefile invokes cvodex.mk or pvodex.mk to make the executable files
# for the examples requested.  It first copies the appropriate version of
# the NVECTOR module header file into place, and then removes this file.
#
# Type 'make' to see usage instructions and list of examples.
#--------------------------------------------------------------------------

#---------------------
# Path to shared files
#---------------------
SHARED_DIR	=	../../shared
INC_SHARED	=	$(SHARED_DIR)/include


#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo 'Usage: make cvodex   - make all CVODE examples')
	@(echo '       make pvodex   - make all PVODE examples')
	@(echo '       make <ex>     - make example <ex>')
	@(echo '       make examples - make all CVODE and PVODE examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(make -f cvodex.mk)
	@(make -f pvodex.mk)

examples: cvodex pvodex

cvodex:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(make -f cvodex.mk examples)
	@(rm $(INC_SHARED)/nvector.h)

pvodex:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(make -f pvodex.mk examples)
	@(rm $(INC_SHARED)/nvector.h)

cvdx cvbx cvkx cvkxb cvdemd cvdemk:
	@(cp $(INC_SHARED)/nvector_serial.h $(INC_SHARED)/nvector.h)
	@(make -f cvodex.mk $@)
	@(rm $(INC_SHARED)/nvector.h)

pvnx pvkx pvkxb:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(make -f pvodex.mk $@)
	@(rm $(INC_SHARED)/nvector.h)

purge:
	@(make -f cvodex.mk purge)
	@(make -f pvodex.mk purge)
