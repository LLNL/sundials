# File Makefile.  Version of 15 January 2002
#--------------------------------------------------------------------------
# Makefile for FPVODE, the PVODE Fortran interface package, and the
# associated examples.
#
# This makefile invokes fpvode.mk or fpvodex.mk.  In the first case, it 
# compiles the interface package and adds to the library file libpvode.a;
# it first copies the parallel version of the NVECTOR module files into
# place, and then removes these files.  In the second case, it makes
# all of the FPVODE examples.
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
	@(echo 'Usage: make fpvode  - generate FPVODE (libpvode.a library)')
	@(echo '       make fpvodex - make all FPVODE examples')
	@(echo '       make listx   - list FPVODE examples')
	@(echo '       make purge   - remove example executable files')
	@(echo)



fpvode:
	@(cp $(INC_SHARED)/nvector_parallel.h $(INC_SHARED)/nvector.h)
	@(cp $(SRC_SHARED)/nvector_parallel.c $(SRC_SHARED)/nvector.c)
	@(make -f fpvode.mk)
	@(rm $(INC_SHARED)/nvector.h $(SRC_SHARED)/nvector.c)

listx:
	@(echo 'List of FPVODE examples:')
	@(echo '  pvdiagnf:  nonstiff diagonal system')
	@(echo '  pvdiagkf:  stiff diagonal system')
	@(echo '  pvdiagkbf: stiff diagonal system with BBD preconditioner')

fpvodex:
	@(make -f fpvodex.mk)

purge:
	@(rm -f pvdiagnf)
	@(rm -f pvdiagkf)
	@(rm -f pvdiagkbf)
