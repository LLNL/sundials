# File Makefile  Version of 15 January 2002
#----------------------------------------------------
# Makefile for the CVODES/PVODES examples
# (serial and parallel versions)
#----------------------------------------------------


#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo 'Usage: make cvodesx   - make all CVODE examples')
	@(echo '       make pvodesx   - make all PVODE examples')
	@(echo '       make <ex>      - make example <ex>')
	@(echo '       make examples  - make all CVODE and PVODE examples')
	@(echo '       make purge     - remove all example executable files')
	@(echo)
	@(make -f cvodesx.mk)
	@(make -f pvodesx.mk)

examples: cvodesx pvodesx

cvodesx:
	@(cp $(SHARED_DIR)/include/nvector_serial.h $(SHARED_DIR)/include/nvector.h)
	@(make -f cvodesx.mk examples)
	@(rm $(SHARED_DIR)/include/nvector.h)

pvodesx:
	@(cp $(SHARED_DIR)/include/nvector_parallel.h $(SHARED_DIR)/include/nvector.h)
	@(make -f pvodesx.mk examples)
	@(rm $(SHARED_DIR)/include/nvector.h)

cvsnx cvsdx cvskx:
	@(cp $(SHARED_DIR)/include/nvector_serial.h $(SHARED_DIR)/include/nvector.h)
	@(make -f cvodesx.mk $@)
	@(rm $(SHARED_DIR)/include/nvector.h)

pvsnx pvskx:
	@(cp $(SHARED_DIR)/include/nvector_parallel.h $(SHARED_DIR)/include/nvector.h)
	@(make -f pvodesx.mk $@)
	@(rm $(SHARED_DIR)/include/nvector.h)

#======================================================

purge:
	@(make -f cvodesx.mk purge)
	@(make -f pvodesx.mk purge)

#======================================================
