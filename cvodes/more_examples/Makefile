# File Makefile  Version of 15 January 2002
#----------------------------------------------------
# Makefile for the CVODES/PVODES more_examples
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
	@(make -f cvodesx.mk)
	@(make -f pvodesx.mk)

advd_s chem_s diur_s food_s heat_s vpol_s:
	@(cp $(SHARED_DIR)/include/nvector_serial.h $(SHARED_DIR)/include/nvector.h)
	@(make -f cvodesx.mk $@)
	@(rm $(SHARED_DIR)/include/nvector.h)

advd_p diur_p:
	@(cp $(SHARED_DIR)/include/nvector_parallel.h $(SHARED_DIR)/include/nvector.h)
	@(make -f pvodesx.mk $@)
	@(rm $(SHARED_DIR)/include/nvector.h)

#======================================================

purge:
	@(make -f cvodesx.mk purge)
	@(make -f pvodesx.mk purge)

#======================================================
