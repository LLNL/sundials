# File Makefile  Version of 4 Feruary 2002
#----------------------------------------------------
# Makefile for the CVODEA examples
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
	@(make -f cvodeax.mk)

chemfwd chemadj ad adadj diurnadj:
	@(cp $(SHARED_DIR)/include/nvector_serial.h $(SHARED_DIR)/include/nvector.h)
	@(make -f cvodeax.mk $@)
	@(rm $(SHARED_DIR)/include/nvector.h)

#======================================================

purge:
	@(make -f cvodeax.mk purge)

#======================================================
