# File Makefile  Version of 7 February 2002

# Makefile for the serial and/or parallel versions 
# of the CVODES/PVODES package

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

all:
	@(echo)
	@(echo 'Usage: make cvodes - generates the serial version of the package')
	@(echo '                    (libcvodes.a library)')
	@(echo '       make cvodea - generates the serial adjoint version of the package')
	@(echo '                    (libcvodea.a library)')
	@(echo '       make pvodes - generates the parallel version of the package')
	@(echo '                    (libpvodes.a library)')
	@(echo)

cvodes:
	@(cp $(SHARED_DIR)/source/nvector_serial.c $(SHARED_DIR)/source/nvector.c)
	@(cp $(SHARED_DIR)/include/nvector_serial.h $(SHARED_DIR)/include/nvector.h)
	@(make -f cvodes.mk)
	@(rm $(SHARED_DIR)/source/nvector.c)
	@(rm $(SHARED_DIR)/include/nvector.h)

cvodea:
	@(cp $(SHARED_DIR)/source/nvector_serial.c $(SHARED_DIR)/source/nvector.c)
	@(cp $(SHARED_DIR)/include/nvector_serial.h $(SHARED_DIR)/include/nvector.h)
	@(make -f cvodea.mk)
	@(rm $(SHARED_DIR)/source/nvector.c)
	@(rm $(SHARED_DIR)/include/nvector.h)

pvodes:
	@(cp $(SHARED_DIR)/source/nvector_parallel.c $(SHARED_DIR)/source/nvector.c)
	@(cp $(SHARED_DIR)/include/nvector_parallel.h $(SHARED_DIR)/include/nvector.h)
	@(make -f pvodes.mk)
	@(rm $(SHARED_DIR)/source/nvector.c)
	@(rm $(SHARED_DIR)/include/nvector.h)

