#--------------------------------------------------------------------------------------
# File        : nvecserial.mk
# Programmers : Radu Serban @ LLNL
# Version of  : 7 March 2002
#--------------------------------------------------------------------------------------
# Makefile for the serial NVECTOR module.
#
#   generate the serial NVECTOR library in sundials/shared/lib
#   add symbolic link to the serial NVECTOR library in sundials/lib 
#   add symbolic link to nvector_serial.h in sundials/include
#
#--------------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../include

# NVECSERIAL library name and location
ARCH             = `uname -s`.`uname -m`
SHARED_LIB_DIR   = ../lib
NVECSER_LIB_NAME = libnvecserial.$(ARCH).a

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

lib:
	@(echo '...Compile NVECTOR serial object files...')
	@($(CC) $(CFLAGS) -c nvector_serial.c)
	@(echo '...Create NVECTOR serial library file...')
	@(ar rc $(SHARED_LIB_DIR)/$(NVECSER_LIB_NAME) nvector_serial.o)
	@(rm -f nvector_serial.o)
	@(echo '...Create symbolic links to serial NVECTOR...')
	@(cd ../../lib;     rm -f $(NVECSER_LIB_NAME);  ln -fs ../shared/lib/$(NVECSER_LIB_NAME) .)
	@(cd ../../include; rm -f nvector_serial.h; ln -fs ../shared/include/nvector_serial.h .)

purge:
	@(rm -f $(SHARED_LIB_DIR)/$(NVECSER_LIB_NAME))
	@(cd ../../lib;     rm -f $(NVECSER_LIB_NAME);)
	@(cd ../../include; rm -f nvector_serial.h;)

#---End of nvecserial.mk---