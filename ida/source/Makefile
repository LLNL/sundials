#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 7 March 2002
#--------------------------------------------------------------------------------
# Makefile for the IDA package.
#
# 'make lib'   generate the IDA library (libida.a) in sundials/ida/lib
#              add a symbolic link to the IDA library in sundials/lib
#              and symbolic links to the IDA header files in sundials/include
#
# 'make purge' remove CVODE library and symbolic links
#--------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../../shared/include
IDA_INC_DIR    = ../include

# CVODE library name and location
ARCH         = `uname -s`.`uname -m`
IDA_LIB_DIR  = ../lib
IDA_LIB_NAME = libida.$(ARCH).a

# Source and object files
SRC    = idaband.c idabbdpre.c ida.c idadense.c idaspgmr.c
OBJ    = ${SRC:.c=.o}

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(IDA_INC_DIR) -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all:
	@(echo)
	@(echo 'Usage: make lib - generate the IDA package library')
	@(echo '                  add symbolic link to IDA library')
	@(echo '                  add symbolic links to IDA header files')
	@(echo)
	@(echo '       make purge - remove IDA library and symbolic links')
	@(echo)

lib:
	@(echo 'IDA MODULE')
	@(echo '...Compile object files...')
	@($(CC) $(CFLAGS) -c $(SRC))
	@(echo '...Create library file...')
	@(ar rc $(IDA_LIB_DIR)/$(IDA_LIB_NAME) $(OBJ))
	@(rm -f $(OBJ))
	@(echo '...Create symbolic links...')
	@(cd ../../lib;     rm -f $(IDA_LIB_NAME); ln -fs ../ida/lib/$(IDA_LIB_NAME) .)
	@(cd ../../include; rm -f idaband.h;   ln -fs ../ida/include/idaband.h .)
	@(cd ../../include; rm -f idabbdpre.h; ln -fs ../ida/include/idabbdpre.h .)
	@(cd ../../include; rm -f idadense.h;  ln -fs ../ida/include/idadense.h .)
	@(cd ../../include; rm -f ida.h;       ln -fs ../ida/include/ida.h .)
	@(cd ../../include; rm -f idaspgmr.h;  ln -fs ../ida/include/idaspgmr.h .)

purge:
	@(echo 'IDA MODULE')
	@(echo '...Remove library file...')
	@(rm -f $(IDA_LIB_DIR)/$(IDA_LIB_NAME))
	@(echo '...Delete symbolic links...')
	@(cd ../../lib;     rm -f $(IDA_LIB_NAME);)
	@(cd ../../include; rm -f idaband.h;)
	@(cd ../../include; rm -f idabbdpre.h;)
	@(cd ../../include; rm -f idadense.h;)
	@(cd ../../include; rm -f ida.h;)
	@(cd ../../include; rm -f idaspgmr.h;)

#---End of Makefile---