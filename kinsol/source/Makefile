#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 5 March 2002
#--------------------------------------------------------------------------------
# Makefile for the KINSOL package.
#
# 'make lib'   generate the KINSOL library (libkinsol.a) in sundials/kinsol/lib
#              add a symbolic link to the KINSOL library in sundials/lib
#              and symbolic links to the KINSOL header files in sundials/include
#
# 'make purge' remove KINSOL library and symbolic links
#--------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../../shared/include
KINSOL_INC_DIR  = ../include

# KINSOL library name and location
ARCH           = `uname -s`.`uname -m`
KINSOL_LIB_DIR  = ../lib
KINSOL_LIB_NAME = libkinsol.$(ARCH).a

# Source and object files
SRC    = kinbbdpre.c kinsol.c kinspgmr.c

OBJ = ${SRC:.c=.o}

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(KINSOL_INC_DIR) -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all:
	@(echo)
	@(echo 'Usage: make lib   - generate the KINSOL package library')
	@(echo '                    add symbolic link to KINSOL library')
	@(echo '                    add symbolic links to KINSOL header files')
	@(echo)
	@(echo '       make purge - remove KINSOL library and symbolic links')
	@(echo)

lib:
	@(echo 'KINSOL MODULE')
	@(echo '...Compile object files...')
	@($(CC) $(CFLAGS) -c $(SRC))
	@(echo '...Create library file...')
	@(ar rc $(KINSOL_LIB_DIR)/$(KINSOL_LIB_NAME) $(OBJ))
	@(rm -f $(OBJ))
	@(echo '...Create symbolic links...')
	@(cd ../../lib;     rm -f $(KINSOL_LIB_NAME); ln -fs ../kinsol/lib/$(KINSOL_LIB_NAME) .)
	@(cd ../../include; rm -f kinbbdpre.h; ln -fs ../kinsol/include/kinbbdpre.h .)
	@(cd ../../include; rm -f kinsol.h;    ln -fs ../kinsol/include/kinsol.h .)
	@(cd ../../include; rm -f kinspgmr.h;  ln -fs ../kinsol/include/kinspgmr.h .)

purge:
	@(echo 'KINSOL MODULE')
	@(echo '...Remove library file...')
	@(rm -f $(KINSOL_LIB_DIR)/$(KINSOL_LIB_NAME))
	@(echo '...Delete symbolic links...')
	@(cd ../../lib;     rm -f $(KINSOL_LIB_NAME))
	@(cd ../../include; rm -f kinbbdpre.h)
	@(cd ../../include; rm -f kinsol.h)
	@(cd ../../include; rm -f kinspgmr.h)

#---End of Makefile---