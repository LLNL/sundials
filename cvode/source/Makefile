#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 7 March 2002
#--------------------------------------------------------------------------------
# Makefile for the CVODE package.
#
# 'make lib' generate the CVODE library (libcvode.a) in sundials/cvode/lib
#            add a symbolic link to the CVODE library in sundials/lib
#            and symbolic links to the CVODE header files in sundials/include
#
# 'make purge' remove CVODE library and symbolic links
#--------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../../shared/include
CVODE_INC_DIR  = ../include

# CVODE library name and location
ARCH           = `uname -s`.`uname -m`
CVODE_LIB_DIR  = ../lib
CVODE_LIB_NAME = libcvode.$(ARCH).a

# Source and object files
SRC    = cvband.c cvbandpre.c cvdense.c cvdiag.c cvode.c cvspgmr.c cvbbdpre.c

OBJ = ${SRC:.c=.o}

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(CVODE_INC_DIR) -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all:
	@(echo)
	@(echo 'Usage: make lib   - generate the CVODE package library')
	@(echo '                    add symbolic link to CVODE library')
	@(echo '                    add symbolic links to CVODE header files')
	@(echo)
	@(echo '       make purge - remove CVODE library and symbolic links')
	@(echo)

lib:
	@(echo 'CVODE MODULE')
	@(echo '...Compile object files...')
	@($(CC) $(CFLAGS) -c $(SRC))
	@(echo '...Create library file...')
	@(ar rc $(CVODE_LIB_DIR)/$(CVODE_LIB_NAME) $(OBJ))
	@(rm -f $(OBJ))
	@(echo '...Create symbolic links...')
	@(cd ../../lib;     rm -f $(CVODE_LIB_NAME);  ln -fs ../cvode/lib/$(CVODE_LIB_NAME) .)
	@(cd ../../include; rm -f cvband.h;    ln -fs ../cvode/include/cvband.h .)
	@(cd ../../include; rm -f cvbandpre.h; ln -fs ../cvode/include/cvbandpre.h .)
	@(cd ../../include; rm -f cvdense.h;   ln -fs ../cvode/include/cvdense.h .)
	@(cd ../../include; rm -f cvdiag.h;    ln -fs ../cvode/include/cvdiag.h .)
	@(cd ../../include; rm -f cvode.h;     ln -fs ../cvode/include/cvode.h .)
	@(cd ../../include; rm -f cvspgmr.h;   ln -fs ../cvode/include/cvspgmr.h .)
	@(cd ../../include; rm -f cvbbdpre.h;  ln -fs ../cvode/include/cvbbdpre.h .)

purge:
	@(echo 'CVODE MODULE')
	@(echo '...Remove library file...')
	@(rm -f $(CVODE_LIB_DIR)/$(CVODE_LIB_NAME))
	@(echo '...Delete symbolic links...')
	@(cd ../../lib;     rm -f $(CVODE_LIB_NAME);)
	@(cd ../../include; rm -f cvband.h;)
	@(cd ../../include; rm -f cvbandpre.h;)
	@(cd ../../include; rm -f cvdense.h;)
	@(cd ../../include; rm -f cvdiag.h;)
	@(cd ../../include; rm -f cvode.h;)
	@(cd ../../include; rm -f cvspgmr.h;)
	@(cd ../../include; rm -f cvbbdpre.h;)

#---End of Makefile---