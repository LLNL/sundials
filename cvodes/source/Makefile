#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 7 March 2002
#--------------------------------------------------------------------------------
# Makefile for the CVODES package.
#
# 'make lib' generate the CVODES library (libcvodes.a) in sundials/cvode/lib
#            add a symbolic link to the CVODES library in sundials/lib
#            and symbolic links to the CVODES header files in sundials/include
#
# 'make purge' remove CVODES library and symbolic links
#--------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../../shared/include
CVODES_INC_DIR  = ../include

# CVODE library name and location
ARCH           = `uname -s`.`uname -m`
CVODES_LIB_DIR  = ../lib
CVODES_LIB_NAME = libcvodes.$(ARCH).a

# Source and object files
SRC    = cvsband.c     cvsdense.c  \
         cvsdiag.c     cvsspgmr.c  \
         cvodes.c      cvodea.c    \
         cvsbandpre.c  cvsbbdpre.c

OBJ = ${SRC:.c=.o}

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(CVODES_INC_DIR) -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------


all:
	@(echo)
	@(echo 'Usage: make lib   - generate the CVODES package library')
	@(echo '                    add symbolic link to CVODES library')
	@(echo '                    add symbolic links to CVODES header files')
	@(echo)
	@(echo '       make purge - remove CVODES library and symbolic links')
	@(echo)

lib:
	@(echo 'CVODES MODULE')
	@(echo '...Compile object files...')
	@($(CC) $(CFLAGS) -c $(SRC))
	@(echo '...Create library file...')
	@(ar rc $(CVODES_LIB_DIR)/$(CVODES_LIB_NAME) $(OBJ))
	@(rm -f $(OBJ))
	@(echo '...Create symbolic links...')
	@(cd ../../lib;     rm -f $(CVODES_LIB_NAME);  ln -fs ../cvodes/lib/$(CVODES_LIB_NAME) .)
	@(cd ../../include; rm -f cvsband.h;    ln -fs ../cvodes/include/cvsband.h .)
	@(cd ../../include; rm -f cvsbandpre.h; ln -fs ../cvodes/include/cvsbandpre.h .)
	@(cd ../../include; rm -f cvsdense.h;   ln -fs ../cvodes/include/cvsdense.h .)
	@(cd ../../include; rm -f cvsdiag.h;    ln -fs ../cvodes/include/cvsdiag.h .)
	@(cd ../../include; rm -f cvodea.h;     ln -fs ../cvodes/include/cvodea.h .)
	@(cd ../../include; rm -f cvodes.h;     ln -fs ../cvodes/include/cvodes.h .)
	@(cd ../../include; rm -f cvsspgmr.h;   ln -fs ../cvodes/include/cvsspgmr.h .)
	@(cd ../../include; rm -f cvsbbdpre.h;  ln -fs ../cvodes/include/cvsbbdpre.h .)

purge:
	@(echo 'CVODES MODULE')
	@(echo '...Remove library file...')
	@(rm -f $(CVODES_LIB_DIR)/$(CVODES_LIB_NAME))
	@(echo '...Delete symbolic links...')
	@(cd ../../lib;     rm -f $(CVODES_LIB_NAME);)
	@(cd ../../include; rm -f cvsband.h;)
	@(cd ../../include; rm -f cvsbandpre.h;)
	@(cd ../../include; rm -f cvsdense.h;)
	@(cd ../../include; rm -f cvsdiag.h;)
	@(cd ../../include; rm -f cvodea.h;)
	@(cd ../../include; rm -f cvodes.h;)
	@(cd ../../include; rm -f cvsspgmr.h;)
	@(cd ../../include; rm -f cvsbbdpre.h;)

#---End of Makefile---

