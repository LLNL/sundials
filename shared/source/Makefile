#--------------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 26 June 2002
#--------------------------------------------------------------------------------------
# Makefile for the SHARED package.
#
# 'make lib'   generates the SHARED library in sundials/shared/lib
#              add symbolic link to the SHARED library in sundials/lib
#              add symbolic links to the SHARED header files in sundials/include
#
# 'make nvec***'  generate the *** NVECTOR library in sundials/shared/lib
#                 add symbolic link to the *** NVECTOR library in sundials/lib 
#                 add symbolic link to nvector_***.h in sundials/include
#
# 'make purge' remove SHARED library and symbolic links
#              invoke nvector makefiles to remove nvector libraries
#
#  Note: SUNDIALS provides the following two NVECTOR modules: nvecserial and nvecparallel.
#        Each NVECTOR module comes with its own makefile (nvec***.mk) 
#        which is self-contained and is called by this top level Makefile.
#
#--------------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# Path to header files
SHARED_INC_DIR = ../include

# SHARED library name and location
ARCH            = `uname -s`.`uname -m`
SHARED_LIB_DIR  = ../lib
SHARED_LIB_NAME = libshared.$(ARCH).a

# Source and object files
SRC = band.c dense.c iterativ.c sundialsmath.c \
      nvector.c smalldense.c spgmr.c

OBJ = ${SRC:.c=.o}

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------------

all:
	@(echo)
	@(echo 'Usage: make lib          - generate the SHARED library')
	@(echo '                           add symbolic link to libshared.a')
	@(echo '                           add symbolic links to SHARED header files')
	@(echo)
	@if test -f nvecserial.mk; \
		then echo '       make nvecserial   - generate the serial NVECTOR library'; \
		     echo '                           add symbolic link to libnvecserial.a'; \
		     echo '                           add symbolic link to nvector_serial.h' ; \
		     echo ; \
	fi
	@if test -f nvecparallel.mk; \
		then echo '       make nvecparallel - generate the parallel NVECTOR library'; \
		     echo '                           add symbolic link to libnvecparallel.a'; \
		     echo '                           add symbolic link to nvector_parallel.h'; \
		     echo; \
	fi
	@(echo '       make purge        - remove libraries and symbolic links')
	@(echo)

lib:
	@(echo 'SHARED MODULE')
	@(echo '...Compile SHARED object files...')
	@($(CC) $(CFLAGS) -c $(SRC))
	@(echo '...Create SHARED library file...')
	@(ar rc $(SHARED_LIB_DIR)/$(SHARED_LIB_NAME) $(OBJ))
	@(rm -f $(OBJ))
	@(echo '...Create symbolic links to SHARED...')
	@(cd ../../lib;     rm -f $(SHARED_LIB_NAME); ln -fs ../shared/lib/$(SHARED_LIB_NAME) .)
	@(cd ../../include; rm -f band.h;             ln -fs ../shared/include/band.h .)
	@(cd ../../include; rm -f dense.h;            ln -fs ../shared/include/dense.h .)
	@(cd ../../include; rm -f fcmixpar.h;         ln -fs ../shared/include/fcmixpar.h .)
	@(cd ../../include; rm -f iterativ.h;         ln -fs ../shared/include/iterativ.h .)
	@(cd ../../include; rm -f sundialsmath.h;     ln -fs ../shared/include/sundialsmath.h .)
	@(cd ../../include; rm -f sundialstypes.h;    ln -fs ../shared/include/sundialstypes.h .)
	@(cd ../../include; rm -f nvector.h;          ln -fs ../shared/include/nvector.h .)
	@(cd ../../include; rm -f smalldense.h;       ln -fs ../shared/include/smalldense.h .)
	@(cd ../../include; rm -f spgmr.h;            ln -fs ../shared/include/spgmr.h .)

nvecserial:
	@if test -f nvecserial.mk; \
		then make -f nvecserial.mk lib; \
	fi

nvecparallel:
	@if test -f nvecparallel.mk; \
		then make -f nvecparallel.mk lib; \
	fi

purge:
	@(echo 'SHARED MODULE')
	@(echo '...Remove library file...')
	@(rm -f $(SHARED_LIB_DIR)/$(SHARED_LIB_NAME))
	@(echo '...Delete symbolic links...')
	@(cd ../../lib;     rm -f $(SHARED_LIB_NAME);)
	@(cd ../../include; rm -f band.h;)
	@(cd ../../include; rm -f dense.h;)
	@(cd ../../include; rm -f fcmixpar.h;)
	@(cd ../../include; rm -f iterativ.h;)
	@(cd ../../include; rm -f sundialsmath.h;)
	@(cd ../../include; rm -f sundialstypes.h;)
	@(cd ../../include; rm -f nvector.h;)
	@(cd ../../include; rm -f smalldense.h;)
	@(cd ../../include; rm -f spgmr.h;)
	@(if test -f nvecserial.mk;  then make -f nvecserial.mk purge;   fi)
	@(if test -f nvecparallel.mk; then make -f nvecparallel.mk purge; fi)

#---End of Makefile---
