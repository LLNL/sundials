#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 28 August 2002
#--------------------------------------------------------------------------
# Makefile for FKINSOL, the Fortran interface package for KINSOL and the
# associated examples.
#
#
# Type 'make' to see usage instructions.
#--------------------------------------------------------------------------

SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../../shared/include
KINSOL_INC_DIR  = ../include

# KINSOL library name and location
ARCH           = `uname -s`.`uname -m`
KINSOL_LIB_DIR  = ../lib
FKINSOL_LIB_NAME = libfkinsol.$(ARCH).a

# Source and object files
SRC = fkinsol.c fkinspgmr01.c fkinspgmr10.c fkinspgmr11.c \
      fkinspgmr20.c fkinspgmr21.c fkinpreco.c fkinpsol.c \
      fkinuatimes.c fkinbbd.c fkinbbdinit1.c

OBJ = ${SRC:.c=.o}

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(KINSOL_INC_DIR) -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all:
	@(echo)
	@(echo 'Makefile for FKINSOL module')
	@(echo)
	@(echo 'Usage: make lib       - create Fortran interface KINSOL library')
	@(echo)
	@(echo '       make purge     - remove Fortran interface KINSOL library')
	@(echo)
	@(echo '       make <ex>      - make example <ex>')
	@(if test -f fkinx_ser.mk; \
		then echo '       make fkinx_ser - make all serial FKINSOL examples'; \
	fi)
	@(if test -f fkinx_par.mk; \
		then echo '       make fkinx_par - make all parallel FKINSOL examples'; \
	fi)
	@(if (test -f fkinx_ser.mk) && (test -f fkinx_par.mk); \
		then echo '       make examples  - make all FKINSOL examples'; \
	fi)
	@(echo '       make purgex    - remove example executables')
	@(echo)
	@(if test -f fkinx_ser.mk; then make -f fkinx_ser.mk; fi)
	@(echo)
	@(if test -f fkinx_par.mk; then make -f fkinx_par.mk; fi)
	@(echo)

#-----------------------------------------------------------------
# FKINSOL library

lib:
	@(echo 'FKINSOL MODULE')
	@(echo '...Compile object files...')
	@($(CC) $(CFLAGS) -c $(SRC))
	@(echo '...Create library file...')
	@(ar rc $(KINSOL_LIB_DIR)/$(FKINSOL_LIB_NAME) $(OBJ))
	@(rm -f $(OBJ))
	@(echo '...Create symbolic links...')
	@(cd ../../lib; rm -f $(FKINSOL_LIB_NAME); ln -s ../kinsol/lib/$(FKINSOL_LIB_NAME) .)

purge:
	@(echo 'FKINSOL MODULE')
	@(echo '...Remove library file...')
	@(rm -f $(KINSOL_LIB_DIR)/$(FKINSOL_LIB_NAME))
	@(cd ../../lib; rm -f $(FKINSOL_LIB_NAME))

#-----------------------------------------------------------------
# FKINSOL Examples

examples: fkinx_ser fkinx_par

fkinx_ser:
	@(if test -f fkinx_ser.mk; then make -f fkinx_ser.mk examples; fi)

fkinx_par:
	@(if test -f fkinx_par.mk; then make -f fkinx_par.mk examples; fi)

kindiagsf:
	@(make -f fkinx_ser.mk $@)

kindiagpf:
	@(make -f fkinx_par.mk $@)

purgex:
	@(if test -f fkinx_ser.mk; then make -f fkinx_ser.mk purge; fi)
	@(if test -f fkinx_par.mk; then make -f fkinx_par.mk purge; fi)

#---End of Makefile---
