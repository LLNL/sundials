#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 11 March 2002
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
KINSOL_LIB_NAME = libkinsol.$(ARCH).a

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
	@(echo 'Usage: make lib       - append Fortran interface to KINSOL library')
	@(echo)
	@(echo '       make purge     - remove Fortran interface from KINSOL library')
	@(echo)
	@(echo '       make <ex>      - make example <ex>')
	@(echo '       make fkinx_ser - make all serial FKINSOL examples')
	@(echo '       make fkinx_par - make all parallel FKINSOL examples')
	@(echo '       make examples  - make all FKINSOL examples')
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
	@(if test -f $(KINSOL_LIB_DIR)/$(KINSOL_LIB_NAME); \
		then \
			echo '...Compile object files...'; \
			$(CC) $(CFLAGS) -c $(SRC); \
			echo '...Append to KINSOL library...'; \
			ar rc $(KINSOL_LIB_DIR)/$(KINSOL_LIB_NAME) $(OBJ); \
			rm -f $(OBJ); \
		else \
			echo 'ERROR: Library KINSOL does not exist'; \
			echo '       Make the KINSOL library first'; \
	fi)

purge:
	@(echo 'FKINSOL MODULE')
	@(if test -f $(KINSOL_LIB_DIR)/$(KINSOL_LIB_NAME); \
		then \
			echo '...Remove fkinsol objects from kinsol library...'; \
			ar d $(KINSOL_LIB_DIR)/$(KINSOL_LIB_NAME) $(OBJ); \
	fi)

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
