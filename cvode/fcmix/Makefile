#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 27 March 2002
#--------------------------------------------------------------------------
# Makefile for FCVODE, the Fortran interface package for CVODE and the
# associated examples.
#
#
# Type 'make' to see usage instructions.
#--------------------------------------------------------------------------

SHELL = /bin/sh

# Paths to header files
SHARED_INC_DIR = ../../shared/include
CVODE_INC_DIR  = ../include

# CVODE library name and location
ARCH           = `uname -s`.`uname -m`
CVODE_LIB_DIR  = ../lib
CVODE_LIB_NAME = libcvode.$(ARCH).a

# Source and object files
SRC = fcvode.c fcvpreco.c fcvpsol.c fcvjtimes.c fcvbbd.c fcvbbdin1.c \
      fcvdense1.c fcvband1.c \
      fcvspgmr01.c fcvspgmr11.c fcvspgmr21.c fcvspgmr10.c fcvspgmr20.c

OBJ = ${SRC:.c=.o}

# Compiler and compiler options
CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(CVODE_INC_DIR) -I$(SHARED_INC_DIR)

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all:
	@(echo)
	@(echo 'Makefile for FCVODE module')
	@(echo)
	@(echo 'Usage: make lib       - append Fortran interface to CVODE library')
	@(echo)
	@(echo '       make purge     - remove Fortran interface from CVODE library')
	@(echo)
	@(echo '       make <ex>      - make example <ex>')
	@(echo '       make fcvx_ser  - make all serial FCVODE examples')
	@(echo '       make fcvx_par  - make all parallel FCVODE examples')
	@(echo '       make examples  - make all FCVODE examples')
	@(echo '       make purgex    - remove example executables')
	@(echo)
	@(if test -f fcvx_ser.mk; then make -f fcvx_ser.mk; fi)
	@(echo)
	@(if test -f fcvx_par.mk; then make -f fcvx_par.mk; fi)
	@(echo)

#-----------------------------------------------------------------
# FCVODE library

lib:
	@(echo 'FCVODE MODULE')
	@(if test -f $(CVODE_LIB_DIR)/$(CVODE_LIB_NAME); \
		then \
			echo '...Compile object files...'; \
			$(CC) $(CFLAGS) -c $(SRC); \
			echo '...Append to CVODE library...'; \
			ar rc $(CVODE_LIB_DIR)/$(CVODE_LIB_NAME) $(OBJ); \
			rm -f $(OBJ); \
		else \
			echo 'ERROR: Library CVODE does not exist'; \
			echo '       Make the CVODE library first'; \
	fi)

purge:
	@(echo 'FCVODE MODULE')
	@(if test -f $(CVODe_LIB_DIR)/$(CVODE_LIB_NAME); \
		then \
			echo '...Remove fcvode objects from cvode library...'; \
			ar d $(CVODE_LIB_DIR)/$(CVODE_LIB_NAME) $(OBJ); \
	fi)

#-----------------------------------------------------------------
# FCVODE Examples

examples: fcvx_ser fcvx_par

fcvx_ser:
	@(if test -f fcvx_ser.mk; then make -f fcvx_ser.mk examples; fi)

fcvx_par:
	@(if test -f fcvx_par.mk; then make -f fcvx_par.mk examples; fi)

cvdensef:
	@(make -f fcvx_ser.mk $@)

pvdiagnf pvdiagkf pvdiagkbf:
	@(make -f fcvx_par.mk $@)

purgex:
	@(if test -f fcvx_ser.mk; then make -f fcvx_ser.mk purge; fi)
	@(if test -f fcvx_par.mk; then make -f fcvx_par.mk purge; fi)

#---End of Makefile---