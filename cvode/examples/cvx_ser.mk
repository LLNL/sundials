# File cvx_ser.mk.  Version of 5 March 2002
#--------------------------------------------------------------------------
# Makefile for the CVODE examples using the serial N_Vector implementation.
#
# For each example requested, this makefile compiles the source and
# generates the executable.
#
# Usage:  make -f cvx_ser.mk          [to display list of examples]
#         make -f cvx_ser.mk <ex>     [to make example <ex>, where <ex> is 
#                                     cvdx, cvbx, cvkx, cvkxb, cvdemd, cvdemk]
#         make -f cvx_ser.mk examples [to make all examples]
#         make -f cvx_ser.mk purge    [to remove all example executable files]
# 
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------

SHELL = /bin/sh

#--------------------------
# Top of SUNDIALS directory
#--------------------------
SUNDIALS_DIR = ../..

#---------------------
# Path to header files
#---------------------
INC_DIR	= $(SUNDIALS_DIR)/include

#----------------------
# Path to library files
#----------------------
LIB_DIR = $(SUNDIALS_DIR)/lib

#-------------
# Architecture
#-------------
ARCH = `uname -s`.`uname -m`

#======================================================
# Machine-dependent variables
#======================================================

#------------------------------
# C compiler and compiler flags
#------------------------------

CC     = gcc
CFLAGS = -Wall -ffloat-store -I$(INC_DIR) -L$(LIB_DIR)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of serial CVODE examples (using the serial NVECTOR module):')
	@(echo '  cvdx  :  dense example')
	@(echo '  cvbx  :  banded example')
	@(echo '  cvkx  :  Krylov example')
	@(echo '  cvkxb :  Krylov example with banded preconditioner')
	@(echo '  cvdemd: demonstration program for direct methods')
	@(echo '  cvdemk: demonstration program for Krylov methods')

cvdx: cvdx.c
	@echo '...Compile cvdx...'
	@$(CC) $(CFLAGS) -o cvdx cvdx.c -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvdx.o

cvbx: cvbx.c
	@echo '...Compile cvbx...'
	@$(CC) $(CFLAGS) -o cvbx cvbx.c -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvbx.o

cvkx: cvkx.c
	@echo '...Compile cvkx...'
	@$(CC) $(CFLAGS) -o cvkx cvkx.c -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvkx.o

cvkxb: cvkxb.c
	@echo '...Compile cvkxb...'
	@$(CC) $(CFLAGS) -o cvkxb cvkxb.c -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvkxb.o

cvdemd: cvdemd.c
	@echo '...Compile cvdemd...'
	@$(CC) $(CFLAGS) -o cvdemd cvdemd.c -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvdemd.o

cvdemk: cvdemk.c
	@echo '...Compile cvdemk...'
	@$(CC) $(CFLAGS) -o cvdemk cvdemk.c -lcvode.$(ARCH) -lshared.$(ARCH) -lnvecserial.$(ARCH) -lm
	@rm -f cvdemk.o

examples: cvdx cvbx cvkx cvkxb cvdemd cvdemk

purge:
	@(rm -f cvdx)
	@(rm -f cvbx)
	@(rm -f cvkx)
	@(rm -f cvkxb)
	@(rm -f cvdemd)
	@(rm -f cvdemk)

#---End of cvx_ser.mk---