# File cvodex.mk.  Version of 14 December 2001
#--------------------------------------------------------------------------
# Makefile for the CVODE examples
#
# For each example requested, this makefile compiles the source, links,
# and removes the .o files.
#
# Usage:  make -f cvodex.mk   [to display list of examples]
#         make -f cvodex.mk <ex>   [to make example <ex>, where <ex> is 
#                                  cvdx, cvbx, cvkx, cvkxb, cvdemd, cvdemk]
#         make -f cvodex.mk examples   [to make all examples]
#         make -f cvodex.mk purge   [to remove all example executable files]
# 
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------


#==========================================================================
# Machine-dependent variables
#==========================================================================

#-----------
# C compiler
#-----------
CC          = gcc
CLINKER     = $(CC)


#======================================================
# Other variables
#======================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR	=	../../shared

#--------------
# Include files
#--------------
INC_SHARED	=	$(SHARED_DIR)/include
INC_SOLVER	=	../include

#---------------
# Compiler flags
#---------------
CFLAGS		=	-I$(INC_SOLVER) -I$(INC_SHARED)

#----------
# Libraries
#----------
LIB		=	../../lib/
LIB_LIST	=	-lcvode -lm
LIBS		=	-L$(LIB) $(LIB_LIST)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of CVODE examples:')
	@(echo '  cvdx:  dense example')
	@(echo '  cvbx:  banded example')
	@(echo '  cvkx:  Krylov example')
	@(echo '  cvkxb  Krylov example with banded preconditioner')
	@(echo '  cvdemd: demonstration program for direct methods')
	@(echo '  cvdemk: demonstration program for Krylov methods')

cvdx: 
	$(CC) $(CFLAGS) -c cvdx.c
	$(CLINKER) -o cvdx cvdx.o $(LIBS)
	@(rm -f cvdx.o)

cvbx: 
	$(CC) $(CFLAGS) -c cvbx.c
	$(CLINKER) -o cvbx cvbx.o $(LIBS)
	@(rm -f cvbx.o)

cvkx: 
	$(CC) $(CFLAGS) -c cvkx.c
	$(CLINKER) -o cvkx cvkx.o $(LIBS)
	@(rm -f cvkx.o)

cvkxb: 
	$(CC) $(CFLAGS) -c cvkxb.c
	$(CLINKER) -o cvkxb cvkxb.o $(LIBS)
	@(rm -f cvkxb.o)

cvdemd: 
	$(CC) $(CFLAGS) -c cvdemd.c
	$(CLINKER) -o cvdemd cvdemd.o $(LIBS)
	@(rm -f cvdemd.o)

cvdemk: 
	$(CC) $(CFLAGS) -c cvdemk.c
	$(CLINKER) -o cvdemk cvdemk.o $(LIBS)
	@(rm -f cvdemk.o)

examples: cvdx cvbx cvkx cvkxb cvdemd cvdemk

purge:
	@(rm -f cvdx)
	@(rm -f cvbx)
	@(rm -f cvkx)
	@(rm -f cvkxb)
	@(rm -f cvdemd)
	@(rm -f cvdemk)
