# File kinxs.mk.  Version of 13 December 2001
#--------------------------------------------------------------------------
# Makefile for the serial KINSOL example
#
# For the example, this makefile compiles the source, links,
# and removes the .o file.
#
# Usage:  make -f kinxs.mk   [to display list of examples]
#         make -f kinxs.mk kinwebs   [to make example kinwebs]
#         make -f kinxs.mk examples   [to make all examples]
#         make -f kinxs.mk purge   [to remove all example executable files]
# 
# Check, and modify if necessary, the machine-dependent variables below.
#--------------------------------------------------------------------------


#==========================================================================
# Machine-dependent variables
#==========================================================================

#-----------
# C compiler
#-----------
CC      = gcc
CLINKER = $(CC)


#======================================================
# Other variables
#======================================================

#---------------------
# Path to shared files
#---------------------
SHARED_DIR = ../../shared

#--------------
# Include files
#--------------
INC_SHARED = $(SHARED_DIR)/include
INC_SOLVER = ../include

#---------------
# Compiler flags
#---------------
CFLAGS     = -I$(INC_SOLVER) -I$(INC_SHARED)

#----------
# Libraries
#----------
LIB        = ../../lib/
LIB_LIST   = -lkinsols -lm
LIBS       = -L$(LIB) $(LIB_LIST)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of serial KINSOL examples:')
	@(echo '  kinwebs:  2-D food web system, block-diagonal preconditioner')

kinwebs: 
	$(CC) $(CFLAGS) -c kinwebs.c
	$(CLINKER) -o kinwebs kinwebs.o $(LIBS)
	@(rm -f kinwebs.o)

examples: kinwebs

purge:
	@(rm -f kinwebs)
