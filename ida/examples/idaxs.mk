# File idaxs.mk.  Version of 17 December 2001
#--------------------------------------------------------------------------
# Makefile for the serial IDA example
#
# For each example, this makefile compiles the source, links,
# and removes the .o file.
#
# Usage:  make -f idaxs.mk   [to display list of examples]
#         make -f idaxs.mk <ex>    [where <ex> is irobx, iheatsb, iheatsk,
#                                   or iwebsb]
#         make -f idaxs.mk examples   [to make all examples]
#         make -f idaxs.mk purge   [to remove all example executable files]
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
LIB      = ../../lib/
LIB_LIST = -lidas -lm
LIBS     = -L$(LIB) $(LIB_LIST)


#======================================================
# Make rules
#======================================================

all:
	@(echo 'List of serial IDA examples:')
	@(echo '  irobx:    3-species Robertson kinetics system')
	@(echo '  iheatsb:  2-D heat equation, banded Jacobian')
	@(echo '  iheatsk:  2-D heat equation, diagonal preconditioner')
	@(echo '  iwebsb:   2-D food web system, banded Jacobian')

irobx:
	$(CC) $(CFLAGS) -c irobx.c
	$(CLINKER) -o irobx irobx.o $(LIBS)
	@(rm -f irobx.o)

iheatsb:
	$(CC) $(CFLAGS) -c iheatsb.c
	$(CLINKER) -o iheatsb iheatsb.o $(LIBS)
	@(rm -f iheatsb.o)

iheatsk:
	$(CC) $(CFLAGS) -c iheatsk.c
	$(CLINKER) -o iheatsk iheatsk.o $(LIBS)
	@(rm -f iheatsk.o)

iwebsb:
	$(CC) $(CFLAGS) -c iwebsb.c
	$(CLINKER) -o iwebsb iwebsb.o $(LIBS)
	@(rm -f iwebsb.o)

examples: irobx iheatsb iheatsk iwebsb

purge:
	@(rm -f irobx iheatsb iheatsk iwebsb)
