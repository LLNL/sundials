# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# File Makefile.fsenskinsols  Version of 07 September 2000 
#
# Makefile for fsenskinsols 
#
# Notes on machine and environment dependencies:
#
# (1) Files from the SensKINSOL solver are accessed as follows:
#     (a) All SensKINSOL header files needed are assumed to be in the
#         directory given by KNC below.
#     (b) The FSENSKINSOL library file libsenskinsols.a is assumed to be
#	  given by KINLIB below.
#
# (2) The C compiler is given by CC below.
#
# (3) Environment-specific compiler flags must be included in CFLAGS below.
# 
# Change these variables as needed for other environments.
#
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 

all: makelib

# Version for the CASC Suns:
CC = cc


# parallel or serial version
VECTOR = serial


KNC = ../include
KINLIB = ../lib/libsens_kinsols.a
CFLAGS = -I$(KNC)/$(VECTOR) -I$(KNC) 

# List of object files needed

OBJS = fsenskinsols.o fsenskinspgmr01.o fsenskinspgmr10.o \
       fsenskinspgmr11.o fsenskinspgmr20.o fsenskinspgmr21.o \
       fkinsols.o fkinspgmr01.o fkinspgmr10.o fkinspgmr11.o  \
       fkinspgmr20.o fkinspgmr21.o fkinpreco.o fkinpsol.o \
       fkinuatimes.o


makelib: clearvect $(OBJS)
	ar rcv $(KINLIB) $(OBJS)
	rm *.o


.c.o:
	$(CC) $(CFLAGS)  -c $*.c


clearvect:
	@if test -f $(KNC)/nvector.h ; then rm $(KNC)/nvector.h ; fi
