#File Makefile.senskinsols  Version of 21 August 2000
#
# Makefile for library libsenskinsols.a for the Sens_KINSOL package
#
# Notes on machine and environment:
#
# (1) Files from the KINSOL solver are accessed as follows:
#     (a) The KINSOL library file is assumed to be
#	  given by KINLIB below.
#     (b) All KINSOL header files needed are assumed to be in the
#         directory given by KININC below.
#
# (2) The C compiler is given by CC below.
#
# (3) Environment-specific compiler flags must be included in CFLAGS below.
# 
# Change these variables as needed for other environments.
# 

#serial case

KINLIB = ../lib/libsens_kinsols.a
KININC= ../include

CC      = cc

# source files

CODES = sens_kinsol.c sens_kinspgmr.c  kinsol.c llnlmath.c kinspgmr.c \
	spgmr.c iterativ.c smalldense.c band.c nvector.c

# object files

OBJS = sens_kinsol.o sens_kinspgmr.o kinsol.o llnlmath.o kinspgmr.o \
	spgmr.o iterativ.o smalldense.o band.o nvector.o


all:
	cp serial/nvector.c nvector.c
	cp ../include/serial/nvector.h ../include/nvector.h
	$(CC) -c $(CODES) -I$(KININC)
	ar rcv $(KINLIB) $(OBJS) 
	rm -f *.o

