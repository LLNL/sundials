# File SensMakefile.idas  Version of 20 March 2001

# Makefile for the serial version of the SensIDA package

# Notes on machine and environment dependencies:
#
# (1) Files from the SensIDA solver are accessed as follows:
#     (a) All SensIDA header files needed are assumed to be in the
#         directory given by SensIDAINC below.
#     (b) The SensIDA library file is given by SensIDALIB below.
#
# (2) The C compiler is given by CC below.
#
# Change these variables as needed for other environments.


SensIDAINC = ../include
SensIDALIB = ../lib/libidas.a

# SUN Sparc version:
CC          = cc -g

# Source files:

CODES = ida.c llnlmath.c idaspgmr.c spgmr.c iterativ.c idadense.c dense.c \
        idaband.c band.c nvector.c sensida.c sensidaspgmr.c sensidaband.c \
	sensidadense.c

# Object files:

OBJS =  ida.o llnlmath.o idaspgmr.o spgmr.o iterativ.o idadense.o dense.o \
        idaband.o band.o nvector.o sensida.o sensidaspgmr.o sensidaband.o \
	sensidadense.o


# Command sequence:
# Note that the proper versions of the files nvector.c and nvector.h are
# copied into the proper directories for the compile, then removed.
# However, leave the appropriate version of nvector.h in /include .
# This is needed by the driver programs in /example .

all:
	cp ./serial/sens_nvector.c nvector.c
	cp $(SensIDAINC)/serial/sens_nvector.h $(SensIDAINC)/nvector.h
	cp $(SensIDAINC)/serial/sens_nvector.h $(SensIDAINC)/serial/nvector.h
	$(CC) -g -c $(CODES) -I$(SensIDAINC)
	ar rcv $(SensIDALIB) $(OBJS)
	rm -f *.o
#	rm $(SensIDAINC)/nvector.h
#	rm ./nvector.c
