#File Makefile.diagsf  Version of 03 Dec 1998

# Makefile for diagsf demo 


# SUN Sparc version:

CC          = cc
CLINKER     = cc
F77         = f77

KINLIB     = ../lib/libkinsols.a
LIBS       =  -lm  -lF77 -lsunmath_mt

KNC = ../include

CFLAGS = 

.c.o:
	$(CC) $(CFLAGS) -c $<

diagsf: diagsf.o 
	$(CLINKER) -o diagsf diagsf.o $(KINLIB) $(LIBS)

diagsf.o: diagsf.f
	$(F77) -c  diagsf.f


