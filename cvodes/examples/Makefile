#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 20 March 2002
#--------------------------------------------------------------------------------
# Makefile for the CVODES examples.
#
# This makefile generates the executables for the CVODES examples requested.
# 
# Type 'make' to see usage instructions and list of examples.
#
#--------------------------------------------------------------------------------

SHELL = /bin/sh

#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo 'Usage: make <ex>     - make example <ex>')
	@(echo '       make cvsx_ser - make all serial CVODES examples')
	@(echo '       make cvsx_par - make all parallel CVODES examples')
	@(echo '       make examples - make all CVODES examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(if test -f cvsx_ser.mk; then make -f cvsx_ser.mk; fi)
	@(echo)
	@(if test -f cvsx_par.mk; then make -f cvsx_par.mk; fi)
	@(echo)

examples: cvsx_ser cvsx_par

cvsx_ser:
	@(if test -f cvsx_ser.mk; then make -f cvsx_ser.mk examples; fi)

cvsx_par:
	@(if test -f cvsx_par.mk; then make -f cvsx_par.mk examples; fi)

cvsdx cvskx cvsnx:
	@(if test -f cvsx_ser.mk; then make -f cvsx_ser.mk $@; fi)

pvskx pvsnx:
	@(if test -f cvsx_par.mk; then make -f cvsx_par.mk $@; fi)

purge:
	@(if test -f cvsx_ser.mk; then make -f cvsx_ser.mk purge; fi)
	@(if test -f cvsx_par.mk; then make -f cvsx_par.mk purge; fi)

#---End of Makefile---