#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 8 March 2002
#--------------------------------------------------------------------------------
# Makefile for the IDA examples.
#
# This makefile generates the executables for the IDA examples requested.
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
	@(echo '       make idax_ser - make all serial IDA examples')
	@(echo '       make idax_par - make all parallel IDA examples')
	@(echo '       make examples - make all IDA examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(if test -f idax_ser.mk; then make -f idax_ser.mk; fi)
	@(echo)
	@(if test -f idax_par.mk; then make -f idax_par.mk; fi)
	@(echo)

examples: idax_ser idax_par

idax_ser:
	@(if test -f idax_ser.mk; then make -f idax_ser.mk examples; fi)

idax_par:
	@(if test -f idax_par.mk; then make -f idax_par.mk examples; fi)

irobx iheatsb iheatsk iwebsb:
	@(if test -f idax_ser.mk; then make -f idax_ser.mk $@; fi)

iheatpk iheatbbd iwebpk iwebbbd:
	@(if test -f idax_par.mk; then make -f idax_par.mk $@; fi)

purge:
	@(if test -f idax_ser.mk; then make -f idax_ser.mk purge; fi)
	@(if test -f idax_par.mk; then make -f idax_par.mk purge; fi)

#---End of Makefile---