#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 7 March 2002
#--------------------------------------------------------------------------
# Makefile for the KINSOL examples (serial and parallel).
#
# This makefile generates the executables for the KINSOL examples requested.
#
# Type 'make' to see usage instructions and list of examples.
#--------------------------------------------------------------------------

SHELL = /bin/sh

#======================================================
# Make rules
#======================================================

all:
	@(echo)
	@(echo 'Usage: make <ex>     - make example <ex>')
	@(echo '       make kinx_ser - make all serial KINSOL examples')
	@(echo '       make kinx_par - make all parallel KINSOL examples')
	@(echo '       make examples - make all KINSOL examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(if test -f kinx_ser.mk; then make -f kinx_ser.mk; fi)
	@(echo)
	@(if test -f kinx_par.mk; then make -f kinx_par.mk; fi)
	@(echo)

examples: kinx_ser kinx_par

kinx_ser:
	@(if test -f kinx_ser.mk; then make -f kinx_ser.mk examples; fi)

kinx_par:
	@(if test -f kinx_par.mk; then make -f kinx_par.mk examples; fi)

kinwebs:
	@(make -f kinx_ser.mk $@)

kinwebp kinwebbbd:
	@(make -f kinx_par.mk $@)

purge:
	@(if test -f kinx_ser.mk; then make -f kinx_ser.mk purge; fi)
	@(if test -f kinx_par.mk; then make -f kinx_par.mk purge; fi)

#---End of Makefile---