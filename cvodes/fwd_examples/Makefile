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
	@(echo '       make cvfx_ser - make all serial CVODES examples')
	@(echo '       make cvfx_par - make all parallel CVODES examples')
	@(echo '       make examples - make all CVODES examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(if test -f cvfx_ser.mk; then make -f cvfx_ser.mk; fi)
	@(echo)
	@(if test -f cvfx_par.mk; then make -f cvfx_par.mk; fi)
	@(echo)

examples: cvfx_ser cvfx_par

cvfx_ser:
	@(if test -f cvfx_ser.mk; then make -f cvfx_ser.mk examples; fi)

cvfx_par:
	@(if test -f cvfx_par.mk; then make -f cvfx_par.mk examples; fi)

cvfdx cvfkx cvfnx:
	@(if test -f cvfx_ser.mk; then make -f cvfx_ser.mk $@; fi)

pvfkx pvfnx:
	@(if test -f cvfx_par.mk; then make -f cvfx_par.mk $@; fi)

purge:
	@(if test -f cvfx_ser.mk; then make -f cvfx_ser.mk purge; fi)
	@(if test -f cvfx_par.mk; then make -f cvfx_par.mk purge; fi)

#---End of Makefile---