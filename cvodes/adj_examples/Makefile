#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 22 March 2002
#--------------------------------------------------------------------------------
# Makefile for the adjoint CVODES examples.
#
# This makefile generates the executables for the adjoint CVODES examples requested.
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
	@(echo '       make cvax_ser - make all serial adjoint CVODES examples')
	@(echo '       make cvax_par - make all parallel adjoint CVODES examples')
	@(echo '       make examples - make all adjoint CVODES examples')
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(if test -f cvax_ser.mk; then make -f cvax_ser.mk; fi)
	@(echo)
	@(if test -f cvax_par.mk; then make -f cvax_par.mk; fi)
	@(echo)

cvadx cvabx cvakx cvakxb:
	@(if test -f cvax_ser.mk; then make -f cvax_ser.mk $@; fi)

pvanx:
	@(if test -f cvax_par.mk; then make -f cvax_par.mk $@; fi)

purge:
	@(if test -f cvax_ser.mk; then make -f cvax_ser.mk purge; fi)
	@(if test -f cvax_par.mk; then make -f cvax_par.mk purge; fi)

#---End of Makefile---
