#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 22 March 2002
#--------------------------------------------------------------------------------
# Makefile for the CVODES adjoint sensitivity examples.
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
	@(echo 'ADJOINT SENSITIVITY EXAMPLES FOR CVODES')
	@(echo)
	@(echo 'Usage: make <ex>     - make example <ex>')
	@(if test -f cvax_ser.mk; \
		then echo '       make cvfx_ser - make all serial CVODES examples'; \
	fi)
	@(if test -f cvax_par.mk; \
		then echo '       make cvfx_par - make all parallel CVODES examples'; \
	fi)
	@(if (test -f cvax_ser.mk) && (test -f cvax_par.mk); \
		then echo '       make examples - make all CVODES examples'; \
	fi)
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(if test -f cvax_ser.mk; then make -f cvax_ser.mk; fi)
	@(echo)
	@(if test -f cvax_par.mk; then make -f cvax_par.mk; fi)
	@(echo)

examples: cvax_ser cvax_par

cvax_ser:
	@(if test -f cvax_ser.mk; then make -f cvax_ser.mk examples; fi)

cvax_par:
	@(if test -f cvax_par.mk; then make -f cvax_par.mk examples; fi)

cvadx cvabx cvakx cvakxb:
	@(if test -f cvax_ser.mk; then make -f cvax_ser.mk $@; fi)

pvanx:
	@(if test -f cvax_par.mk; then make -f cvax_par.mk $@; fi)

purge:
	@(if test -f cvax_ser.mk; then make -f cvax_ser.mk purge; fi)
	@(if test -f cvax_par.mk; then make -f cvax_par.mk purge; fi)

#---End of Makefile---
