#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 27 February 2002
#--------------------------------------------------------------------------------
# Makefile for the CVODE examples.
#
# This makefile generates the executables for the CVODE examples requested.
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
	@(if test -f cvx_ser.mk; \
		then echo '       make cvx_ser  - make all serial CVODE examples'; \
	fi)
	@(if test -f cvx_par.mk; \
		then echo '       make cvx_par  - make all parallel CVODE examples'; \
	fi)
	@(if (test -f cvx_ser.mk) && (test -f cvx_par.mk); \
		then echo '       make examples - make all CVODE examples'; \
	fi)
	@(echo '       make purge    - remove all example executable files')
	@(echo)
	@(if test -f cvx_ser.mk; then make -f cvx_ser.mk; fi)
	@(echo)
	@(if test -f cvx_par.mk; then make -f cvx_par.mk; fi)
	@(echo)

examples: cvx_ser cvx_par

cvx_ser:
	@(if test -f cvx_ser.mk; then make -f cvx_ser.mk examples; fi)

cvx_par:
	@(if test -f cvx_par.mk; then make -f cvx_par.mk examples; fi)

cvbx cvdx cvdemd cvkx cvkxb cvdemk:
	@(if test -f cvx_ser.mk; then make -f cvx_ser.mk $@; fi)

pvnx pvkx pvkxb:
	@(if test -f cvx_par.mk; then make -f cvx_par.mk $@; fi)

purge:
	@(if test -f cvx_ser.mk; then make -f cvx_ser.mk purge; fi)
	@(if test -f cvx_par.mk; then make -f cvx_par.mk purge; fi)

#---End of Makefile---