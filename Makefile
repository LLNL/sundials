#--------------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
# Version of  : 8 March 2002
#--------------------------------------------------------------------------------------
# Makefile for SUNDIALS suite
#
# 'make install' recursively generates libraries for all existing SUNDIALS modules
#
# 'make purge'   recursively purges all existing SUNDIALS modules
#
#--------------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh


#--------------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------------

all:
	@(echo)
	@(echo 'Usage: make install - generate libraries for all existing SUNDIALS modules')
	@(echo '                      create symbolic links for all existing SUNDIALS modules')
	@(echo '       make purge   - remove libraries for all existing SUNDIALS modules')
	@(echo '                      remove symbolic links for all existing SUNDIALS modules')
	@(echo)

install: shared_lib cvode_lib ida_lib kinsol_lib cvodes_lib

purge: shared_purge cvode_purge ida_purge kinsol_purge cvode_purge

#--------------------------------------------------------------------------------

shared_lib:
	@if test -d shared; \
		then cd shared/source; make lib; make nvecserial; make nvecparallel; \
		else echo 'SHARED module does not exist!!! ERROR'; \
	fi

cvode_lib:
	@if test -d cvode; \
		then cd cvode/source; make lib; \
	fi

ida_lib:
	@if test -d ida; \
		then cd ida/source; make lib; \
	fi

kinsol_lib:
	@if test -d kinsol; \
		then cd kinsol/source; make lib; \
	fi
	@if test -d kinsol/fcmix; \
		then cd kinsol/fcmix; make lib; \
	fi

cvodes_lib:
	@if test -d cvodes; \
		then cd cvodes/source; make lib; \
	fi

#--------------------------------------------------------------------------------

shared_purge:
		@if test -d shared; \
		then cd shared/source; make purge; \
		else echo 'SHARED module does not exist!!! ERROR'; \
	fi

cvode_purge:
	@if test -d cvode; \
		then cd cvode/source; make purge; \
	fi

ida_purge:
	@if test -d ida; \
		then cd ida/source; make purge; \
	fi

kinsol_purge:
	@if test -d kinsol/fcmix; \
		then cd kinsol/fcmix; make purge; \
	fi
	@if test -d kinsol; \
		then cd kinsol/source; make purge; \
	fi

cvodes_purge:
	@if test -d cvodes; \
		then cd cvodes/source; make purge; \
	fi

#--------------------------------------------------------------------------------
