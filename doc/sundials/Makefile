# -----------------------------------------------------------------
# $Revision$
# $Date$
# -----------------------------------------------------------------
# Programmer(s): Carol Woodward @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see sundials/kinsol/LICENSE.
# -----------------------------------------------------------------
# Makefile for SUNDIALS Installation Guide
# -----------------------------------------------------------------

QUIET =

SHARED_DOC = ../sundials

UG_FILES = \
install_guide.tex

UG_SHARED_FILES = \
install.tex \
ug.tex \
biblio.bib \
texi2dvi 

UG_SHARED_FIGURES = \
doc_logo.pdf

UG_CMAKE_FIGURES = \
ccmakedefault.png \
ccmakeprefix.png \
testrunnertop.png\
testrunnerbot.png
#ccmakeempty.pdf \
#cmaketest.pdf

#-------------------------------------------------------------------------------

all:
	@echo "Usage:"
	@echo "  make in_g   - compile LaTeX files to generate PDF file for Instalation Guide"
	@echo "  make clean  - remove all temporary files"
	@echo "  make purge  - remove temporary files + PDF files"
	@echo ""

#-------------------------------------------------------------------------------

in_g: 	copy_ug_figs copy_ug_cmake_figs install_guide.pdf

copy_ug:
	@for i in ${UG_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

copy_ug_figs:
	@for i in ${UG_SHARED_FIGURES} ; do \
	cp ${SHARED_DOC}/figures/$${i} . ;\
	done

copy_ug_cmake_figs:
	@for i in ${UG_CMAKE_FIGURES} ; do \
	cp ${SHARED_DOC}/figures/cmake/$${i} . ;\
	done

install_guide.pdf: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	xelatex ${QUIET} install_guide.tex
	xelatex ${QUIET} install_guide.tex


#-------------------------------------------------------------------------------

clean:
#	rm -rf ${UG_SHARED_FILES}
	rm -rf ${UG_SHARED_FIGURES}
	rm -rf ${UG_CMAKE_FIGURES}
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out

purge:  clean
	rm -f install_guide.ps install_guide.pdf
