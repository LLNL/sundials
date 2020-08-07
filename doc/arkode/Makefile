# -----------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2020, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# Builds user guide and examples documents from Sphinx generated
# LaTeX files. To generate new LaTeX files run make latexpdf in the
# guide and examples directories respectively.
# -----------------------------------------------------------------

QUIET =
TEX_COMPILER = pdflatex

UG_ROOT = guide/build/latex
UG_FILE = ark_guide

EX_ROOT = examples/build/latex
EX_FILE = ark_examples

#-------------------------------------------------------------------------------

all:
	@echo "Usage:"
	@echo "  make ug    - generate PDF file for User's Guide"
	@echo "  make ex    - generate PDF file for Examples Guide"
	@echo "  make clean - remove all temporary files"
	@echo "  make purge - remove temporary files + PDF files"
	@echo ""

#-------------------------------------------------------------------------------

ug: ug_pdf

ug_pdf: ${UG_ROOT}/${UG_FILE}.tex
	cd ${UG_ROOT}; \
	cp ../../source/figs/* .; \
	cp ../../source/sunmatrix/figs/* .; \
	${TEX_COMPILER} ${QUIET} ${UG_FILE}.tex; \
	makeindex ${UG_FILE}.idx; \
	${TEX_COMPILER} ${QUIET} ${UG_FILE}.tex; \
	${TEX_COMPILER} ${QUIET} ${UG_FILE}.tex; \
	cd -; \
	cp ${UG_ROOT}/${UG_FILE}.pdf .

#-------------------------------------------------------------------------------

ex: ex_pdf

ex_pdf: ${EX_ROOT}/${EX_FILE}.tex
	cd ${EX_ROOT}; \
	cp ../../source/figs/* .; \
	${TEX_COMPILER} ${QUIET} ${EX_FILE}.tex; \
	makeindex ${EX_FILE}.idx; \
	${TEX_COMPILER} ${QUIET} ${EX_FILE}.tex; \
	${TEX_COMPILER} ${QUIET} ${EX_FILE}.tex; \
	cd -; \
	cp ${EX_ROOT}/${EX_FILE}.pdf .

#-------------------------------------------------------------------------------

clean:
	cd ${UG_ROOT}; \
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out
	cd ${EX_ROOT}; \
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out

purge:  clean
	rm -f ${UG_FILE}.pdf
	rm -f ${EX_FILE}.pdf
	cd ${UG_ROOT}; rm -f ${UG_FILE}.pdf
	cd ${EX_ROOT}; rm -f ${EX_FILE}.pdf
