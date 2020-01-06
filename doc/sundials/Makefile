# -----------------------------------------------------------------
# Programmer(s): Carol Woodward, Cody Balos @ LLNL
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
# Makefile for SUNDIALS Installation Guide
# -----------------------------------------------------------------

QUIET =

USER_GUIDE = install_guide
SHARED_DIR = ../sundials

# Use pdflatex as default compiler
TEX_COMPILER = pdflatex

# names of guide files
UG_FILES = \
${USER_GUIDE}.tex

UG_SHARED_FILES = \
install.tex \
ug.tex \
biblio.bib 

UG_SHARED_FIGURES = \
doc_logo_blue.pdf \
warning.pdf

UG_CMAKE_FIGURES = \
ccmakedefault.png \
ccmakeprefix.png \
testrunnertop.png \
testrunnerbot.png
#ccmakeempty.pdf \
#cmaketest.pdf

UG_SHARED = $(addprefix $(SHARED_DIR)/, ${UG_SHARED_FILES})
UG_SHARED += $(addprefix $(SHARED_DIR)/figures/, ${UG_SHARED_FIGURES})
UG_SHARED += $(addprefix $(SHARED_DIR)/figures/cmake/, ${UG_CMAKE_FIGURES})

export TEXINPUTS=./:${SHARED_DIR}//:${SHARED_DIR}/figures//:
export BIBINPUTS=./:${SHARED_DIR}//:
export BSTINPUTS=./:${SHARED_DIR}//:

#-------------------------------------------------------------------------------

all:
	@echo "Usage:"
	@echo "  make in_g     - compile LaTeX files to generate PDF file for Installation Guide"
	@echo "  make clean    - remove all temporary files"
	@echo "  make purge    - remove temporary files + PDF files"
	@echo "  make in_g_pdf - use pdflatex to generate PDF file for Installation Guide"
	@echo "  make in_g_xel - use xelatex to generate PDF file for Installation Guide"
	@echo ""

#-------------------------------------------------------------------------------

in_g: $(USER_GUIDE).pdf

in_g_pdf: use_pdflatex in_g

in_g_xel: use_xellatex in_g

${USER_GUIDE}.pdf: ${UG_FILES} ${UG_SHARED}
	${TEX_COMPILER} ${QUIET} ${USER_GUIDE}.tex
	${TEX_COMPILER} ${QUIET} ${USER_GUIDE}.tex

#-------------------------------------------------------------------------------

use_pdflatex:
	$(eval TEX_COMPILER := pdflatex)

clean:
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out

purge:  clean
	rm -f ${USER_GUIDE}.pdf
