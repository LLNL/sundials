# -----------------------------------------------------------------
# $Revision: 1.4 $
# $Date: 2012-04-26 00:06:10 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see sundials/ida/LICENSE.
# -----------------------------------------------------------------
# Makefile for IDA documents
# -----------------------------------------------------------------

QUIET =

SHARED_DOC = ../sundials

HTML_DIR = ../../html/documentation

UG_FILES = \
ida_guide.tex \
ida_intro.tex \
ida_install.tex \
ida_math.tex \
ida_usage.tex \
ida_fcmix.tex \
ida_nvec.tex \
ida_org.tex \
ida_newlinsolv.tex \
#ida_linsolv.tex \
ida_constants.tex

UG_SHARED_FILES = \
install.tex \
linsolv.tex \
nvec_generic.tex \
nvec_parallel.tex \
nvec_serial.tex \
nvec_openmp.tex \
nvec_pthreads.tex \
spils_constants.tex \
sun_org.tex \
types.tex \
ug.tex \
biblio.bib \
texi2dvi
#ug_html.tex \
#processHTML

UG_SHARED_FIGURES = \
bandmat.eps \
cscmat.eps \
sunorg1.eps \
sunorg2.eps \
doc_logo.eps

UG_CMAKE_FIGURES = \
ccmakedefault.eps \
ccmakeempty.eps \
ccmakeprefix.eps \
testrunnertop.eps\
testrunnerbot.eps\
cmaketest.eps

EX_FILES = \
ida_examples.tex \
ida_ex_intro.tex \
ida_ex_serial.tex \
ida_ex_parallel.tex \
ida_ex_codes.tex

EX_SHARED_FILES = \
ug.tex \
biblio.bib \
texi2dvi
#ug_html.tex \
#processHTML

#HTML_BUTTONS =\
#back.png \
#contents.png \
#index.png \
#next.png \
#next_g.png \
#prev.png \
#prev_g.png \
#up.png \
#up_g.png \
#warning.png

#-------------------------------------------------------------------------------

all:
	@echo "Usage:"
	@echo "  make ug     - compile LaTeX files to generate DVI file for IDA User's Guide"
	@echo "  make ug_pdf - compile LaTeX files to generate PS + PDF files for IDA User's Guide"
	@echo "  make ex     - compile LaTeX files to generate DVI file for IDA Examples Guide"
	@echo "  make ex_pdf - compile LaTeX files to generate PS + PDF files for IDA Examples Guide"
	@echo "  make clean  - remove all temporary files"
	@echo "  make purge  - remove temporary files + PS + PDF files"
	@echo ""

install: ug_pdf ex_pdf
#install: ug_pdf ex_pdf ug_html ex_html
#	rm -rf $(HTML_DIR)/ida_guide
#	rm -rf $(HTML_DIR)/ida_examples
#	cp ida_guide.ps $(HTML_DIR)
	cp ida_guide.pdf $(HTML_DIR)
#	cp -r ida_guide $(HTML_DIR)
#	cp ida_examples.ps $(HTML_DIR)
	cp ida_examples.pdf $(HTML_DIR)
#	cp -r ida_examples $(HTML_DIR)

#-------------------------------------------------------------------------------

ug: 	copy_ug copy_ug_figs copy_ug_cmake_figs ida_guide.aux

ug_pdf: ug ida_guide.ps

#ug_html: ug
#	cp ug_html.tex ug.tex
#	latex2html -split 3 -link 3 -info 1 -no_math -html_version 3.2,math -show_section_numbers -short_index -local_icons -white -notransparent ida_guide.tex
#	@for i in ${HTML_BUTTONS} ; do \
#	cp ${SHARED_DOC}/buttons/$${i} ida_guide/ ;\
#	done
#	@cp ${SHARED_DOC}/l2h_styles.css ida_guide/ida_guide.css
#	./processHTML ida_guide

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

ida_guide.aux: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} ida_guide.tex

ida_guide.ps: ida_guide.dvi
	dvips ${QUIET} -Ppdf ida_guide
	ps2pdf ida_guide.ps

#-------------------------------------------------------------------------------

ex:	ug copy_ex ida_examples.aux

ex_pdf: ex ida_examples.ps

#ex_html: ex
#	cp ug_html.tex ug.tex
#	latex2html -split 4 -link 4 -info 1 -no_math -html_version 3.2,math -show_section_numbers -local_icons -white -notransparent ida_examples.tex
#	@for i in ${HTML_BUTTONS} ; do \
#	cp ${SHARED_DOC}/buttons/$${i} ida_examples/ ;\
#	done
#	@cp ${SHARED_DOC}/l2h_styles.css ida_examples/ida_examples.css
#	./processHTML ida_examples

copy_ex:
	@for i in ${EX_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

ida_examples.aux: ${EX_FILES} ${EX_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} ida_examples.tex

ida_examples.ps: ida_examples.dvi
	dvips ${QUIET} -Ppdf ida_examples
	ps2pdf ida_examples.ps

#-------------------------------------------------------------------------------

clean:
	rm -rf ${UG_SHARED_FILES}
	rm -rf ${UG_SHARED_FIGURES}
	rm -rf ${UG_CMAKE_FIGURES}
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out
purge:  clean
	rm -f ida_guide.ps ida_guide.pdf
	rm -f ida_examples.ps ida_examples.pdf
	rm -rf ida_guide ida_examples
