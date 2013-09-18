# -----------------------------------------------------------------
# $Revision: 1.4 $
# $Date: 2012-04-26 00:05:57 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see sundials/cvode/LICENSE.
# -----------------------------------------------------------------
# Makefile for CVODE documents
# -----------------------------------------------------------------

QUIET =

SHARED_DOC = ../sundials

HTML_DIR = ../../html/documentation

UG_FILES = \
cv_guide.tex \
cv_intro.tex \
cv_install.tex \
cv_math.tex \
cv_usage.tex \
cv_fcmix.tex \
cv_nvec.tex \
cv_org.tex \
cv_newlinsolv.tex \
cv_linsolv.tex \
cv_constants.tex

UG_SHARED_FILES = \
install.tex \
linsolv.tex \
nvec_generic.tex \
nvec_parallel.tex \
nvec_serial.tex \
nvec_openmp.tex \
nvec_pthreads.tex \
sun_org.tex \
types.tex \
ug.tex \
biblio.bib \
bandmat.eps \
sunorg1.eps \
sunorg2.eps \
doc_logo.eps \
#ug_html.tex \
#processHTML

EX_FILES = \
cv_examples.tex \
cv_ex_intro.tex \
cv_ex_serial.tex \
cv_ex_parallel.tex \
cv_ex_fortran.tex \
cv_ex_codes.tex \
cv_ex_tests.tex   

EX_SHARED_FILES = \
ug.tex \
biblio.bib \
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
	@echo "  make ug     - compile LaTeX files to generate DVI file for CVODE User's Guide"
	@echo "  make ug_pdf - compile LaTeX files to generate PS + PDF files for CVODE User's Guide"
	@echo "  make ex     - compile LaTeX files to generate DVI file for CVODE Examples Guide"
	@echo "  make ex_pdf - compile LaTeX files to generate PS + PDF files for CVODE Examples Guide"
	@echo "  make clean  - remove all temporary files"
	@echo "  make purge  - remove temporary files + PS + PDF files"
	@echo ""

install: ug_pdf ex_pdf
#install: ug_pdf ex_pdf ug_html ex_html
#	rm -rf $(HTML_DIR)/cv_guide
#	rm -rf $(HTML_DIR)/cv_examples
#	cp cv_guide.ps $(HTML_DIR)
	cp cv_guide.pdf $(HTML_DIR)
#	cp -r cv_guide $(HTML_DIR)
#	cp cv_examples.ps $(HTML_DIR)
	cp cv_examples.pdf $(HTML_DIR)
#	cp -r cv_examples $(HTML_DIR)

#-------------------------------------------------------------------------------

ug: 	copy_ug cv_guide.aux

ug_pdf: ug cv_guide.ps

#ug_html: ug
#	cp ug_html.tex ug.tex
#	latex2html -split 3 -link 3 -info 1 -no_math -html_version 3.2,math -short_index -show_section_numbers -local_icons -white -notransparent cv_guide.tex
#	@for i in ${HTML_BUTTONS} ; do \
#	cp ${SHARED_DOC}/buttons/$${i} cv_guide/ ;\
#	done
#	@cp ${SHARED_DOC}/l2h_styles.css cv_guide/cv_guide.css
#	./processHTML cv_guide

copy_ug:
	@for i in ${UG_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

cv_guide.aux: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} cv_guide.tex

cv_guide.ps: cv_guide.dvi
	dvips ${QUIET} -Ppdf cv_guide
	ps2pdf cv_guide.ps

#-------------------------------------------------------------------------------

ex:	ug copy_ex cv_examples.aux

ex_pdf: ex cv_examples.ps

#ex_html: ex
#	cp ug_html.tex ug.tex
#	latex2html -split 4 -link 4 -info 1 -no_math -html_version 3.2,math -show_section_numbers -local_icons -white -notransparent cv_examples.tex
#	@for i in ${HTML_BUTTONS} ; do \
#	cp ${SHARED_DOC}/buttons/$${i} cv_examples/ ;\
#	done
#	@cp ${SHARED_DOC}/l2h_styles.css cv_examples/cv_examples.css
#	./processHTML cv_examples

copy_ex:
	@for i in ${EX_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

cv_examples.aux: ${EX_FILES} ${EX_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} cv_examples.tex

cv_examples.ps: cv_examples.dvi
	dvips ${QUIET} -Ppdf cv_examples
	ps2pdf cv_examples.ps

#-------------------------------------------------------------------------------

clean:
	rm -rf ${UG_SHARED_FILES}
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out

purge:  clean
	rm -f cv_guide.ps cv_guide.pdf
	rm -f cv_examples.ps cv_examples.pdf
	rm -rf cv_guide cv_examples
