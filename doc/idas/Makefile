# -----------------------------------------------------------------
# $Revision: 1.7 $
# $Date: 2012-04-26 00:06:12 $
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

# make sure doc/texi2dvi in PATH
export PATH := ../:$(PATH)

SHARED_DOC = ../sundials

HTML_DIR = ../../html/documentation

UG_FILES = \
idas_guide.tex \
idas_intro.tex \
idas_install.tex \
idas_math.tex \
idas_sim.tex \
idas_fwd.tex \
idas_adj.tex \
idas_nvec.tex \
idas_org.tex \
idas_newlinsolv.tex \
idas_linsolv.tex \
idas_constants.tex

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
idas_examples.tex \
idas_ex_intro.tex \
idas_ex_fwd.tex \
idas_ex_adj.tex \
idas_ex_codes.tex

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
	@echo "  make ug     - compile LaTeX files to generate DVI file for IDAS User's Guide"
	@echo "  make ug_pdf - compile LaTeX files to generate PS + PDF files for IDAS User's Guide"
	@echo "  make ex     - compile LaTeX files to generate DVI file for IDAS Examples Guide"
	@echo "  make ex_pdf - compile LaTeX files to generate PS + PDF files for IDAS Examples Guide"
	@echo "  make clean  - remove all temporary files"
	@echo "  make purge  - remove temporary files + PS + PDF files"
	@echo ""

install: ug_pdf ex_pdf
#install: ug_pdf ex_pdf ug_html ex_html
#	rm -rf $(HTML_DIR)/idas_guide
#	rm -rf $(HTML_DIR)/idas_examples
#	cp idas_guide.ps $(HTML_DIR)
	cp idas_guide.pdf $(HTML_DIR)
#	cp -r idas_guide $(HTML_DIR)
#	cp idas_examples.ps $(HTML_DIR)
	cp idas_examples.pdf $(HTML_DIR)
#	cp -r idas_examples $(HTML_DIR)

#-------------------------------------------------------------------------------

ug: 	copy_ug idas_guide.aux

ug_pdf: ug idas_guide.ps

#ug_html: ug
#	cp ug_html.tex ug.tex
#	latex2html -split 3 -link 3 -info 1 -no_math -html_version 3.2,math -show_section_numbers -short_index -local_icons -white -notransparent idas_guide.tex
#	@for i in ${HTML_BUTTONS} ; do \
#	cp ${SHARED_DOC}/buttons/$${i} idas_guide/ ;\
#	done
#	@cp ${SHARED_DOC}/l2h_styles.css idas_guide/idas_guide.css
#	./processHTML idas_guide

copy_ug:
	@for i in ${UG_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

idas_guide.aux: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} idas_guide.tex



idas_guide.ps: idas_guide.dvi
	dvips ${QUIET} -Ppdf idas_guide
	ps2pdf idas_guide.ps

#-------------------------------------------------------------------------------

ex:	ug copy_ex idas_examples.aux

ex_pdf: ex idas_examples.ps

#ex_html: ex
#	cp ug_html.tex ug.tex
#	latex2html -split 4 -link 4 -info 1 -no_math -html_version 3.2,math -show_section_numbers -local_icons -white -notransparent idas_examples.tex
#	@for i in ${HTML_BUTTONS} ; do \
#	cp ${SHARED_DOC}/buttons/$${i} idas_examples/ ;\
#	done
#	@cp ${SHARED_DOC}/l2h_styles.css idas_examples/idas_examples.css
#	./processHTML idas_examples

copy_ex:
	@for i in ${EX_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

idas_examples.aux: ${EX_FILES} ${EX_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} idas_examples.tex

idas_examples.ps: idas_examples.dvi
	dvips ${QUIET} -Ppdf idas_examples
	ps2pdf idas_examples.ps

#-------------------------------------------------------------------------------

clean:
	rm -rf ${UG_SHARED_FILES}
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out
purge:  clean
	rm -f idas_guide.ps idas_guide.pdf
	rm -f idas_examples.ps idas_examples.pdf
	rm -rf idas_guide idas_examples
