# -----------------------------------------------------------------
# $Revision: 1.6 $
# $Date: 2004-11-11 02:06:52 $
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

SHARED_DOC = ../../shared/doc

UG_FILES = \
ida_guide.tex \
ida_intro.tex \
ida_install.tex \
ida_math.tex \
ida_usage.tex \
ida_nvec.tex \
ida_org.tex \
ida_newlinsolv.tex \
ida_linsolv.tex

UG_SHARED_FILES = \
band.tex \
install.tex \
dense.tex \
nvec_generic.tex \
nvec_parallel.tex \
nvec_serial.tex \
spgmr.tex \
sun_org.tex \
types.tex \
ug.tex \
ug_html.tex \
biblio.bib \
bandmat.eps \
sunorg1.eps \
sunorg2.eps \
llnlCoverPage.tex \
Logo_for_papers.ps \
Rule_and_address.ps \
chngpage.sty

EX_FILES = \
ida_examples.tex \
ida_ex_intro.tex \
ida_ex_serial.tex \
ida_ex_parallel.tex \
ida_ex_codes.tex

EX_SHARED_FILES = \
llnlCoverPage.tex \
Logo_for_papers.ps \
Rule_and_address.ps \
chngpage.sty \
ug.tex \
ug_html.tex \
biblio.bib

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

#-------------------------------------------------------------------------------

ug: 	copy_ug ida_guide.aux

ug_pdf: ug ida_guide.ps

ug_html: ug
	cp ug_html.tex ug.tex
	latex2html -split 3 -link 3 -info 1 -no_math -html_version 3.2,math -short_index -local_icons -white ida_guide.tex

copy_ug:
	@for i in ${UG_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

ida_guide.aux: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} ida_guide.tex

ida_guide.ps: ida_guide.dvi
	dvips ${QUIET} -Ppdf ida_guide
	ps2pdf ida_guide.ps

#-------------------------------------------------------------------------------

ex:	copy_ex ida_examples.aux

ex_pdf: ex ida_examples.ps

ex_html: ex
	cp ug_html.tex ug.tex
	latex2html -split 4 -link 4 -info 1 -no_math -html_version 3.2,math -local_icons -white ida_examples.tex

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
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
purge:  clean
	rm -f ida_guide.ps ida_guide.pdf
	rm -f ida_examples.ps ida_examples.pdf

