# -----------------------------------------------------------------
# $Revision: 1.12 $
# $Date: 2005-07-12 23:56:10 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see sundials/kinsol/LICENSE.
# -----------------------------------------------------------------
# Makefile for KINSOL documents
# -----------------------------------------------------------------

QUIET =

SHARED_DOC = ../../shared/doc

UG_FILES = \
kin_guide.tex \
kin_intro.tex \
kin_install.tex \
kin_math.tex \
kin_use.tex \
kin_nvec.tex \
kin_org.tex \
kin_newlinsolv.tex \
kin_linsolv.tex \
kin_constants.tex

UG_SHARED_FILES = \
band.tex \
install.tex \
dense.tex \
nvec_generic.tex \
nvec_parallel.tex \
nvec_serial.tex \
spgmr.tex \
spbcg.tex \
sptfqmr.tex \
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
kin_examples.tex \
kin_ex_intro.tex \
kin_ex_c.tex \
kin_ex_fortran.tex \
kin_ex_codes.tex

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
	@echo "  make ug     - compile LaTeX files to generate DVI file for KINSOL User's Guide"
	@echo "  make ug_pdf - compile LaTeX files to generate PS + PDF files for KINSOL User's Guide"
	@echo "  make ex     - compile LaTeX files to generate DVI file for KINSOL Examples Guide"
	@echo "  make ex_pdf - compile LaTeX files to generate PS + PDF files for KINSOL Examples Guide"
	@echo "  make clean  - remove all temporary files"
	@echo "  make purge  - remove temporary files + PS + PDF files"
	@echo ""

#-------------------------------------------------------------------------------

ug: 	copy_ug kin_guide.aux

ug_pdf: ug kin_guide.ps

ug_html: ug
	cp ug_html.tex ug.tex
	latex2html -split 3 -link 3 -info 1 -no_math -html_version 3.2,math -short_index -show_section_numbers -local_icons -white -notransparent kin_guide.tex

copy_ug:
	@for i in ${UG_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

kin_guide.aux: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} kin_guide.tex

kin_guide.ps: kin_guide.dvi
	dvips ${QUIET} -Ppdf kin_guide
	ps2pdf kin_guide.ps

#-------------------------------------------------------------------------------

ex:	ug copy_ex kin_examples.aux

ex_pdf: ex kin_examples.ps

ex_html: ex
	cp ug_html.tex ug.tex
	latex2html -split 4 -link 4 -info 1 -no_math -html_version 3.2,math -show_section_numbers -local_icons -white -notransparent kin_examples.tex

copy_ex:
	@for i in ${EX_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

kin_examples.aux: ${EX_FILES} ${EX_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} kin_examples.tex

kin_examples.ps: kin_examples.dvi
	dvips ${QUIET} -Ppdf kin_examples
	ps2pdf kin_examples.ps

#-------------------------------------------------------------------------------

clean:
	rm -rf ${UG_SHARED_FILES}
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
purge:  clean
	rm -f kin_guide.ps kin_guide.pdf
	rm -f kin_examples.ps kin_examples.pdf
