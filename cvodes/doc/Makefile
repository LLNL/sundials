# -----------------------------------------------------------------
# $Revision: 1.18 $
# $Date: 2005-11-29 23:15:30 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see sundials/cvodes/LICENSE.
# -----------------------------------------------------------------
# Makefile for CVODES documents
# -----------------------------------------------------------------

QUIET =

SHARED_DOC = ../../shared/doc

UG_FILES = \
cvs_guide.tex \
cvs_intro.tex \
cvs_install.tex \
cvs_math.tex \
cvs_org.tex \
cvs_sim.tex \
cvs_fwd.tex \
cvs_adj.tex \
cvs_nvec.tex \
cvs_newlinsolv.tex \
cvs_linsolv.tex \
cvs_constants.tex

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
cvs_examples.tex \
cvs_ex_intro.tex \
cvs_fwd_ex.tex \
cvs_adj_ex.tex \
cvs_ex_codes.tex \
cvs_ex_tests.tex   

EX_SHARED_FILES = \
llnlCoverPage.tex \
Logo_for_papers.ps \
Rule_and_address.ps \
chngpage.sty \
ug.tex \
ug_html.tex \
biblio.bib

all:
	@echo "Usage:"
	@echo "  make ug     - compile LaTeX files to generate DVI file for CVODES User's Guide"
	@echo "  make ug_pdf - compile LaTeX files to generate PS + PDF files for CVODES User's Guide"
	@echo "  make ex     - compile LaTeX files to generate DVI file for CVODES Examples Guide"
	@echo "  make ex_pdf - compile LaTeX files to generate PS + PDF files for CVODES Examples Guide"
	@echo "  make clean  - remove all temporary files"
	@echo "  make purge  - remove temporary files + PS + PDF files"
	@echo ""

#-------------------------------------------------------------------------------

ug: 	copy_ug cvs_guide.aux

ug_pdf: ug cvs_guide.ps

ug_html: ug
	cp ug_html.tex ug.tex
	latex2html -split 3 -link 3 -info 1 -no_math -html_version 3.2,math -short_index -show_section_numbers -local_icons -white -notransparent cvs_guide.tex

copy_ug:
	@for i in ${UG_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

cvs_guide.aux: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} cvs_guide.tex

cvs_guide.ps: cvs_guide.dvi
	dvips ${QUIET} -Ppdf cvs_guide
	ps2pdf cvs_guide.ps

#-------------------------------------------------------------------------------

ex:	ug copy_ex cvs_examples.aux

ex_pdf: ex cvs_examples.ps

ex_html: ex
	cp ug_html.tex ug.tex
	latex2html -split 4 -link 4 -info 1 -no_math -html_version 3.2,math -show_section_numbers -local_icons -white -notransparent cvs_examples.tex

copy_ex:
	@for i in ${EX_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

cvs_examples.aux: ${EX_FILES} ${EX_SHARED_FILES:%=${SHARED_DOC}/%}
	texi2dvi ${QUIET} cvs_examples.tex

cvs_examples.ps: cvs_examples.dvi
	dvips ${QUIET} -Ppdf cvs_examples
	ps2pdf cvs_examples.ps

#-------------------------------------------------------------------------------

clean:
	rm -rf ${UG_SHARED_FILES}
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi *.out
purge:  clean
	rm -f cvs_guide.ps cvs_guide.pdf
	rm -f cvs_examples.ps cvs_examples.pdf
