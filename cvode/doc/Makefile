# -----------------------------------------------------------------
# $Revision: 1.9 $
# $Date: 2004-11-08 17:48:30 $
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

SHARED_DOC = ../../shared/doc

UG_FILES = \
cv_guide.tex \
cv_intro.tex \
cv_install.tex \
cv_math.tex \
cv_usage.tex \
cv_nvec.tex \
cv_org.tex \
cv_newlinsolv.tex \
cv_linsolv.tex

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
biblio.bib \
bandmat.eps \
sunorg1.eps \
sunorg2.eps \
llnlCoverPage.tex \
Logo_for_papers.ps \
Rule_and_address.ps \
chngpage.sty

EX_FILES = \
cv_examples.tex \
cv_ex_intro.tex \
cv_ex_serial.tex \
cv_ex_parallel.tex \
cv_ex_fortran.tex \
cv_ex_codes.tex \
cv_ex_tests.tex   

EX_SHARED_FILES = \
llnlCoverPage.tex \
Logo_for_papers.ps \
Rule_and_address.ps \
chngpage.sty \
ug.tex \
biblio.bib

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

#-------------------------------------------------------------------------------

ug: 	copy_ug cv_guide.aux

ug_pdf: ug 
	dvips -Ppdf cv_guide
	ps2pdf cv_guide.ps

copy_ug:
	@for i in ${UG_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

cv_guide.aux: ${UG_FILES} ${UG_SHARED_FILES:%=${SHARED_DOC}/%}
	latex cv_guide
	latex cv_guide
	bibtex cv_guide
	makeindex -s cv_indx.ist cv_guide
	latex cv_guide
	latex cv_guide

#-------------------------------------------------------------------------------

ex:	copy_ex cv_examples.aux

ex_pdf: ex 
	dvips -Ppdf cv_examples
	ps2pdf cv_examples.ps

copy_ex:
	@for i in ${EX_SHARED_FILES} ; do \
	cp ${SHARED_DOC}/$${i} . ;\
	done

cv_examples.aux: ${EX_FILES} ${EX_SHARED_FILES:%=${SHARED_DOC}/%}
	latex cv_examples
	bibtex cv_examples
	latex cv_examples
	latex cv_examples

#-------------------------------------------------------------------------------

clean:
	rm -rf ${UG_SHARED_FILES}
	rm -rf ${EX_SHARED_FILES}
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
purge:  clean
	rm -f cv_guide.ps cv_guide.pdf
	rm -f cv_examples.ps cv_examples.pdf
