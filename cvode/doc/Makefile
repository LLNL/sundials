# -----------------------------------------------------------------
# $Revision: 1.8 $
# $Date: 2004-10-21 19:04:36 $
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

all:
	@echo "Usage:"
	@echo "  make ug - compile LaTeX files to generate CVODE User's Guide"
	@echo "  make ex - compile LaTeX files to generate CVODE Examples Guide"
	@echo "  make clean - remove all files"
	@echo ""

ug:
	cp $(SHARED_DOC)/band.tex .
	cp $(SHARED_DOC)/install.tex .
	cp $(SHARED_DOC)/dense.tex .
	cp $(SHARED_DOC)/nvec_generic.tex .
	cp $(SHARED_DOC)/nvec_parallel.tex .
	cp $(SHARED_DOC)/nvec_serial.tex .
	cp $(SHARED_DOC)/spgmr.tex .
	cp $(SHARED_DOC)/sun_org.tex .
	cp $(SHARED_DOC)/types.tex .
	cp $(SHARED_DOC)/ug.tex .
	cp $(SHARED_DOC)/biblio.bib .
	cp $(SHARED_DOC)/bandmat.eps .
	cp $(SHARED_DOC)/sunorg1.eps .
	cp $(SHARED_DOC)/sunorg2.eps .
	cp $(SHARED_DOC)/llnlCoverPage.tex .
	cp $(SHARED_DOC)/Logo_for_papers.ps .
	cp $(SHARED_DOC)/Rule_and_address.ps .
	cp $(SHARED_DOC)/chngpage.sty .
	latex cv_guide
	latex cv_guide
	bibtex cv_guide
	makeindex -s cv_indx.ist cv_guide
	latex cv_guide
	latex cv_guide
#	dvips -Ppdf -o cv_guide.ps cv_guide
#	ps2pdf cv_guide.ps cv_guide.pdf

ex:
	cp $(SHARED_DOC)/llnlCoverPage.tex .
	cp $(SHARED_DOC)/Logo_for_papers.ps .
	cp $(SHARED_DOC)/Rule_and_address.ps .
	cp $(SHARED_DOC)/chngpage.sty .
	cp $(SHARED_DOC)/ug.tex .
	cp $(SHARED_DOC)/biblio.bib .
	latex cv_examples
	bibtex cv_examples
	latex cv_examples
	latex cv_examples
#	dvips -Ppdf -o cv_examples.ps cv_examples
#	ps2pdf cv_examples.ps cv_examples.pdf

clean:
	rm -f bandmat.eps
	rm -f band.tex
	rm -f biblio.bib
	rm -f dense.tex
	rm -f install.tex
	rm -f nvec_generic.tex
	rm -f nvec_parallel.tex
	rm -f nvec_serial.tex
	rm -f spgmr.tex
	rm -f sunorg1.eps
	rm -f sunorg2.eps
	rm -f sun_org.tex
	rm -f types.tex
	rm -f ug.tex
	rm -f llnlCoverPage.tex Logo_for_papers.ps Rule_and_address.ps chngpage.sty
	rm -f *~ *.aux *.bbl *.blg *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
	rm -f cv_guide.ps cv_guide.pdf
	rm -f cv_examples.ps cv_examples.pdf
