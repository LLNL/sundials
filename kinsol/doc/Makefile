SHARED_DOC = ../../shared/doc
all:
	@echo 'Usage: make ug   - compile LaTeX files to generate kinsol user guide'
	@echo 'Usage: make ex   - compile LaTeX files to generate kinsol examples doc'
	@echo '       make clean - remove temporary files (*.aux, *.log, etc.)'
	@echo

ug:
	cp $(SHARED_DOC)/install.tex .
	cp $(SHARED_DOC)/nvec_generic.tex .
	cp $(SHARED_DOC)/nvec_parallel.tex .
	cp $(SHARED_DOC)/nvec_serial.tex .
	cp $(SHARED_DOC)/spgmr.tex .
	cp $(SHARED_DOC)/sun_org.tex .
	cp $(SHARED_DOC)/types.tex .
	cp $(SHARED_DOC)/ug.tex .
	cp $(SHARED_DOC)/biblio.bib .
	cp $(SHARED_DOC)/sunorg1.eps .
	cp $(SHARED_DOC)/sunorg2.eps .
	cp $(SHARED_DOC)/llnlCoverPage.tex .
	cp $(SHARED_DOC)/Logo_for_papers.ps .
	cp $(SHARED_DOC)/Rule_and_address.ps .
	cp $(SHARED_DOC)/chngpage.sty .
	latex kin_guide
	bibtex kin_guide
	latex kin_guide
	makeindex -s kin_indx.ist kin_guide
	latex kin_guide
	latex kin_guide
#	dvips -Ppdf -o kin_guide.ps kin_guide
#	ps2pdf kin_guide.ps kin_guide.pdf
ex:
	cp $(SHARED_DOC)/llnlCoverPage.tex .
	cp $(SHARED_DOC)/Logo_for_papers.ps .
	cp $(SHARED_DOC)/Rule_and_address.ps .
	cp $(SHARED_DOC)/chngpage.sty .
	cp $(SHARED_DOC)/ug.tex .
	cp $(SHARED_DOC)/biblio.bib .
	latex kin_examples
	bibtex kin_examples
	latex kin_examples
	latex kin_examples
#	dvips -Ppdf -o kin_examples.ps kin_examples
#	ps2pdf kin_examples.ps kin_examples.pdf
clean:
	rm -f biblio.bib
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
	rm -f kin_guide.ps kin_guide.pdf
	rm -f kin_examples.ps kin_examples.pdf

olddoc:
	latex kinudoc
	latex kinudoc
	dvips -Ppdf -o kinudoc.ps kinudoc
	ps2pdf kinudoc.ps kinudoc.pdf

