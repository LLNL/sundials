SHARED_DOC = ../../shared/doc
all:
	@echo 'Usage: make ug   - compile LaTeX files to generate ida user guide'
	@echo 'Usage: make ex   - compile LaTeX files to generate ida examples doc'
	@echo '       make clean - remove temporary files (*.aux, *.log, etc.)'
	@echo
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
	latex ida_guide
	latex ida_guide
	bibtex ida_guide
	makeindex -s ida_indx.ist ida_guide
	latex ida_guide
	latex ida_guide
#	dvips -Ppdf -o ida_guide.ps ida_guide
#	ps2pdf ida_guide.ps ida_guide.pdf
ex:
	cp $(SHARED_DOC)/llnlCoverPage.tex .
	cp $(SHARED_DOC)/Logo_for_papers.ps .
	cp $(SHARED_DOC)/Rule_and_address.ps .
	cp $(SHARED_DOC)/chngpage.sty .
	cp $(SHARED_DOC)/ug.tex .
	cp $(SHARED_DOC)/biblio.bib .
	latex ida_examples
	bibtex ida_examples
	latex ida_examples
	latex ida_examples
#	dvips -Ppdf -o ida_examples.ps ida_examples
#	ps2pdf ida_examples.ps ida_examples.pdf
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
	rm -f ida_guide.ps ida_guide.pdf
	rm -f ida_examples.ps ida_examples.pdf

