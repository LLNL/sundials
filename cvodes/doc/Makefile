SHARED_DOC = ../../shared/doc
all:
	@echo 'Usage: make doc   - compile LaTeX files to generate cvs_guide.ps and cvs_guide.pdf'
	@echo '       make clean - remove temporary files (*.aux, *.log, etc.)'
	@echo
doc:
	cp $(SHARED_DOC)/band.tex .
	cp $(SHARED_DOC)/dense.tex .
	cp $(SHARED_DOC)/nvec_generic.tex .
	cp $(SHARED_DOC)/nvec_parallel.tex .
	cp $(SHARED_DOC)/nvec_serial.tex .
	cp $(SHARED_DOC)/spgmr.tex .
	cp $(SHARED_DOC)/sun_org.tex .
	cp $(SHARED_DOC)/types.tex .
	cp $(SHARED_DOC)/bandmat.eps .
	cp $(SHARED_DOC)/sunorg.eps .
	latex cvs_guide
	latex cvs_guide
	makeindex -s cvs_indx.ist cvs_guide
	latex cvs_guide
	latex cvs_guide
	dvips -Ppdf -o cvs_guide.ps cvs_guide
	ps2pdf cvs_guide.ps cvs_guide.pdf
clean:
	rm -f bandmat.eps
	rm -f band.tex
	rm -f dense.tex
	rm -f nvec_generic.tex
	rm -f nvec_parallel.tex
	rm -f nvec_serial.tex
	rm -f spgmr.tex
	rm -f sunorg.eps
	rm -f sun_org.tex
	rm -f types.tex
	rm -f *~ *.aux *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
#	rm -f cvs_guide.ps cvs_guide.pdf

