SHARED_DOC = ../../shared/doc
all:
	@echo 'Usage: make doc   - compile LaTeX files to generate idaudoc.ps and idaudoc.pdf'
	@echo '       make clean - remove temporary files (*.aux, *.log, etc.)'
	@echo
doc:
	latex idaudoc
	latex idaudoc
	dvips -Ppdf -o idaudoc.ps idaudoc
	ps2pdf idaudoc.ps idaudoc.pdf
clean:
	rm -f *~ *.aux *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
#	rm -f idaudoc.ps idaudoc.pdf

