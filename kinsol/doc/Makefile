SHARED_DOC = ../../shared/doc
all:
	@echo 'Usage: make doc   - compile LaTeX files to generate kinudoc.ps and kinudoc.pdf'
	@echo '       make clean - remove temporary files (*.aux, *.log, etc.)'
	@echo
doc:
	latex kinudoc
	latex kinudoc
	dvips -Ppdf -o kinudoc.ps kinudoc
	ps2pdf kinudoc.ps kinudoc.pdf
clean:
	rm -f *~ *.aux *.log *.toc *.lot *.lof *.ind *.ilg *.idx *.dvi
#	rm -f kinudoc.ps kinudoc.pdf

