# Makefile for Sphinx LaTeX output
ARK_UG = ARKode.pdf
ARK_EX = ARKode_example.pdf
SHARED_DOC = ../sundials
HTML_DIR = ../../html/documentation
# make sure doc/texi2dvi in PATH

export PATH := ../:$(PATH)

ALLDOCS = $(basename $(wildcard *.tex))
ALLPDF = $(addsuffix .pdf,$(ALLDOCS))
ALLDVI = $(addsuffix .dvi,$(ALLDOCS))

# Prefix for archive names
ARCHIVEPRREFIX =
# Additional LaTeX options
LATEXOPTS =

usage:
	@echo "Usage:"
	@echo "  make ug     - compile LaTeX files to generate PDF file for ARKode User's Guide"
	@echo "  make ex     - compile LaTeX files to generate PDF file for ARKode Examples Guide"
	@echo "  make clean  - remove all temporary files"
	@echo "  make purge  - remove temporary files + PDF files"

ug: ug_pdf

ug_pdf: $(ARK_UG)

ex: ex_pdf

ex_pdf: $(ARK_EX)

purge: clean
	rm -f $(ARK_UG)
	rm -f $(ARK_EX)

install: ug_pdf ex_pdf
	cp $(ARK_UG) $(HTML_DIR)
	cp $(ARK_EX) $(HTML_DIR)



all: $(ALLPDF)
all-pdf: $(ALLPDF)
all-dvi: $(ALLDVI)
all-ps: all-dvi
	for f in *.dvi; do dvips $$f; done

all-pdf-ja:
	for f in *.pdf *.png *.gif *.jpg *.jpeg; do extractbb $$f; done
	for f in *.tex; do platex -kanji=utf8 $(LATEXOPTS) $$f; done
	for f in *.tex; do platex -kanji=utf8 $(LATEXOPTS) $$f; done
	for f in *.tex; do platex -kanji=utf8 $(LATEXOPTS) $$f; done
	-for f in *.idx; do mendex -U -f -d "`basename $$f .idx`.dic" -s python.ist $$f; done
	for f in *.tex; do platex -kanji=utf8 $(LATEXOPTS) $$f; done
	for f in *.tex; do platex -kanji=utf8 $(LATEXOPTS) $$f; done
	for f in *.dvi; do dvipdfmx $$f; done

zip: all-$(FMT)
	mkdir $(ARCHIVEPREFIX)docs-$(FMT)
	cp $(ALLPDF) $(ARCHIVEPREFIX)docs-$(FMT)
	zip -q -r -9 $(ARCHIVEPREFIX)docs-$(FMT).zip $(ARCHIVEPREFIX)docs-$(FMT)
	rm -r $(ARCHIVEPREFIX)docs-$(FMT)

tar: all-$(FMT)
	mkdir $(ARCHIVEPREFIX)docs-$(FMT)
	cp $(ALLPDF) $(ARCHIVEPREFIX)docs-$(FMT)
	tar cf $(ARCHIVEPREFIX)docs-$(FMT).tar $(ARCHIVEPREFIX)docs-$(FMT)
	rm -r $(ARCHIVEPREFIX)docs-$(FMT)

bz2: tar
	bzip2 -9 -k $(ARCHIVEPREFIX)docs-$(FMT).tar

# The number of LaTeX runs is quite conservative, but I don't expect it
# to get run often, so the little extra time won't hurt.
%.dvi: %.tex
	latex $(LATEXOPTS) '$<'
	latex $(LATEXOPTS) '$<'
	latex $(LATEXOPTS) '$<'
	-makeindex -s python.ist '$(basename $<).idx'
	latex $(LATEXOPTS) '$<'
	latex $(LATEXOPTS) '$<'

%.pdf: %.tex
	pdflatex $(LATEXOPTS) '$<'
	pdflatex $(LATEXOPTS) '$<'
	pdflatex $(LATEXOPTS) '$<'
	-makeindex -s python.ist '$(basename $<).idx'
	pdflatex $(LATEXOPTS) '$<'
	pdflatex $(LATEXOPTS) '$<'

clean:
	rm -f *.dvi *.log *.ind *.aux *.toc *.syn *.idx *.out *.ilg *.pla *~

.PHONY: all all-pdf all-dvi all-ps clean
.PHONY: all-pdf-ja

