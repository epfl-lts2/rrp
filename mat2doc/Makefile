PDFNAME = experiments.pdf

all:	$(PDFNAME)

%.pdf: %.tex Makefile $(shell ls *.bib)
	pdflatex $*
	-bibtex $*
	pdflatex $*
	pdflatex $*

%.pdf: %.lyx Makefile $(shell ls *.bib)
	lyx -f --export pdf2 $<

%.tex: %.lyx Makefile $(shell ls *.bib)
	lyx -f --export latex $<

texclean:	
	-rm *.aux *.bbl *.blg *.log *.dvi *.out *.toc

clean:	texclean
	-rm $(PDFNAME)

.PHONY: all texclean clean
