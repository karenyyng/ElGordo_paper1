PROJECT=ElGordo1
TEX=pdflatex
BIBTEX=bibtex
BUILDTEX=$(TEX) $(PROJECT).tex

# tell make that all, clean and test do not explicitly depend on files
.PHONY: all clean test

# tell make to test if finally product exist first
# avoid rebuilding if no changes are made to the tex files
test: $(PROJECT).pdf 

all:
	$(BUILDTEX)
	$(BIBTEX) $(PROJECT)
	$(BIBTEX) $(PROJECT)
	$(BUILDTEX)
	$(BUILDTEX)

clean-all:
	rm -f *.dvi *.log *.bak *.aux *.bbl *.blg *.idx *.ps *.eps *.pdf *.toc *.out *~

clean:
	rm -f *.log *.bak *.aux *.bbl *.blg *.idx *.toc *.out *~
