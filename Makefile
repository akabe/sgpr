FILE=changes
LATEX=pdflatex -halt-on-error -interaction=nonstopmode

all:
	$(LATEX) $(FILE)
	$(LATEX) $(FILE)

clean:
	rm -rf *.acn *.acr *.alg *.aux *.bbl *.blg *.dvi *.fdb_latexmk *.glg *.glo *.gls *.idx *.ilg *.ind *.ist *.lof *.log *.lot *.maf *.mtc *.mtc0 *.nav *.nlo *.out *.pdfsync *.ps *.snm *.synctex.gz *.toc *.vrb *.xdy *.tdo
