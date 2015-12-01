HOCKING-supervised-error-bands.pdf: HOCKING-supervised-error-bands.tex figure-labels-5.png refs.bib
	rm -f *.aux *.bbl
	pdflatex HOCKING-supervised-error-bands
	bibtex HOCKING-supervised-error-bands
	pdflatex HOCKING-supervised-error-bands
	pdflatex HOCKING-supervised-error-bands
figure-labels-5.png: figure-labels.R data.by.type.RData
	R --no-save < $<
data.by.type.RData: data.by.type.R pairs.txt
	R --no-save < $<