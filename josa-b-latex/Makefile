LATEXMK ?= latexmk
LATEXFLAGS = -pdf -interaction=nonstopmode -halt-on-error

.PHONY: all

# Main uses the publisher's BibTeX style; the supplement uses biber.
# Build sequentially because each document imports labels from the other.
all:
	$(LATEXMK) $(LATEXFLAGS) supplement.tex
	python3 export-references.py supplement
	$(LATEXMK) $(LATEXFLAGS) fwm-and-coupler.tex
	python3 export-references.py fwm-and-coupler
	$(LATEXMK) $(LATEXFLAGS) supplement.tex
	python3 export-references.py supplement
	$(LATEXMK) $(LATEXFLAGS) fwm-and-coupler.tex
	python3 export-references.py fwm-and-coupler
