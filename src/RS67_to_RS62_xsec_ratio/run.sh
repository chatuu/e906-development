#!/bin/bash

rm -f *.pdf *.tex ./comparison_plots/*.pdf
python3 generatePlots.py
python3 generate_slides.py
pdflatex -interaction=nonstopmode -halt-on-error CrossSection_Comparison.tex
