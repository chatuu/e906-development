#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "Cleaning up old LaTeX artifacts..."
# Safely delete ONLY the generated LaTeX source, PDFs, and build clutter. 
# This preserves your actual plot PDFs!
#rm -f slides.tex slides.pdf slides.aux slides.log slides.out slides.nav slides.snm slides.toc

echo "Running Python script to generate new .tex files..."
python3 generate_latex.py

# Check if the slides.tex file was generated successfully
if [[ ! -f slides.tex ]]; then
    echo "Error: slides.tex was not generated."
    exit 1
fi

echo "Compiling slides.tex (Pass 1)..."
pdflatex -interaction=nonstopmode slides.tex

echo "Compiling slides.tex (Pass 2 for total slide counts/references)..."
pdflatex -interaction=nonstopmode slides.tex

echo "Cleaning up LaTeX build clutter..."
rm -f *.aux *.log *.out *.nav *.snm *.toc

echo "========================================"
echo "✔ Success! slides.pdf is ready."
echo "========================================"