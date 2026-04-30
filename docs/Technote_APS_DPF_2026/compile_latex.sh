#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "Cleaning up old LaTeX artifacts..."
# Safely delete ONLY the generated LaTeX source, PDFs, and build clutter. 
# This preserves your actual plot PDFs!
rm -f slides.tex document.tex slides.pdf document.pdf \
      slides.aux document.aux slides.log document.log \
      slides.out document.out slides.nav slides.snm \
      slides.toc document.toc

echo "Running Python script to generate new .tex files..."
python3 generate_latex.py

# Check if the .tex files were generated successfully
if [[ ! -f slides.tex || ! -f document.tex ]]; then
    echo "Error: .tex files were not generated."
    exit 1
fi

echo "Compiling slides.tex (Pass 1)..."
pdflatex -interaction=nonstopmode slides.tex > /dev/null

echo "Compiling slides.tex (Pass 2 for formatting/references)..."
pdflatex -interaction=nonstopmode slides.tex > /dev/null

echo "Compiling document.tex (Pass 1)..."
pdflatex -interaction=nonstopmode document.tex > /dev/null

echo "Compiling document.tex (Pass 2 for formatting/references)..."
pdflatex -interaction=nonstopmode document.tex > /dev/null

echo "Cleaning up LaTeX build clutter..."
rm -f *.aux *.log *.out *.nav *.snm *.toc

echo "========================================"
echo "✔ Success! slides.pdf and document.pdf are ready."
echo "========================================"