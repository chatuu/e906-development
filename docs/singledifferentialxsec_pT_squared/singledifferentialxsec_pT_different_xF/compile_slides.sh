#!/bin/bash

# 1. Generate the .tex file
echo "Generating LaTeX slides..."
python3 generate_slides.py

# 2. Compile LaTeX
# We compile twice to ensure proper indexing/navigation
echo "Compiling LaTeX (Run 1)..."
pdflatex -interaction=nonstopmode slides.tex

echo "Compiling LaTeX (Run 2)..."
pdflatex -interaction=nonstopmode slides.tex

# 3. Cleanup auxiliary files
echo "Cleaning up..."
rm *.aux *.log *.nav *.out *.snm *.toc

echo "Done! File 'slides.pdf' has been generated."