#!/bin/bash

# Define the target tex file
TEX_FILE="dy_cross_section_slides.tex"

# Ensure the python script is run first to generate the tex file
echo "Generating LaTeX source..."
python3 generate_slides.py

echo "Compiling $TEX_FILE..."

# First pass (generates layout and aux files)
pdflatex -interaction=nonstopmode $TEX_FILE

# Second pass (resolves table of contents and navigation numbering)
pdflatex -interaction=nonstopmode $TEX_FILE

# Clean up auxiliary build files
echo "Cleaning up compilation logs..."
rm -f *.aux *.log *.nav *.out *.snm *.toc

echo "Done! Presentation saved as dy_cross_section_slides.pdf"