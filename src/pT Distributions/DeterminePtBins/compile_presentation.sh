#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

TEX_FILE="comparison_slides.tex"

echo "=== Step 1: Generating LaTeX source code ==="
python3 generate_slides.py

echo -e "\n=== Step 2: Compiling PDF (Pass 1) ==="
pdflatex -interaction=nonstopmode $TEX_FILE

echo -e "\n=== Step 3: Compiling PDF (Pass 2) ==="
# Running pdflatex a second time resolves total page count and layout references
pdflatex -interaction=nonstopmode $TEX_FILE

echo -e "\n=== Clean up ==="
# Optional: Remove auxiliary files created by LaTeX to keep your directory clean
rm -f *.aux *.log *.nav *.out *.snm *.toc

echo -e "\n=== Done! ==="
echo "Your presentation is ready: comparison_slides.pdf"