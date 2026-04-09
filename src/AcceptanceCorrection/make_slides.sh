#!/bin/bash

# Exit immediately if any command fails
set -e

# Define file names
PYTHON_SCRIPT="generate_slides.py"
TEX_FILE="acceptance_slides.tex"

echo "========================================"
echo "1. Generating LaTeX file via Python..."
echo "========================================"
python3 "$PYTHON_SCRIPT"

echo ""
echo "========================================"
echo "2. Compiling LaTeX to PDF..."
echo "========================================"
# The -interaction=nonstopmode flag ensures it won't hang waiting for user input if there's a LaTeX warning
pdflatex -interaction=nonstopmode "$TEX_FILE"

echo ""
echo "========================================"
echo "Done! Your slides are ready: acceptance_slides.pdf"
echo "========================================"