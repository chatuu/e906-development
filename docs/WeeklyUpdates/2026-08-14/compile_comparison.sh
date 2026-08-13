#!/bin/bash

# Define filenames
PY_SCRIPT="generate_comparison_slides.py"
TEX_FILE="Unfolding_Comparison.tex"
PDF_FILE="Unfolding_Comparison.pdf"

echo "=========================================="
echo "1. Generating LaTeX source via Python..."
echo "=========================================="
python3 $PY_SCRIPT

if [ ! -f "$TEX_FILE" ]; then
    echo "Error: LaTeX file was not generated."
    exit 1
fi

echo -e "\n=========================================="
echo "2. Compiling PDF (Pass 1)..."
echo "=========================================="
pdflatex -interaction=nonstopmode $TEX_FILE > pdflatex.log

if [ $? -ne 0 ]; then
    echo "FAILED: LaTeX Compilation Error on Pass 1."
    echo "Here are the last 15 lines of the error log:"
    echo "------------------------------------------"
    tail -n 15 pdflatex.log
    exit 1
fi

echo -e "\n=========================================="
echo "3. Compiling PDF (Pass 2)..."
echo "=========================================="
pdflatex -interaction=nonstopmode $TEX_FILE > pdflatex.log

if [ $? -ne 0 ]; then
    echo "FAILED: LaTeX Compilation Error on Pass 2."
    echo "Here are the last 15 lines of the error log:"
    echo "------------------------------------------"
    tail -n 15 pdflatex.log
    exit 1
fi

# Clean up auxiliary files
echo -e "\n=========================================="
echo "4. Cleaning up..."
echo "=========================================="
rm -f *.aux *.log *.nav *.out *.snm *.toc

if [ -f "$PDF_FILE" ]; then
    echo "SUCCESS: $PDF_FILE has been generated."
else
    echo "FAILED: PDF was not found."
fi