#!/bin/bash

TEX_FILE="hodoscope_slides.tex"
OUTPUT_PDF="hodoscope_slides.pdf"

if [ ! -f "$TEX_FILE" ]; then
    echo "Error: $TEX_FILE not found. Please run the python slide generator script first."
    exit 1
fi

echo "Compiling $TEX_FILE using pdflatex..."

pdflatex -interaction=nonstopmode "$TEX_FILE"

if [ $? -eq 0 ]; then
    echo "----------------------------------------"
    echo "Success! Presentation generated: $OUTPUT_PDF"
    echo "Cleaning up auxiliary files..."
    rm -f *.aux *.log *.nav *.out *.snm *.toc
else
    echo "----------------------------------------"
    echo "LaTeX compilation encountered errors. Check the console output above."
fi
