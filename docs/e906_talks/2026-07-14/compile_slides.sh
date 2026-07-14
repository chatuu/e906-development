#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

TEX_FILE="presentation.tex"

echo "1. Executing Python script to generate Beamer framework..."
python3 generate_slides.py

if [ -f "$TEX_FILE" ]; then
    echo "2. Compiling LaTeX document to PDF..."
    # Running twice handles Beamer's internal reference compilation properly
    pdflatex -interaction=nonstopmode $TEX_FILE
    pdflatex -interaction=nonstopmode $TEX_FILE
    
    echo "Done. Compiled presentation is available as presentation.pdf."
else
    echo "Error: Python script failed to generate $TEX_FILE."
    exit 1
fi