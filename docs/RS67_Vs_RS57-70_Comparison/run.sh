#!/bin/bash

echo "Cleaning old files..."
rm -f presentation.aux presentation.log presentation.nav presentation.out presentation.snm presentation.toc presentation.tex

echo "Generating LaTeX..."
python3 generate_slides.py

echo "Compiling (Pass 1)..."
# Removed > /dev/null and added -halt-on-error so it tells us exactly what's wrong
pdflatex -interaction=nonstopmode -halt-on-error presentation.tex

# If Pass 1 fails, the script will exit before running Pass 2
if [ $? -ne 0 ]; then
    echo "========================================================"
    echo "LaTeX compilation failed! Check the error message above."
    echo "========================================================"
    exit 1
fi

echo "Compiling (Pass 2)..."
pdflatex -interaction=nonstopmode -halt-on-error presentation.tex

echo "Success! presentation.pdf is ready."