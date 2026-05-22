#!/bin/bash

# Define the working directory
#WORK_DIR=~/github/e906-development/src/xsec_pT

# Navigate to directory and exit if it fails
#cd "$WORK_DIR" || { echo "Failed to navigate to $WORK_DIR"; exit 1; }

echo "Generating LaTeX code..."
python3 generate_slides.py

# Check if the tex file was successfully generated
if [ ! -f yield_slides.tex ]; then
    echo "Error: yield_slides.tex was not created."
    exit 1
fi

echo "Compiling Beamer slides with pdflatex..."
pdflatex -interaction=nonstopmode yield_slides.tex

echo "Done! The compiled slides: yield_slides.pdf are ready"