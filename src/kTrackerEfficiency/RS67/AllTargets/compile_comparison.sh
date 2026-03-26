#!/bin/bash

TEX_FILE="efficiency_comparison.tex"
PDF_FILE="efficiency_comparison.pdf"

echo "=========================================="
echo "Generating LaTeX comparison script..."
echo "=========================================="
python3 generate_comparison_tex.py

# Check if the python script successfully created the .tex file
if [ -f "$TEX_FILE" ]; then
    echo ""
    echo "=========================================="
    echo "Compiling LaTeX document..."
    echo "=========================================="
    
    # Run pdflatex with nonstopmode so it doesn't hang on minor warnings
    pdflatex -interaction=nonstopmode "$TEX_FILE"
    
    # Run a second time to ensure the Table of Contents generates correctly
    pdflatex -interaction=nonstopmode "$TEX_FILE"
    
    echo ""
    echo "=========================================="
    if [ -f "$PDF_FILE" ]; then
        echo "Success! Comparison document saved as: $PDF_FILE"
    else
        echo "Error: PDF compilation failed. Check the .log file."
    fi
    echo "=========================================="
    
    # Optional: Clean up auxiliary LaTeX files (uncomment if you want to keep the directory clean)
    # rm -f *.aux *.log *.toc
else
    echo "Error: $TEX_FILE was not generated. Check the Python script."
fi
