#!/bin/bash

REPORT_TEX="efficiency_corrections_comparison.tex"
SLIDES_TEX="efficiency_corrections_slides.tex"

echo "=========================================="
echo "Generating LaTeX comparison scripts..."
echo "=========================================="
python3 generate_eff_comparison.py

echo ""
echo "=========================================="
echo "Compiling LaTeX Report (Article)..."
echo "=========================================="
if [ -f "$REPORT_TEX" ]; then
    # Run pdflatex twice for Table of Contents and List of Figures
    pdflatex -interaction=nonstopmode "$REPORT_TEX"
    pdflatex -interaction=nonstopmode "$REPORT_TEX"
    echo "Report compiled successfully."
else
    echo "Error: $REPORT_TEX was not generated."
fi

echo ""
echo "=========================================="
echo "Compiling LaTeX Slides (Beamer)..."
echo "=========================================="
if [ -f "$SLIDES_TEX" ]; then
    # Run pdflatex twice to ensure Beamer formatting and slide numbers resolve
    pdflatex -interaction=nonstopmode "$SLIDES_TEX"
    pdflatex -interaction=nonstopmode "$SLIDES_TEX"
    echo "Slides compiled successfully."
else
    echo "Error: $SLIDES_TEX was not generated."
fi

echo ""
echo "=========================================="
echo "Done! Outputs generated:"
echo " - efficiency_corrections_comparison.pdf (Report)"
echo " - efficiency_corrections_slides.pdf (Presentation)"
echo "=========================================="