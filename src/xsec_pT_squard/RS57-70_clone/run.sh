#!/bin/bash

# Exit immediately if any command exits with a non-zero status
set -e

echo "Cleaning up old output files..."
# Removes all previous run artifacts to ensure a clean slate
rm -f *.pdf *.csv *.tex *.root

echo "Running Drell-Yan cross-section analysis..."
python3 main.py

echo "Applying SeaQuest logo to CrossSection PDFs..."
# Only process if CrossSection PDFs were generated successfully
if ls CrossSection_*.pdf 1> /dev/null 2>&1; then
    ls -1 CrossSection_*.pdf | xargs -I {} python3 addlogo.py {} SeaQuestLogo.png
    echo "Logos applied successfully."
    
    echo "Deleting original un-logoed CrossSection PDFs..."
    # Finds and deletes CrossSection PDFs that DO NOT have '_with_logo' in the name
    find . -maxdepth 1 -type f -name "CrossSection_*.pdf" ! -name "*_with_logo.pdf" -delete
else
    echo "Warning: No CrossSection_*.pdf files found to process."
fi

echo "Workflow complete!"