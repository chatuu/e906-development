#!/bin/bash

# Exit immediately if any command fails
set -e

echo "==========================================================="
echo " Starting Drell-Yan Cross-Section & Unfolding Pipeline"
echo "==========================================================="
echo ""

# Step 1: Run the main analysis and bootstrapping
echo ">>> [1/2] Executing main.py..."
echo "    (Extracting kinematics, building matrices, unfolding, and bootstrapping)"
python3 main.py

echo ""
echo ">>> main.py completed successfully."
echo ""

# Step 2: Extract and calculate the systematics
echo ">>> [2/2] Executing extract_unfolding_systematics.py..."
echo "    (Processing 100 toys, building histograms, extracting mean systematics)"
python3 extract_unfolding_systematics.py

echo ""
echo "==========================================================="
echo " ✔ Pipeline execution finished successfully!"
echo " Key outputs generated:"
echo "  - All_XSec_Objects.root"
echo "  - DY_ResponseMatrices.root"
echo "  - Unfolding_Sys_Hists.root"
echo "  - Unfolding_Systematics_Bootstrap.csv"
echo "  - Response_Matrices_PDFs/ (Directory)"
echo "  - Systematics_Histograms/ (Directory)"
echo "==========================================================="