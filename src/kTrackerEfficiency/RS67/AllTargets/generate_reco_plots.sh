#!/bin/bash

# Define arrays for the targets and their corresponding input files
TARGETS=("LH2" "LD2" "Flask")
FILES=(
    "merged_RS67_3089_LH2_recoeff.root"
    "merged_RS67_3089_LD2_recoeff.root"
    "merged_RS67_3089_Flask_recoeff.root"
)

# Loop over the arrays and execute the python script
for i in "${!TARGETS[@]}"; do
    TARGET="${TARGETS[$i]}"
    FILE="${FILES[$i]}"
    
    # Check if the file exists before running
    if [ ! -f "$FILE" ]; then
        echo "Error: File '$FILE' not found. Skipping $TARGET..."
        echo "--------------------------------------------------"
        continue
    fi

    echo "Running script for $TARGET..."
    python3 generate_reco_hists.py --input "$FILE" --target "$TARGET"
    
done

echo "All targets processed successfully."
