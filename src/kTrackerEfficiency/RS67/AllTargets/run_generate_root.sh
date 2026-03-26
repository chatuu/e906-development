#!/bin/bash

# Define arrays for the targets, input files, and output files
TARGETS=("LH2" "LD2" "Flask")

INPUTS=(
    "/root/github/e906-development/ROOTFiles/Chatura/All_LH2_Data_Mix.root"
    "/root/github/e906-development/ROOTFiles/Chatura/All_LD2_Data_Mix.root"
    "/root/github/e906-development/ROOTFiles/Chatura/All_Flask_Data_Mix.root"
)

OUTPUTS=(
    "merged_RS67_3089_LH2_recoeff.root"
    "merged_RS67_3089_LD2_recoeff.root"
    "merged_RS67_3089_Flask_recoeff.root"
)

# Loop over the arrays and execute the python script
for i in "${!TARGETS[@]}"; do
    TARGET="${TARGETS[$i]}"
    INPUT_FILE="${INPUTS[$i]}"
    OUTPUT_FILE="${OUTPUTS[$i]}"
    
    # Check if the input file actually exists before running
    if [ ! -f "$INPUT_FILE" ]; then
        echo "Error: Input file '$INPUT_FILE' not found. Skipping $TARGET..."
        echo "--------------------------------------------------"
        continue
    fi

    echo "Executing python script for $TARGET..."
    python3 GenerateROOTFiles.py \
        --input "$INPUT_FILE" \
        --output "$OUTPUT_FILE" \
        --target "$TARGET"
        
    echo "Finished $TARGET."
    echo "--------------------------------------------------"
done

echo "All files processed successfully."
