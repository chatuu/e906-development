#!/bin/bash

# Define the input directory and arrays for Roadsets/Targets
INPUT_DIR="/root/github/e906-development/src/kTrackerEfficiency/RS57-70/RS70"

ROADSETS=("RS70")
TARGETS=("LH2" "LD2" "Flask")

# Loop over the arrays and execute the python script
for RS in "${ROADSETS[@]}"; do
    for TARGET in "${TARGETS[@]}"; do
        
        INPUT_FILE="${INPUT_DIR}/merged_${RS}_${TARGET}_recoeff.root"
        OUTPUT_FILE="merged_${RS}_${TARGET}_recoeff_hodoeff.root"
        
        # Check if the input file actually exists before running
        if [ ! -f "$INPUT_FILE" ]; then
            echo "Error: Input file '$INPUT_FILE' not found. Skipping $TARGET for $RS..."
            echo "--------------------------------------------------"
            continue
        fi

        echo "Executing python script for Roadset: $RS, Target: $TARGET..."
        python3 AddHodoEffVarsToROOTFile.py \
            --input "$INPUT_FILE" \
            --output "$OUTPUT_FILE" \
            --target "$TARGET"
            
        echo "Finished $TARGET for $RS."
        echo "--------------------------------------------------"
        
    done
done

echo "All files processed successfully."