#!/bin/bash

# This script runs the ROOT macro 'ProcessSingleXFBin.C' in parallel
# for all 17 xF bins, and then merges the output into a single file.

# --- CONFIGURATION ---
N_BINS=17
PROCESS_SCRIPT="Efficiency_xF.C"
MERGE_SCRIPT="merge_files.C"
FINAL_OUTPUT_FILE="AllGraphs_D2_Efficiency.root"


# --- STEP 1: PARALLEL PROCESSING ---

# Check if the ROOT macro file exists
if [ ! -f "$PROCESS_SCRIPT" ]; then
    echo "Error: Processing script '$PROCESS_SCRIPT' not found!"
    exit 1
fi

echo "### STEP 1: Starting parallel processing for $N_BINS bins..."
echo "----------------------------------------------------------"

# Loop from 0 to (N_BINS - 1) and launch jobs in the background
for i in $(seq 0 $((N_BINS - 1)))
do
   root -l -b -q "${PROCESS_SCRIPT}(${i})" &
done

# Wait for all background jobs to complete
echo
echo "All jobs launched. Waiting for completion..."
wait
echo "----------------------------------------------------------"
echo "✅ All parallel jobs have finished."
echo


# --- STEP 2: MERGING RESULTS ---

# Check if merge script exists
if [ ! -f "$MERGE_SCRIPT" ]; then
    echo "Error: Merge script '$MERGE_SCRIPT' not found!"
    exit 1
fi

echo "### STEP 2: Merging individual ROOT files..."
echo "----------------------------------------------------------"
root -l -b -q "${MERGE_SCRIPT}(\"${FINAL_OUTPUT_FILE}\")"
echo


# --- STEP 3: CLEANUP ---
echo "### STEP 3: Cleaning up intermediate files..."
echo "----------------------------------------------------------"
rm -f D2_occ/D2_Efficiency_xF_*.root
# You can uncomment the lines below to also remove the image files
# rm -f D2_occ/D2_Efficiency_xF_*.png
# rm -f D2_occ/D2_Efficiency_xF_*.eps
echo "Cleanup complete."
echo

echo "🎉 All done! Final combined results are in '$FINAL_OUTPUT_FILE'."
