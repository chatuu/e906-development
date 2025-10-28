#!/bin/bash

# This script compiles and runs the efficiency calculation in parallel
# for all 11 mass bins, and then merges the results into a single file.

# Define the C++ source file and the output executable name
SOURCE_FILE="Efficiency_Mass_Parallel.cpp"
EXECUTABLE="runEfficiency"
MERGE_SCRIPT="merge_graphs.C"

# Number of bins (NBINS = 11, so index is 0-10)
MAX_BIN_INDEX=10

# Output directory for plots, ROOT files, and logs
OUTPUT_DIR="D2_occ"
LOG_DIR="${OUTPUT_DIR}/logs"

# --- 1. Compilation ---
echo "Compiling ${SOURCE_FILE}..."
g++ -Wall -O2 ${SOURCE_FILE} -o ${EXECUTABLE} $(root-config --cflags --glibs)

if [ $? -ne 0 ]; then
    echo "Compilation failed! Exiting."
    exit 1
fi
echo "Compilation successful. Executable is '${EXECUTABLE}'."

# --- 2. Create Output Directories ---
echo "Creating output directories..."
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# --- 3. Parallel Execution ---
echo "Starting parallel jobs for bins 0 to ${MAX_BIN_INDEX}..."
for i in $(seq 0 ${MAX_BIN_INDEX}); do
    echo "  -> Launching job for mass bin ${i}"
    ./${EXECUTABLE} ${i} > "${LOG_DIR}/log_bin_${i}.txt" 2>&1 &
done

# --- 4. Wait ---
echo "All jobs launched. Waiting for all background jobs to complete..."
wait
echo "All parallel jobs finished."

# --- 5. Merging Step (NEW) ---
echo "Merging individual ROOT files..."
if [ -f "$MERGE_SCRIPT" ]; then
    # Run the merge macro in batch mode (-b)
    root -l -b -q "${MERGE_SCRIPT}()"
    echo "Merging complete. Final file is 'EfficiencyCurves_All_Mass_bins.root'"
else
    echo "Error: Merge script '$MERGE_SCRIPT' not found! Skipping merge."
fi

# --- 6. Cleanup (Optional - uncomment to enable) ---
# echo "Cleaning up intermediate .root files..."
# rm ${OUTPUT_DIR}/D2_Efficiency_M*.root
# echo "Cleanup complete."

echo "Script finished."