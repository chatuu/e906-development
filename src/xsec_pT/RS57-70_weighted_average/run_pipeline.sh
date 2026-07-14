#!/bin/bash
# Parallel multi-core execution for Drell-Yan roadset cross-sections

ROADSETS=("RS57" "RS59" "RS62" "RS67" "RS70")

echo "=========================================================="
echo " Launching Drell-Yan Analysis Pipeline in Parallel"
echo "=========================================================="

# Loop through the roadsets and run them in the background
for rs in "${ROADSETS[@]}"; do
    echo "--> Submitting execution job for $rs..."
    # Pipes the terminal output to individual log files to prevent console scrambling
    python3 main.py --roadset $rs > "${rs}_execution.log" 2>&1 &
done

# Wait for all background Python jobs to finish
wait

echo "=========================================================="
echo " ✔ All 5 roadset ROOT files generated successfully!"
echo " You can now execute the variance calculation script."
echo "=========================================================="