#!/bin/bash

# Define the number of parallel jobs you want to run.
# A good starting point is the number of cores on your machine.
MAX_JOBS=$(nproc)

# Define the bin ranges from the C++ script
XF_BINS=17
MASS_BINS=11

echo "Starting plot generation with up to $MAX_JOBS parallel jobs..."

# Loop over all bins and launch the ROOT script for each one
for ixf in $(seq 0 $((XF_BINS - 1))); do
    for imass in $(seq 0 $((MASS_BINS - 1))); do
        echo "Queueing job for xF bin $ixf, Mass bin $imass"
        
        # Run the ROOT macro in batch mode, passing the bin numbers.
        # The '&' at the end runs the command in the background.
        root -l -b -q "Efficiency_single_bin.C($ixf, $imass)" &

        # Control the number of parallel jobs
        # When the number of background jobs reaches MAX_JOBS, wait for one to finish.
        if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
            wait -n
        fi
    done
done

# Wait for all remaining background jobs to finish
echo "Waiting for the last batch of jobs to complete..."
wait

echo "All tasks are complete! 🎉"
