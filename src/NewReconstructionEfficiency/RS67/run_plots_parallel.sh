#!/bin/bash

# Name of your python script
PYTHON_SCRIPT="generate_reco_eff_plots.py"

# Array of arguments for each job
declare -a COMMANDS=(
    "-i ./LH2/merged_RS67_LH2_recoeff_unfolding_new.root -t LH2 --tree result -n ./GlobalEfficiencyCurve/rs67_lh2_eff_D1.npz -o ./reco_eff_hists_LH2.root"
    "-i ./LD2/merged_RS67_LD2_recoeff_unfolding_new.root -t LD2 --tree result -n ./GlobalEfficiencyCurve/rs67_ld2_eff_D1.npz -o ./reco_eff_hists_LD2.root"
    "-i ./Flask/merged_RS67_Empty_recoeff_unfolding_new.root -t Flask --tree result -n ./GlobalEfficiencyCurve/rs67_avg_eff_D1.npz -o ./reco_eff_hists_Flask.root"
    "-i ./LH2/merged_RS67_LH2_recoeff_unfolding_new.root -t LH2 --tree result_mix -n ./GlobalEfficiencyCurve/rs67_lh2_eff_D1.npz -o ./reco_eff_hists_LH2_mix.root"
    "-i ./LD2/merged_RS67_LD2_recoeff_unfolding_new.root -t LD2 --tree result_mix -n ./GlobalEfficiencyCurve/rs67_ld2_eff_D1.npz -o ./reco_eff_hists_LD2_mix.root"
    "-i ./Flask/merged_RS67_Empty_recoeff_unfolding_new.root -t Flask --tree result_mix -n ./GlobalEfficiencyCurve/rs67_avg_eff_D1.npz -o ./reco_eff_hists_Flask_mix.root"
)

echo "======================================================"
echo "Starting parallel generation of reconstruction plots..."
echo "======================================================"

# Loop over the commands and launch each in the background
for ARGS in "${COMMANDS[@]}"; do
    echo "Submitting: python $PYTHON_SCRIPT $ARGS"
    # The '&' at the end sends the process to the background
    python3 $PYTHON_SCRIPT $ARGS &
done

echo "------------------------------------------------------"
echo "All 6 jobs submitted. Waiting for them to finish..."
echo "------------------------------------------------------"

# 'wait' pauses the script until all background jobs (&) have completed
wait

echo "======================================================"
echo "All plots generated successfully!"
echo "======================================================"
