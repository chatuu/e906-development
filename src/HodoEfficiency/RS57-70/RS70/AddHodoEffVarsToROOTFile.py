import sys
import os

# Add the parent directory to the Python path to import functions.py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions import *

import gc
import numpy as np
import pandas as pd
import ROOT
import argparse

# Disable ROOT graphical processing to save resources
ROOT.gROOT.SetBatch(True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--xf_min', type=float, help="Minimum xF to keep (optional)")
    parser.add_argument('--xf_max', type=float, help="Maximum xF to keep (optional)")
    parser.add_argument('--input', type=str, required=True, help="Path to input ROOT file")
    parser.add_argument('--output', type=str, required=True, help="Name of output ROOT file")
    parser.add_argument('--target', type=str, required=True, choices=["LH2", "LD2", "Flask"], help="Target type (LH2, LD2, Flask)")
    args = parser.parse_args()

    xsec_file_path = args.input
    
    print(f"\n==========================================")
    print(f"Processing Target: {args.target}")
    print(f"Input File: {xsec_file_path}")
    print(f"Output File: {args.output}")
    print(f"==========================================")

    # --- Data Filtering Logic ---
    # 1. Get the base Chuck Cuts as a C++ string
    chuck_cuts_str = get_e906_chuck_cuts_string()
    
    # 2. Combine with optional xF cuts
    final_filter_str = chuck_cuts_str
    if args.xf_min is not None and args.xf_max is not None:
        print(f"Config: Filtering data for xF in range [{args.xf_min}, {args.xf_max}]")
        final_filter_str = f"({chuck_cuts_str}) && (xF >= {args.xf_min} && xF < {args.xf_max})"
    else:
        print("Config: No additional xF cuts provided. Using base Chuck Cuts.")

    # --- Processing Loop ---
    # We want to process both 'result' and 'result_mix'
    trees_to_process = ["result", "result_mix"]
    
    # First tree RECREATEs the file, subsequent trees UPDATE it
    file_mode = "RECREATE"

    for tree_name in trees_to_process:
        print(f"\n--- Processing Tree: {tree_name} ---")
        
        # 1. Dynamically get ALL column names from the tree
        try:
            # We open a temporary RDataFrame just to read the column list
            rdf_temp = ROOT.RDataFrame(tree_name, xsec_file_path)
            all_cols_vec = rdf_temp.GetColumnNames()
            # Convert std::vector<string> to a python list of strings
            xsec_cols = [str(c) for c in all_cols_vec]
            print(f"Detected {len(xsec_cols)} columns in {tree_name}.")
        except Exception as e:
            print(f"Error reading column names for {tree_name}: {e}")
            continue

        print(f"Loading and filtering data via RDataFrame from: {xsec_file_path}")
        
        # This now handles BOTH loading and the Chuck Cuts at the C++ level
        df_passed = load_data_from_root(xsec_file_path, tree_name, xsec_cols, filter_str=final_filter_str)
        
        if df_passed.empty:
            print(f"Skipping {tree_name} (no events passed cuts or tree not found).")
            continue

        print(f"Data successfully loaded and cut. Rows remaining: {len(df_passed)}")

        # Calculate efficiencies (adds 'hodoeff' and 'hodoeff_error')
        df_passed = calculate_hodo_efficiency_columns(df_passed)
        
        # Save to ROOT file
        save_dataframe_to_root(df_passed, args.output, tree_name, mode=file_mode)
        
        # Ensure next iteration updates the file rather than overwriting
        file_mode = "UPDATE"
        
        # Cleanup
        del df_passed
        gc.collect()
            
    print("\nProcessing complete.")

if __name__ == "__main__":
    main()