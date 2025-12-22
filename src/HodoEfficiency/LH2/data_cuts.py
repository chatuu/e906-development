from functions import *
import numpy as np
import pandas as pd
import ROOT
import argparse
import sys
import gc

# Disable ROOT graphical processing to save resources
ROOT.gROOT.SetBatch(True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--xf_min', type=float)
    parser.add_argument('--xf_max', type=float)
    args = parser.parse_args()

    mass_bins = [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7]

    xsec_file_path = "/root/github/e906-development/ROOTFiles/MixedEvents/merged_RS67_3089LH2.root"
    xsec_cols = list(set(get_cut_columns() + ['runID', 'spillID', 'eventID']))
    
    # --- MEMORY OPTIMIZATION STEP ---
    # Construct a filter string so we ONLY load the xF bin we currently care about
    xf_filter_str = ""
    if args.xf_min is not None and args.xf_max is not None:
        xf_filter_str = f"xF >= {args.xf_min} && xF < {args.xf_max}"
    
    df_xsec = load_data_from_root(xsec_file_path, "result", xsec_cols, xf_filter_str)
    
    if df_xsec.empty:
        return

    # Apply cuts to this smaller dataframe
    cut_mask = apply_data_cuts(df_xsec)
    df_xsec_passed = df_xsec[cut_mask].copy()

    # Clear memory of the raw dataframe
    del df_xsec
    gc.collect()

    # Load Hugo Tree (Roads)
    # We can't filter this easily by xF (it doesn't have it), so we load IDs + Roads
    hugo_file_path = "/root/github/e906-development/ROOTFiles/Hugo/roadset57_70_R008_2111v42_tmp_noPhys.root"
    hugo_cols = ['runID', 'spillID', 'eventID', 'posRoad', 'negRoad']
    
    df_hugo = load_data_from_root(hugo_file_path, "Tree", hugo_cols)
    if df_hugo.empty:
        return

    # Merge
    df_merged = pd.merge(df_xsec_passed, df_hugo, on=['runID', 'spillID', 'eventID'], how='inner')

    # Iterate Mass Bins
    if args.xf_min is not None:
        # We are processing one specific xF bin (passed via CLI)
        # We iterate through all Mass bins for this slice
        for i in range(len(mass_bins) - 1):
            process_single_bin(df_merged, mass_bins[i], mass_bins[i+1], args.xf_min, args.xf_max)
    else:
        # Legacy mode (no args)
        pass

if __name__ == "__main__":
    main()