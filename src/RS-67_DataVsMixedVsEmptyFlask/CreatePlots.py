# -*- coding: utf-8 -*-
"""
Plot mass distribution for specific datasets after a target cut,
by re-processing ROOT files in parallel.
"""

import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import re
from collections import OrderedDict # To ensure cut order
import concurrent.futures
import os # For os.getpid() and os.cpu_count()
import time

# --- Configuration ---

TARGET_LABELS_FOR_GOOD_SPILLS_CUT_CONFIG = {"Data", "Data(RS67)", "empty flask (RS67)"}
GOOD_SPILLS_FILE_PATH = "../../res/GoodSpills/GoodSpills.txt"
TARGET_CUT_NAME = "D1 + D2 + D3 < 1000"

try:
    # Adjust this path based on the script's location relative to 'res'
    base_path_mixed_events = "../../res/ROOTFiles/MixedEvents/"
except NameError:
    print("Warning: base_path_mixed_events not defined, using placeholder './'", file=sys.stderr)
    base_path_mixed_events = "./"

FILE_PATHS_CONFIG = {
    "Data(RS67)": {
        "path": base_path_mixed_events + "merged_RS67_3089LH2.root",
        "tree": "result"
    },
    "Mixed(RS67)": {
        "path": base_path_mixed_events + "merged_RS67_3089LH2.root",
        "tree": "result_mix"
    },
    "empty flask (RS67)": {
        "path": base_path_mixed_events + "merged_RS67_3089flask.root",
        "tree": "result"
    },
    "empty flask (RS67) mixed": {
        "path": base_path_mixed_events + "merged_RS67_3089flask.root",
        "tree": "result_mix"
    }
}

# Full ordered dictionary of cuts
cuts_dict_full = OrderedDict([
    ("Good Spills Cut", "spillID # Custom handled"), 
    ("z1_v and z2_v within -320cm to -5cm", "z1_v > -320 and z1_v < -5 and z2_v > -320 and z2_v < -5"),
    ("xt^2 + (yt - beamOffset)^2 < 320cm^2", "((x1_t**2 + (y1_t - beamOffset)**2 < 320)) and ((x2_t**2 + (y2_t - beamOffset)**2 < 320))"),
    ("xd^2 + (yd - beamOffset)^2 within 16cm^2 to 1100cm^2", "((x1_d**2 + (y1_d - beamOffset)**2 > 16) and (x1_d**2 + (y1_d - beamOffset)**2 < 1100)) and ((x2_d**2 + (y2_d - beamOffset)**2 > 16) and (x2_d**2 + (y2_d - beamOffset)**2 < 1100))"),
    ("abs(abs(px_st1 - px_st3) - 0.416) < 0.008GeV", "abs(abs(px1_st1 - px1_st3) - 0.416) < 0.008 and abs(abs(px2_st1 - px2_st3) - 0.416) < 0.008"),
    ("abs(py_st1 - py_st3) < 0.008GeV", "abs(py1_st1 - py1_st3) < 0.008 and abs(py2_st1 - py2_st3) < 0.008"),
    ("abs(pz_st1 - pz_st3) < 0.08GeV", "abs(pz1_st1 - pz1_st3) < 0.08 and abs(pz2_st1 - pz2_st3) < 0.08"),
    ("chisq1_target and chisq2_target < 15", "chisq1_target < 15 and chisq2_target < 15"),
    ("nHits Track1 and Track2 > 13", "nHits1 > 13 and nHits2 > 13"),
    ("chisq/NDF < 12", "(chisq1/(nHits1-5)) < 12 and (chisq2/(nHits2-5)) < 12"),
    ("pz1_st1 and pz2_st1 within 9 to 75", "pz1_st1 > 9 and pz1_st1 < 75 and pz2_st1 > 9 and pz2_st1 < 75"),
    ("y1_st1/y1_st3 and y2_st1/y2_st3 < 1", "(y1_st1/y1_st3) < 1 and (y2_st1/y2_st3) < 1"),
    ("y1_st1 x y1_st3 and y2_st1 x y2_st3 > 0", "(y1_st1 * y1_st3) > 0 and (y2_st1 * y2_st3) > 0"),
    ("abs(py1_st1) and abs(py2_st1) > 0.02", "abs(py1_st1) > 0.02 and abs(py2_st1) > 0.02"),
    ("dz within -280 to -5", "dz > -280 and dz < -5"),
    ("dx within -0.25 to 0.25", "dx > -0.25 and dx < 0.25"),
    ("dy - beamOffset within -0.22 to 0.22", "abs(dy - beamOffset) < 0.22"),
    ("dx^2 + (dy - beamOffset)^2 < 0.06","(dx**2 + (dy - beamOffset)**2) < 0.06"),
    ("abs(dpx) < 1.8", "abs(dpx) < 1.8"),
    ("abs(dpy) < 2", "abs(dpy) < 2"),
    ("dpz within 38 to 116", "dpz > 38 and dpz < 116"),
    ("dpx^2 + dpy^2 < 5", "(dpx**2 + dpy**2) < 5"),
    ("mass within 4.2 to 8.8", "mass > 4.2 and mass < 8.8"),
    ("xF within -0.1 to 0.95", "xF > -0.1 and xF < 0.95"),
    ("xT within -0.1 to 0.58", "xT > -0.1 and xT < 0.58"),
    ("cosTheta within -0.5 to 0.5", "costh > -0.5 and costh < 0.5"),
    ("abs(trackSeparation) < 270", "abs(trackSeparation) < 270"),
    ("nhits1 + nhits2 > 29", "(nHits1 + nHits2) > 29"),
    ("nhits1st1 + nhits2st1 > 8", "(nHits1St1 + nHits2St1) > 8"),
    ("y1_st3 x y2_st3 < 0", "y1_st3 * y2_st3 < 0"),
    ("chisq_dimuon < 18", "chisq_dimuon < 18"),
    ("abs(x1_st1 + x2_st1) < 42", "abs(x1_st1 + x2_st1) < 42"),
    ("chisq Target within 2", "abs(chisq1_target + chisq2_target - chisq_dimuon) < 2"),
    ("D1 < 400", "D1 < 400"),
    ("D2 < 400", "D2 < 400"),
    ("D3 < 400", "D3 < 400"),
    ("D1 + D2 + D3 < 1000", "D1 + D2 + D3 < 1000"), # TARGET_CUT_NAME
    ("intensity within 0 to 80000", "intensityP > 0 and intensityP < 80000")
])

# List of variables needed for all cuts up to TARGET_CUT_NAME + 'mass'
variables_list = [
    "spillID", "mass",
    "nHits1", "nHits2", "nHits1St1", "nHits2St1",
    "chisq1", "chisq2", "chisq_dimuon", "chisq1_target", "chisq2_target",
    "dx", "dy", "dz", "dpx", "dpy", "dpz",
    "D1", "D2", "D3", "xF", "xT", "xB", "costh", "intensityP", # Keep all for potential future use, ensures subset works
    "runID", "trackSeparation",
    "x1_t", "x2_t", "y1_t", "y2_t", "x1_d", "x2_d", "y1_d", "y2_d", "z1_v", "z2_v",
    "x1_st1", "x2_st1", "y1_st1", "y2_st1", "y1_st3", "y2_st3",
    "py1_st1", "py2_st1", "py1_st3", "py2_st3",
    "pz1_st1", "pz2_st1", "pz1_st3", "pz2_st3",
    "px1_st1", "px2_st1", "px1_st3", "px2_st3",
    # "ReWeight" is intentionally omitted here to avoid reading it unless necessary
]

# --- Helper Functions ---
def load_good_spills(file_path):
    # (Function definition is unchanged from previous version)
    good_spills = set()
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                stripped_line = line.strip()
                if stripped_line:
                    try:
                        good_spills.add(int(stripped_line))
                    except ValueError:
                        print(f"⚠️ Warning: Non-integer value '{stripped_line}' in '{file_path}' line {line_num}. Skipping.", file=sys.stderr)
        if good_spills:
            print(f"✅ Loaded {len(good_spills)} good spill IDs from {file_path}")
        else:
            print(f"⚠️ No valid good spill IDs loaded from {file_path}.", file=sys.stderr)
        return good_spills
    except FileNotFoundError:
        print(f"❌ Error: Good spills file '{file_path}' not found.", file=sys.stderr)
        return set()
    except Exception as e:
        print(f"❌ Error loading good spills file '{file_path}': {e}", file=sys.stderr)
        return set()

def read_tree(file_path, tree_name, variables_to_read):
    # Modified: Does not automatically add 'ReWeight' unless in variables_to_read
    local_variables = list(set(variables_to_read + ['mass', 'spillID', 'runID', 'dy'])) # Ensure these are requested
    try:
        with uproot.open(file_path) as file:
            if tree_name not in file:
                print(f"❌ Error: Tree '{tree_name}' not found in '{file_path}'.", file=sys.stderr)
                return pd.DataFrame()
            tree = file[tree_name]
            available_branches = tree.keys()
            branches_to_fetch = [var for var in local_variables if var in available_branches]
            missing_branches = [var for var in local_variables if var not in available_branches]
            if missing_branches:
                # Filter out 'ReWeight' from missing message if it wasn't explicitly requested
                missing_to_report = [b for b in missing_branches if b != 'ReWeight' or 'ReWeight' in variables_to_read]
                if missing_to_report:
                    print(f"ℹ️ Note: Branches not found in {file_path} ({tree_name}): {missing_to_report}. Added as NaN.", file=sys.stderr)

            if not branches_to_fetch:
                print(f"❌ Error: No specified variables to read are available in tree '{tree_name}'.", file=sys.stderr)
                return pd.DataFrame()

            df = tree.arrays(branches_to_fetch, library="pd")

            # Ensure all requested variables are columns
            for var in local_variables:
                if var not in df.columns:
                    df[var] = np.nan

            # Basic type conversion for essential columns
            for col in ['runID', 'dy', 'spillID', 'mass']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # If ReWeight was explicitly requested and read, convert it
            if 'ReWeight' in variables_to_read and 'ReWeight' in df.columns:
                 df['ReWeight'] = pd.to_numeric(df['ReWeight'], errors='coerce').fillna(1.0)
            elif 'ReWeight' in variables_to_read and 'ReWeight' not in df.columns:
                df['ReWeight'] = 1.0 # Add default if requested but missing

            return df
    except Exception as e:
        print(f"❌ Error reading tree '{tree_name}' from '{file_path}': {e}", file=sys.stderr)
        return pd.DataFrame()

def add_beam_offset(df):
    if not df.empty: df['beamOffset'] = 1.6
    return df

def get_dataframe_after_cuts(df_initial, cuts_to_apply_ordered, good_spills_set_param, dataset_label_param, target_labels_good_spills_param):
    # (Function definition is unchanged from previous version)
    if df_initial.empty: return df_initial
    df = df_initial.copy()
    if 'beamOffset' not in df.columns: df['beamOffset'] = 1.6
    if 'spillID' in df.columns: df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')
    else: 
        if "Good Spills Cut" in cuts_to_apply_ordered: df['spillID'] = np.nan

    cumulative_mask = pd.Series(True, index=df.index)
    for cut_name, cut_string in cuts_to_apply_ordered.items():
        if not cumulative_mask.any(): break
        current_df_for_cut = df[cumulative_mask]
        previous_cumulative_mask_state = cumulative_mask.copy()
        try:
            mask_this_step = pd.Series(True, index=current_df_for_cut.index)
            if not current_df_for_cut.empty:
                if cut_name == "Good Spills Cut":
                    if dataset_label_param in target_labels_good_spills_param:
                        if good_spills_set_param is not None and 'spillID' in current_df_for_cut.columns and bool(good_spills_set_param):
                            current_spillID_numeric = pd.to_numeric(current_df_for_cut['spillID'], errors='coerce')
                            mask_this_step = current_spillID_numeric.isin(good_spills_set_param) & current_spillID_numeric.notna()
                else:
                    mask_this_step = current_df_for_cut.eval(cut_string, engine='python')
            aligned_mask_for_this_step = pd.Series(False, index=df.index)
            if not current_df_for_cut.empty:
                aligned_mask_for_this_step.loc[current_df_for_cut.index[mask_this_step]] = True
            cumulative_mask = cumulative_mask & aligned_mask_for_this_step
        except Exception as e:
            print(f"⚠️ Error applying cut '{cut_name}' for '{dataset_label_param}': {e}. Skipping this cut.", file=sys.stderr)
            cumulative_mask = previous_cumulative_mask_state
    return df[cumulative_mask]

# --- Worker function for parallel processing ---
def process_single_dataset_for_mass(label, file_path, tree_name, variables_list_arg,
                                    cuts_to_apply_arg, good_spills_arg,
                                    target_labels_good_spills_arg):
    """
    Processes a single dataset: reads ROOT file, applies cuts, extracts mass data.
    Designed to be run in a separate process.
    """
    # This function needs access to the helper functions defined above it.
    try:
        df_initial = read_tree(file_path, tree_name, variables_list_arg)
        if df_initial.empty:
            return {'label': label, 'mass_data': pd.Series(dtype='float64'), 'error': False, 'message': 'Initial DataFrame empty'}
            
        df_with_offset = add_beam_offset(df_initial)
        
        for var in variables_list_arg: # Ensure columns exist before applying cuts
            if var not in df_with_offset.columns:
                df_with_offset[var] = np.nan
        
        df_filtered = get_dataframe_after_cuts(df_with_offset, cuts_to_apply_arg, 
                                               good_spills_set_param=good_spills_arg,
                                               dataset_label_param=label,
                                               target_labels_good_spills_param=target_labels_good_spills_arg)
        
        mass_series = pd.Series(dtype='float64')
        if df_filtered.empty: pass
        elif 'mass' not in df_filtered.columns: pass
        else: mass_series = df_filtered['mass'].dropna()
        
        return {'label': label, 'mass_data': mass_series, 'error': False, 'message': ''}

    except Exception as e:
        print(f"❌ Error processing dataset {label} in PID {os.getpid()}: {e}", file=sys.stderr, flush=True)
        return {'label': label, 'mass_data': pd.Series(dtype='float64'), 'error': True, 'message': str(e)}


# --- Main Plotting Logic ---
def generate_mass_plots_after_cut_parallel(datasets_config, full_cuts_ordered_dict, target_cut_name_local, 
                                           variables_list_local, good_spills_file_path_local, 
                                           target_labels_good_spills_local,
                                           output_filename_base="mass_dist_after_cut_parallel"): # Base name
    """
    Generates mass distribution plots by processing datasets in parallel.
    """
    start_time = time.time()
    print("Loading good spills list...")
    sys.stdout.flush()
    good_spills = load_good_spills(good_spills_file_path_local)
    
    cuts_to_apply = OrderedDict()
    found_target_cut = False
    for name, condition in full_cuts_ordered_dict.items():
        cuts_to_apply[name] = condition
        if name == target_cut_name_local:
            found_target_cut = True
            break
    
    if not found_target_cut:
        print(f"❌ Error: Target cut '{target_cut_name_local}' not found. Cannot proceed.", file=sys.stderr)
        return

    print(f"\nWill apply cuts up to and including: '{target_cut_name_local}'")
    sys.stdout.flush()
    
    mass_data_collections = {}
    datasets_to_plot_labels = list(datasets_config.keys())
    
    futures = []
    print(f"Starting parallel processing with up to {os.cpu_count()} workers...")
    sys.stdout.flush()
    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        for label in datasets_to_plot_labels:
            config = datasets_config[label]
            futures.append(executor.submit(process_single_dataset_for_mass,
                                           label, config["path"], config["tree"],
                                           variables_list_local, cuts_to_apply, good_spills,
                                           target_labels_good_spills_local))

        processed_count = 0
        total_events_processed = {label: 0 for label in datasets_to_plot_labels} # Track events found
        for future in concurrent.futures.as_completed(futures):
            processed_count +=1
            try:
                result = future.result()
                current_label = result['label']
                mass_data_collections[current_label] = result['mass_data']
                total_events_processed[current_label] = len(result['mass_data'])
                if result['error']:
                    print(f"Error reported for dataset {current_label}: {result['message']}", file=sys.stderr)
                print(f"Processed {current_label} ({processed_count}/{len(futures)}). Events after cuts: {total_events_processed[current_label]}")
                sys.stdout.flush()
            except Exception as e:
                print(f"❌ Exception processing a future result: {e}", file=sys.stderr)
                sys.stdout.flush()

    processing_time = time.time() - start_time
    print(f"\nData processing finished in {processing_time:.2f} seconds.")

    if not mass_data_collections or all(v is None or v.empty for v in mass_data_collections.values()):
        print("No mass data collected from any dataset after cuts. Exiting plot generation.", file=sys.stderr)
        return

    print("Generating mass distribution plot...")
    sys.stdout.flush()
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(12, 7))
    
    # Updated plotting parameters
    mass_range_plot = (4, 10)
    bin_width = 0.1
    bins_plot = int(np.ceil((mass_range_plot[1] - mass_range_plot[0]) / bin_width)) # Use ceil for safety

    for label, mass_values in mass_data_collections.items():
        num_entries = total_events_processed.get(label, 0) # Get count from earlier processing
        if mass_values is not None and not mass_values.empty:
            plt.hist(mass_values, bins=bins_plot, range=mass_range_plot, 
                     histtype='step', linewidth=1.5, label=f"{label} (Entries: {num_entries})")
        else:
            plt.hist([], bins=bins_plot, range=mass_range_plot, label=f"{label} (Entries: 0)")

    plt.xlabel("Dimuon Mass (GeV/$c^2$)", fontsize=12)
    plt.ylabel(f"Events / ({bin_width} GeV/$c^2$)", fontsize=12)
    plt.title(f"Dimuon Mass Distribution after cut: '{target_cut_name_local}'", fontsize=14)
    plt.xlim(mass_range_plot) # Set x-axis limits explicitly
    
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    
    # Optional Log scale - check if data exists before setting limits
    if any(n > 0 for n in total_events_processed.values()):
        plt.yscale('log')
        current_ylim = plt.ylim()
        # Set lower limit reasonably for log scale (e.g., 0.5) if data exists
        plt.ylim(bottom=0.5, top=current_ylim[1] * 1.5 if current_ylim[1] > 1 else 100) # Adjust top limit if needed
    else: # No data, keep linear scale or set arbitrary log scale
        plt.ylim(bottom=0, top=10) # Linear scale if no data


    plt.tight_layout()
    
    output_filename_pdf = output_filename_base + ".pdf" # Save as PDF
    try:
        plt.savefig(output_filename_pdf, format='pdf', dpi=150) # Specify format
        print(f"Plot successfully saved as '{output_filename_pdf}'")
    except Exception as e:
        print(f"Error saving plot: {e}", file=sys.stderr)
    
    # plt.show() # Uncomment to display plot interactively

if __name__ == '__main__':
    print("Starting parallel mass distribution plotting script...")
    sys.stdout.flush()
    
    datasets_to_process_and_plot = {
        "Data(RS67)": FILE_PATHS_CONFIG["Data(RS67)"],
        "Mixed(RS67)": FILE_PATHS_CONFIG["Mixed(RS67)"],
        "empty flask (RS67)": FILE_PATHS_CONFIG["empty flask (RS67)"],
        "empty flask (RS67) mixed": FILE_PATHS_CONFIG["empty flask (RS67) mixed"],
    }
    
    # Base name for the output file, extension (.pdf) is added later
    output_plot_name_base = "mass_dist_after_Dcut_parallel_4to10GeV" 
    
    generate_mass_plots_after_cut_parallel(
        datasets_to_process_and_plot, 
        cuts_dict_full, 
        TARGET_CUT_NAME,
        variables_list,
        GOOD_SPILLS_FILE_PATH,
        TARGET_LABELS_FOR_GOOD_SPILLS_CUT_CONFIG,
        output_filename_base=output_plot_name_base
    )
    print("Script finished.")