# -*- coding: utf-8 -*-
"""
Plot mass distribution for specific datasets and their combination
after a target cut, by re-processing ROOT files in parallel.

This script is designed to:
1. Load a predefined list of "good spills" from a text file.
2. Read data from specified ROOT files for several datasets.
3. Apply a sequence of event selection cuts up to a defined target cut.
4. Extract the 'mass' variable for events passing these cuts for each dataset.
5. Perform the data processing for each dataset in parallel to utilize multiple CPU cores.
6. Generate and save a plot comparing the mass distributions of these datasets,
   including a calculated "Corrected Data (RS67)" histogram.
7. The y-axis of the plot is set to a logarithmic scale.

The script includes functions for loading data, applying beam offsets,
applying selection cuts, and orchestrating the parallel processing and plotting.
Key configurations like file paths, cut definitions, and variables to be read
are defined as module-level constants.
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

#: Set of dataset labels for which the "Good Spills Cut" should be actively applied.
TARGET_LABELS_FOR_GOOD_SPILLS_CUT_CONFIG = {"Data", "Data(RS67)", "empty flask (RS67)"}

#: Path to the text file containing the list of good spill IDs.
GOOD_SPILLS_FILE_PATH = "../../res/GoodSpills/GoodSpills.txt"

#: The name of the cut in `cuts_dict_full` after which the mass distribution will be plotted.
TARGET_CUT_NAME = "D1 + D2 + D3 < 1000"

#: Factor used in the calculation of the "Corrected Data (RS67)" histogram.
#: Formula: Data(RS67) - Mixed(RS67) - FACTOR * (empty flask (RS67) - empty flask (RS67) mixed)
COMBINED_HIST_FACTOR = 16.1122 / 3.66242

try:
    #: Base path to the directory containing ROOT files for mixed events and empty flask data.
    #: Adjust this path based on the script's location relative to the 'res' directory.
    base_path_mixed_events = "../../res/ROOTFiles/MixedEvents/"
except NameError:
    print("Warning: base_path_mixed_events not defined, using placeholder './'", file=sys.stderr)
    base_path_mixed_events = "./"

#: Configuration dictionary mapping dataset labels to their ROOT file paths and TTree names.
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

#: Full ordered dictionary of all event selection cuts.
#: Keys are descriptive cut names, and values are pandas query strings or special identifiers.
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

#: List of variables required from ROOT TTrees for applying all defined cuts and for plotting 'mass'.
variables_list = [
    "spillID", "mass",
    "nHits1", "nHits2", "nHits1St1", "nHits2St1",
    "chisq1", "chisq2", "chisq_dimuon", "chisq1_target", "chisq2_target",
    "dx", "dy", "dz", "dpx", "dpy", "dpz",
    "D1", "D2", "D3", "xF", "xT", "xB", "costh", "intensityP",
    "runID", "trackSeparation",
    "x1_t", "x2_t", "y1_t", "y2_t", "x1_d", "x2_d", "y1_d", "y2_d", "z1_v", "z2_v",
    "x1_st1", "x2_st1", "y1_st1", "y2_st1", "y1_st3", "y2_st3",
    "py1_st1", "py2_st1", "py1_st3", "py2_st3",
    "pz1_st1", "pz2_st1", "pz1_st3", "pz2_st3",
    "px1_st1", "px2_st1", "px1_st3", "px2_st3",
    # "ReWeight" is intentionally omitted from this list for data processing.
]

# --- Helper Functions ---
def load_good_spills(file_path):
    """Loads good spill IDs from a text file.

    Each line in the file is expected to contain one integer spillID.
    Lines that are empty or cannot be converted to int are ignored, and warnings are printed.

    :param file_path: Path to the text file containing good spill IDs.
    :type file_path: str
    :return: A set of integer spill IDs. Returns an empty set if the file
             is not found, cannot be read, or contains no valid integer spill IDs.
    :rtype: set[int]
    """
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
    """Reads specified variables from a ROOT TTree into a pandas DataFrame.

    Handles missing files, TTrees, and variables gracefully. Missing requested variables
    are added as columns filled with NaN. Essential variables like 'mass', 'spillID',
    'runID', 'dy' are always attempted to be read and converted to numeric types.
    'ReWeight' is not automatically requested unless present in `variables_to_read`.

    :param file_path: Path to the ROOT file.
    :type file_path: str
    :param tree_name: Name of the TTree within the ROOT file.
    :type tree_name: str
    :param variables_to_read: List of variable names (branches) to read.
    :type variables_to_read: list[str]
    :return: A pandas DataFrame containing the requested variables. Returns an empty
             DataFrame on critical errors (e.g., file not found, tree not found).
    :rtype: pd.DataFrame
    """
    local_variables = list(set(variables_to_read + ['mass', 'spillID', 'runID', 'dy']))
    try:
        with uproot.open(file_path) as file:
            if tree_name not in file:
                print(f"❌ Error: Tree '{tree_name}' not found in '{file_path}'.", file=sys.stderr)
                df_empty = pd.DataFrame()
                for var_name in local_variables: df_empty[var_name] = np.nan # Ensure schema
                return df_empty
            tree = file[tree_name]
            available_branches = tree.keys()
            branches_to_fetch = [var for var in local_variables if var in available_branches]
            
            missing_in_tree = [var for var in local_variables if var not in available_branches]
            if missing_in_tree:
                missing_to_report = [b for b in missing_in_tree if b != 'ReWeight' or 'ReWeight' in variables_to_read]
                if missing_to_report:
                     print(f"ℹ️ Note: Branches not found in {file_path} ({tree_name}): {missing_to_report}. Will be added as NaN.", file=sys.stderr)

            if not branches_to_fetch:
                print(f"❌ Error: None of the primary variables to read are available in tree '{tree_name}'. Creating empty DF.", file=sys.stderr)
                df = pd.DataFrame()
                for var_name in local_variables: df[var_name] = np.nan
                return df
                
            df = tree.arrays(branches_to_fetch, library="pd")

            for var in local_variables:
                if var not in df.columns:
                    df[var] = np.nan
            
            for col in ['runID', 'dy', 'spillID', 'mass']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            if 'ReWeight' in variables_to_read: # Only handle ReWeight if explicitly requested
                if 'ReWeight' in df.columns:
                    df['ReWeight'] = pd.to_numeric(df['ReWeight'], errors='coerce').fillna(1.0)
                else: # If requested but missing from tree
                    df['ReWeight'] = 1.0 

            return df
    except Exception as e:
        print(f"❌ Error reading tree '{tree_name}' from '{file_path}': {e}", file=sys.stderr)
        df_err = pd.DataFrame()
        for var_name in local_variables: df_err[var_name] = np.nan
        return df_err

def add_beam_offset(df):
    """Adds a 'beamOffset' column with a fixed value (1.6) to the DataFrame.

    Modifies the DataFrame in place if not empty.

    :param df: Input pandas DataFrame.
    :type df: pd.DataFrame
    :return: The DataFrame with the 'beamOffset' column added or updated.
    :rtype: pd.DataFrame
    """
    if not df.empty:
        df['beamOffset'] = 1.6
    return df

def get_dataframe_after_cuts(df_initial, cuts_to_apply_ordered, good_spills_set_param, 
                             dataset_label_param, target_labels_good_spills_param):
    """Applies a dictionary of sequential cuts to a DataFrame and returns the filtered DataFrame.

    Handles a special "Good Spills Cut" conditionally based on `dataset_label_param`.
    If any cut string evaluation fails, that specific cut is skipped for the current dataset,
    and processing continues with data that passed previous cuts.

    :param df_initial: The initial DataFrame to which cuts will be applied.
    :type df_initial: pd.DataFrame
    :param cuts_to_apply_ordered: An ordered dictionary where keys are cut names (str)
                                  and values are pandas query strings or special identifiers.
    :type cuts_to_apply_ordered: collections.OrderedDict[str, str]
    :param good_spills_set_param: A set of good spill IDs for the "Good Spills Cut".
                                  Can be None if the cut is not used or spills are not loaded.
    :type good_spills_set_param: set[int] | None
    :param dataset_label_param: The label of the current dataset being processed. Used to
                                determine if conditional cuts like "Good Spills Cut" apply.
    :type dataset_label_param: str
    :param target_labels_good_spills_param: A set of dataset labels for which the
                                            "Good Spills Cut" should be actively applied.
    :type target_labels_good_spills_param: set[str]
    :return: A pandas DataFrame containing only the rows that passed all applied cuts.
             Returns an empty DataFrame if the initial DataFrame is empty or all events are cut.
    :rtype: pd.DataFrame
    """
    if df_initial.empty: return df_initial
    df = df_initial.copy()
    if 'beamOffset' not in df.columns: df['beamOffset'] = 1.6
    if 'spillID' in df.columns: df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')
    else: 
        if "Good Spills Cut" in cuts_to_apply_ordered: df['spillID'] = np.nan

    cumulative_mask = pd.Series(True, index=df.index)
    for cut_name, cut_string in cuts_to_apply_ordered.items():
        if not cumulative_mask.any(): break
        current_df_for_cut = df[cumulative_mask]; previous_cumulative_mask_state = cumulative_mask.copy()
        try:
            mask_this_step = pd.Series(True, index=current_df_for_cut.index)
            if not current_df_for_cut.empty:
                if cut_name == "Good Spills Cut":
                    if dataset_label_param in target_labels_good_spills_param:
                        if good_spills_set_param is not None and 'spillID' in current_df_for_cut.columns and bool(good_spills_set_param):
                            current_spillID_numeric = pd.to_numeric(current_df_for_cut['spillID'], errors='coerce')
                            mask_this_step = current_spillID_numeric.isin(good_spills_set_param) & current_spillID_numeric.notna()
                else: mask_this_step = current_df_for_cut.eval(cut_string, engine='python')
            aligned_mask_for_this_step = pd.Series(False, index=df.index)
            if not current_df_for_cut.empty: aligned_mask_for_this_step.loc[current_df_for_cut.index[mask_this_step]] = True
            cumulative_mask = cumulative_mask & aligned_mask_for_this_step
        except Exception as e: print(f"⚠️ Error applying cut '{cut_name}' for '{dataset_label_param}': {e}. Skipping cut.", file=sys.stderr); cumulative_mask = previous_cumulative_mask_state
    return df[cumulative_mask]

# --- Worker function for parallel processing ---
def process_single_dataset_for_mass(label, file_path, tree_name, variables_list_arg,
                                    cuts_to_apply_arg, good_spills_arg,
                                    target_labels_good_spills_arg):
    """Processes a single dataset: reads ROOT file, applies cuts, extracts mass data.

    This function is designed to be executed in a separate process by
    `concurrent.futures.ProcessPoolExecutor`.

    :param label: A string label identifying the dataset (e.g., "Data(RS67)").
    :type label: str
    :param file_path: The path to the ROOT file for this dataset.
    :type file_path: str
    :param tree_name: The name of the TTree within the ROOT file.
    :type tree_name: str
    :param variables_list_arg: A list of variable names to be read from the TTree.
    :type variables_list_arg: list[str]
    :param cuts_to_apply_arg: An ordered dictionary of cuts to be applied sequentially.
    :type cuts_to_apply_arg: collections.OrderedDict[str, str]
    :param good_spills_arg: A set of good spill IDs for the "Good Spills Cut".
    :type good_spills_arg: set[int]
    :param target_labels_good_spills_arg: A set of dataset labels for which the
                                          "Good Spills Cut" is active.
    :type target_labels_good_spills_arg: set[str]
    :return: A dictionary containing the dataset label, a pandas Series of mass data,
             an error flag (bool), and an error message (str).
    :rtype: dict
    """
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
        if not df_filtered.empty and 'mass' in df_filtered.columns:
            mass_series = df_filtered['mass'].dropna()
        
        return {'label': label, 'mass_data': mass_series, 'error': False, 'message': ''}

    except Exception as e:
        print(f"❌ Error processing dataset {label} in PID {os.getpid()}: {e}", file=sys.stderr, flush=True)
        return {'label': label, 'mass_data': pd.Series(dtype='float64'), 'error': True, 'message': str(e)}


# --- Main Plotting Logic ---
def generate_mass_plots_after_cut_parallel(datasets_config, full_cuts_ordered_dict, target_cut_name_local, 
                                           variables_list_local, good_spills_file_path_local, 
                                           target_labels_good_spills_local,
                                           output_filename_base="mass_dist_after_cut_parallel"):
    """Orchestrates parallel processing of datasets and generates mass distribution plots.

    Includes a combined "Corrected Data (RS67)" histogram. The y-axis is set to a log scale.

    :param datasets_config: Dictionary mapping dataset labels to their file path and TTree name.
    :type datasets_config: dict
    :param full_cuts_ordered_dict: OrderedDict of all possible cuts.
    :type full_cuts_ordered_dict: collections.OrderedDict[str, str]
    :param target_cut_name_local: Name of the cut after which mass distributions are plotted.
    :type target_cut_name_local: str
    :param variables_list_local: List of variables to read from ROOT files.
    :type variables_list_local: list[str]
    :param good_spills_file_path_local: Path to the file containing good spill IDs.
    :type good_spills_file_path_local: str
    :param target_labels_good_spills_local: Set of dataset labels for which "Good Spills Cut" is active.
    :type target_labels_good_spills_local: set[str]
    :param output_filename_base: Base name for the output PDF plot file (without extension).
    :type output_filename_base: str
    :return: None. The function generates and saves a plot.
    :rtype: None
    """
    start_time_main = time.time()
    print("Loading good spills list..."); sys.stdout.flush()
    good_spills = load_good_spills(good_spills_file_path_local)
    
    cuts_to_apply = OrderedDict()
    found_target_cut = False
    for name, condition in full_cuts_ordered_dict.items():
        cuts_to_apply[name] = condition
        if name == target_cut_name_local: found_target_cut = True; break
    if not found_target_cut: print(f"❌ Error: Target cut '{target_cut_name_local}' not found.", file=sys.stderr); return

    print(f"\nWill apply cuts up to and including: '{target_cut_name_local}'"); sys.stdout.flush()
    
    mass_data_collections = {}; datasets_to_plot_labels = list(datasets_config.keys()); futures = []
    num_workers = os.cpu_count(); print(f"Starting parallel processing with up to {num_workers} workers..."); sys.stdout.flush()
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        for label in datasets_to_plot_labels:
            config = datasets_config[label]
            futures.append(executor.submit(process_single_dataset_for_mass, label, config["path"], config["tree"], variables_list_local, cuts_to_apply, good_spills, target_labels_good_spills_local))
        
        processed_count = 0; total_events_processed = {label: 0 for label in datasets_to_plot_labels}
        for future in concurrent.futures.as_completed(futures):
            processed_count += 1
            try:
                result = future.result(); current_label = result['label']
                mass_data_collections[current_label] = result['mass_data']
                total_events_processed[current_label] = len(result['mass_data'])
                if result['error']: print(f"Error reported for {current_label}: {result['message']}", file=sys.stderr)
                print(f"Processed {current_label} ({processed_count}/{len(futures)}). Events for mass plot: {total_events_processed[current_label]}"); sys.stdout.flush()
            except Exception as e: print(f"❌ Exception processing future result: {e}", file=sys.stderr); sys.stdout.flush()
    
    processing_time = time.time() - start_time_main; print(f"\nData processing finished in {processing_time:.2f} seconds.")
    if not mass_data_collections or all(v is None or v.empty for v in mass_data_collections.values()):
        print("No mass data collected. Exiting plot generation.", file=sys.stderr); return
    
    print("Generating mass distribution plot..."); sys.stdout.flush()
    plt.style.use('seaborn-v0_8-whitegrid'); plt.figure(figsize=(12, 7))
    
    mass_range_plot = (4, 10); bin_width = 0.1
    bins_plot = int(np.ceil((mass_range_plot[1] - mass_range_plot[0]) / bin_width))
    
    hist_counts_collection = {}; bin_edges_common = None
    datasets_for_combination = ["Data(RS67)", "Mixed(RS67)", "empty flask (RS67)", "empty flask (RS67) mixed"]
    for label in datasets_for_combination:
        mass_values = mass_data_collections.get(label)
        if mass_values is not None and not mass_values.empty:
            counts, edges = np.histogram(mass_values, bins=bins_plot, range=mass_range_plot)
            hist_counts_collection[label] = counts
            if bin_edges_common is None: bin_edges_common = edges
        else:
            hist_counts_collection[label] = np.zeros(bins_plot)
            if bin_edges_common is None: _, bin_edges_common = np.histogram([], bins=bins_plot, range=mass_range_plot)

    data_rs67_h = hist_counts_collection.get("Data(RS67)", np.zeros(bins_plot))
    mixed_rs67_h = hist_counts_collection.get("Mixed(RS67)", np.zeros(bins_plot))
    ef_rs67_h = hist_counts_collection.get("empty flask (RS67)", np.zeros(bins_plot))
    ef_mixed_h = hist_counts_collection.get("empty flask (RS67) mixed", np.zeros(bins_plot))
    
    combined_counts = data_rs67_h - mixed_rs67_h - COMBINED_HIST_FACTOR * (ef_rs67_h - ef_mixed_h)
    
    for label in datasets_to_plot_labels:
        mass_values = mass_data_collections.get(label)
        num_entries = total_events_processed.get(label, 0)
        if mass_values is not None and not mass_values.empty:
            plt.hist(mass_values, bins=bins_plot, range=mass_range_plot, histtype='step', linewidth=1.5, label=f"{label} (Entries: {num_entries})")
        else:
            plt.hist([], bins=bins_plot, range=mass_range_plot, histtype='step', label=f"{label} (Entries: 0)")

    if bin_edges_common is not None:
        plt.step(bin_edges_common[:-1], combined_counts, where='post', 
                 label="Corrected Data (RS67)", color='black', linewidth=2.0, linestyle='--') # Changed label

    plt.xlabel("Dimuon Mass (GeV/$c^2$)", fontsize=12)
    plt.ylabel(f"Events / ({bin_width:.1f} GeV/$c^2$)", fontsize=12)
    plt.title(f"Dimuon Mass Distribution after cut: '{target_cut_name_local}'", fontsize=14)
    plt.xlim(mass_range_plot); plt.xticks(fontsize=10); plt.yticks(fontsize=10)
    plt.legend(fontsize=9); plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    
    plt.yscale('log') # Set y-axis to log scale
    current_ylim = plt.ylim()
    # Adjust y-limits for log scale, ensuring bottom is positive
    bottom_limit = 0.1 # Default small positive number for log scale
    # Find minimum positive count across all histograms for a more dynamic lower limit
    all_positive_counts_for_ylim = []
    for label in datasets_to_plot_labels:
        if hist_counts_collection.get(label) is not None:
            all_positive_counts_for_ylim.extend(hist_counts_collection.get(label)[hist_counts_collection.get(label) > 0])
    if combined_counts[combined_counts > 0].size > 0: # Consider positive part of combined for y_lim
        all_positive_counts_for_ylim.extend(combined_counts[combined_counts > 0])

    if all_positive_counts_for_ylim:
        min_val = np.min(all_positive_counts_for_ylim)
        bottom_limit = max(0.05, min_val * 0.5) # Adjust based on data, but not too small
    
    plt.ylim(bottom=bottom_limit, top=current_ylim[1] * 1.5 if current_ylim[1] > bottom_limit else bottom_limit * 100)

    plt.tight_layout()
    output_filename_pdf = output_filename_base + ".pdf"
    try:
        plt.savefig(output_filename_pdf, format='pdf', dpi=150)
        print(f"Plot successfully saved as '{output_filename_pdf}'")
    except Exception as e: print(f"Error saving plot: {e}", file=sys.stderr)
    # plt.show()

if __name__ == '__main__':
    print("Starting parallel mass distribution plotting script with combined histogram (log scale)...")
    sys.stdout.flush()
    datasets_to_process_and_plot = {
        "Data(RS67)": FILE_PATHS_CONFIG["Data(RS67)"],
        "Mixed(RS67)": FILE_PATHS_CONFIG["Mixed(RS67)"],
        "empty flask (RS67)": FILE_PATHS_CONFIG["empty flask (RS67)"],
        "empty flask (RS67) mixed": FILE_PATHS_CONFIG["empty flask (RS67) mixed"],
    }
    output_plot_name_base = "mass_dist_corrected_after_Dcut_4to10GeV_logY" 
    generate_mass_plots_after_cut_parallel(
        datasets_to_process_and_plot, cuts_dict_full, TARGET_CUT_NAME,
        variables_list, GOOD_SPILLS_FILE_PATH,
        TARGET_LABELS_FOR_GOOD_SPILLS_CUT_CONFIG,
        output_filename_base=output_plot_name_base
    )
    print("Script finished.")