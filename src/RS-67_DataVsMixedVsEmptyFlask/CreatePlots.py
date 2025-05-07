# -*- coding: utf-8 -*-
"""
Plot mass distribution for specific datasets after a target cut,
by re-processing ROOT files in parallel.

This script is designed to:
1. Load a predefined list of "good spills" from a text file.
2. Read data from specified ROOT files for several datasets.
3. Apply a sequence of event selection cuts up to a defined target cut.
4. Extract the 'mass' variable for events passing these cuts for each dataset.
5. Perform the data processing for each dataset in parallel to utilize multiple CPU cores.
6. Generate and save a plot comparing the mass distributions of these datasets.

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
    # "ReWeight" is intentionally omitted to avoid auto-reading unless explicitly needed by a future MC analysis.
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
    # Ensure essential analysis variables are always part of the request list
    local_variables = list(set(variables_to_read + ['mass', 'spillID', 'runID', 'dy']))
    try:
        with uproot.open(file_path) as file:
            if tree_name not in file:
                print(f"❌ Error: Tree '{tree_name}' not found in '{file_path}'.", file=sys.stderr)
                return pd.DataFrame()
            tree = file[tree_name]
            available_branches = tree.keys()
            
            # Determine which of the requested variables are actually available in the TTree
            branches_to_fetch = [var for var in local_variables if var in available_branches]
            missing_in_tree = [var for var in local_variables if var not in available_branches]
            
            if missing_in_tree:
                # Report only truly missing variables (excluding ReWeight if not explicitly asked for)
                missing_to_report = [b for b in missing_in_tree if b != 'ReWeight' or 'ReWeight' in variables_to_read]
                if missing_to_report:
                     print(f"ℹ️ Note: Branches not found in {file_path} ({tree_name}): {missing_to_report}. Will be added as NaN columns if not already present.", file=sys.stderr)

            if not branches_to_fetch: # No relevant branches found to read
                print(f"❌ Error: None of the primary variables to read are available in tree '{tree_name}' of file '{file_path}'.", file=sys.stderr)
                # Still return a DataFrame with expected columns for consistency downstream, filled with NaN
                df = pd.DataFrame()
                for var in local_variables: df[var] = np.nan
                return df
                
            df = tree.arrays(branches_to_fetch, library="pd")

            # Ensure all initially requested local_variables are columns, even if not in TTree (filled with NaN)
            for var in local_variables:
                if var not in df.columns:
                    df[var] = np.nan
            
            # Convert essential and plotting variables to numeric types
            for col in ['runID', 'dy', 'spillID', 'mass']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Handle ReWeight only if it was explicitly in variables_to_read and thus in local_variables
            if 'ReWeight' in variables_to_read:
                if 'ReWeight' in df.columns:
                    df['ReWeight'] = pd.to_numeric(df['ReWeight'], errors='coerce').fillna(1.0)
                else: # If requested but missing from tree (and thus branches_to_fetch)
                    df['ReWeight'] = 1.0 # Add as a column of 1.0s

            return df
    except Exception as e:
        print(f"❌ Error reading tree '{tree_name}' from '{file_path}': {e}", file=sys.stderr)
        # Return an empty DataFrame with expected columns on error for robust downstream handling
        df = pd.DataFrame()
        for var in local_variables: df[var] = np.nan # Ensure schema
        return df


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
    if df_initial.empty:
        return df_initial # Return empty if input is empty

    df = df_initial.copy() # Work on a copy to avoid modifying the original DataFrame passed in

    # Ensure essential columns for cuts are present and correctly typed
    if 'beamOffset' not in df.columns:
        # This should ideally be added before this function, but as a fallback:
        print("Warning: 'beamOffset' column not found in df for get_dataframe_after_cuts. Adding default.", file=sys.stderr)
        df['beamOffset'] = 1.6
    
    if "Good Spills Cut" in cuts_to_apply_ordered:
        if 'spillID' in df.columns:
            df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')
        else: # Add spillID as NaN if missing and Good Spills Cut is to be applied
            print("Warning: 'spillID' column missing, adding as NaN for Good Spills Cut.", file=sys.stderr)
            df['spillID'] = np.nan

    cumulative_mask = pd.Series(True, index=df.index)

    for cut_name, cut_string in cuts_to_apply_ordered.items():
        if not cumulative_mask.any():  # Optimization: if no events are left, stop applying cuts
            # print(f"ℹ️ No events remaining before applying cut '{cut_name}' for {dataset_label_param}. Further cuts skipped.", file=sys.stderr)
            break
        
        # Preserve the mask state before attempting the current cut
        previous_cumulative_mask_state = cumulative_mask.copy()
        
        try:
            # DataFrame view for the current cut, based on events that passed previous cuts
            current_df_for_cut = df[previous_cumulative_mask_state] 
            
            # By default, assume this cut step passes all events it sees
            mask_this_step_on_current_df = pd.Series(True, index=current_df_for_cut.index)

            if not current_df_for_cut.empty: # Only evaluate if there's data
                if cut_name == "Good Spills Cut":
                    if dataset_label_param in target_labels_good_spills_param: # Check if cut applies to this dataset
                        if good_spills_set_param is not None and 'spillID' in current_df_for_cut.columns and bool(good_spills_set_param):
                            current_spillID_numeric = pd.to_numeric(current_df_for_cut['spillID'], errors='coerce')
                            mask_this_step_on_current_df = current_spillID_numeric.isin(good_spills_set_param) & current_spillID_numeric.notna()
                        else: # Conditions for applying Good Spills Cut not met (e.g., empty spill list)
                            if not bool(good_spills_set_param): print(f"ℹ️ 'Good Spills Cut' for '{dataset_label_param}' ineffective (spill list empty/not loaded). Passing.", file=sys.stderr)
                            if 'spillID' not in current_df_for_cut.columns: print(f"⚠️ 'Good Spills Cut' for '{dataset_label_param}' skipped ('spillID' missing). Passing.", file=sys.stderr)
                    # If not a target dataset, mask_this_step_on_current_df remains True (passes all)
                else: # Regular cut string evaluation
                    mask_this_step_on_current_df = current_df_for_cut.eval(cut_string, engine='python')
            
            # Align the mask from this step (which is on current_df_for_cut's index) 
            # back to the original DataFrame's index (df.index) before combining.
            # Initialize with False on the original index.
            aligned_mask_for_this_step_on_original_index = pd.Series(False, index=df.index) 
            if not current_df_for_cut.empty: # Only try to align if there was data to filter
                aligned_mask_for_this_step_on_original_index.loc[current_df_for_cut.index[mask_this_step_on_current_df]] = True
            
            # Update the main cumulative_mask
            cumulative_mask = previous_cumulative_mask_state & aligned_mask_for_this_step_on_original_index

        except Exception as e:
            print(f"⚠️ Error applying cut '{cut_name}' for dataset '{dataset_label_param}': {e}. This specific cut will be skipped (no filtering effect for this step).", file=sys.stderr)
            # Revert cumulative_mask to its state before this failed cut attempt.
            # This means the failed cut has no filtering effect.
            cumulative_mask = previous_cumulative_mask_state 
            
    return df[cumulative_mask]


# --- Worker function for parallel processing ---
def process_single_dataset_for_mass(label, file_path, tree_name, variables_list_arg,
                                    cuts_to_apply_arg, good_spills_arg,
                                    target_labels_good_spills_arg):
    """Processes a single dataset to extract mass data after applying specified cuts.

    This function is designed to be executed in a separate process by
    `concurrent.futures.ProcessPoolExecutor`. It reads data from a ROOT file,
    applies beam offset corrections, filters events based on a series of cuts,
    and then extracts the 'mass' variable from the surviving events.

    :param label: A string label identifying the dataset (e.g., "Data(RS67)").
    :type label: str
    :param file_path: The path to the ROOT file for this dataset.
    :type file_path: str
    :param tree_name: The name of the TTree within the ROOT file.
    :type tree_name: str
    :param variables_list_arg: A list of variable names to be read from the TTree.
                               Must include 'mass' and all variables required by the cuts.
    :type variables_list_arg: list[str]
    :param cuts_to_apply_arg: An ordered dictionary of cuts to be applied sequentially.
                              Keys are cut names, values are pandas query strings.
    :type cuts_to_apply_arg: collections.OrderedDict[str, str]
    :param good_spills_arg: A set of good spill IDs for the "Good Spills Cut".
    :type good_spills_arg: set[int]
    :param target_labels_good_spills_arg: A set of dataset labels for which the
                                          "Good Spills Cut" is active.
    :type target_labels_good_spills_arg: set[str]
    :return: A dictionary containing the dataset label, a pandas Series of mass data,
             an error flag (bool), and an error message (str).
             Example: ``{'label': 'Data(RS67)', 'mass_data': pd.Series([...]), 'error': False, 'message': ''}``
    :rtype: dict
    """
    try:
        # print(f"--- [Process PID: {os.getpid()}] Starting processing for: {label} ---", flush=True) # Optional debug
        
        df_initial = read_tree(file_path, tree_name, variables_list_arg)
        
        if df_initial.empty:
            # print(f"Info: DataFrame for {label} is empty after reading. No mass data.", file=sys.stderr, flush=True)
            return {'label': label, 'mass_data': pd.Series(dtype='float64'), 'error': False, 'message': 'Initial DataFrame empty'}
            
        df_with_offset = add_beam_offset(df_initial)
        
        # Ensure all variables needed by cuts are columns in the DataFrame.
        # This is a safeguard; read_tree should ideally create them as NaN if missing from file.
        for var in variables_list_arg:
            if var not in df_with_offset.columns:
                df_with_offset[var] = np.nan
        
        df_filtered = get_dataframe_after_cuts(df_with_offset, cuts_to_apply_arg, 
                                               good_spills_set_param=good_spills_arg,
                                               dataset_label_param=label,
                                               target_labels_good_spills_param=target_labels_good_spills_arg)
        
        mass_series = pd.Series(dtype='float64') # Default to empty series
        if not df_filtered.empty and 'mass' in df_filtered.columns:
            mass_series = df_filtered['mass'].dropna()
        # elif df_filtered.empty:
            # print(f"Info: DataFrame for {label} is empty after applying cuts. No mass data.", file=sys.stderr, flush=True)
        # elif 'mass' not in df_filtered.columns:
            # print(f"Warning: 'mass' column not found in filtered DataFrame for {label}.", file=sys.stderr, flush=True)
        
        # print(f"--- [Process PID: {os.getpid()}] Finished processing for: {label}. Found {len(mass_series)} events with mass data.", flush=True)
        return {'label': label, 'mass_data': mass_series, 'error': False, 'message': ''}

    except Exception as e:
        # Capture any unexpected errors during the processing of this dataset
        print(f"❌ Unhandled Error processing dataset {label} in PID {os.getpid()}: {e}", file=sys.stderr, flush=True)
        return {'label': label, 'mass_data': pd.Series(dtype='float64'), 'error': True, 'message': str(e)}


# --- Main Plotting Logic ---
def generate_mass_plots_after_cut_parallel(datasets_config, full_cuts_ordered_dict, target_cut_name_local, 
                                           variables_list_local, good_spills_file_path_local, 
                                           target_labels_good_spills_local,
                                           output_filename_base="mass_dist_after_cut_parallel"):
    """Orchestrates parallel processing of datasets and generates mass distribution plots.

    This function performs the following steps:
    1. Loads the set of good spill IDs.
    2. Determines the sequence of cuts to apply up to the specified `target_cut_name_local`.
    3. Submits processing tasks for each dataset to a `ProcessPoolExecutor` for parallel execution.
       Each task involves reading data, applying cuts, and extracting mass values.
    4. Collects the mass data from all processed datasets.
    5. Generates a histogram plot comparing the mass distributions.
    6. Saves the plot to a PDF file.

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
        print(f"❌ Error: Target cut '{target_cut_name_local}' not found in the cuts dictionary. Cannot proceed.", file=sys.stderr)
        return

    print(f"\nWill apply cuts up to and including: '{target_cut_name_local}'")
    sys.stdout.flush()
    
    mass_data_collections = {}
    datasets_to_plot_labels = list(datasets_config.keys())
    
    futures = []
    # Use max_workers=None to use os.cpu_count() worker processes.
    num_workers = os.cpu_count()
    print(f"Starting parallel processing with up to {num_workers} workers...")
    sys.stdout.flush()

    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        for label in datasets_to_plot_labels:
            config = datasets_config[label]
            futures.append(executor.submit(process_single_dataset_for_mass,
                                           label, config["path"], config["tree"],
                                           variables_list_local, cuts_to_apply, good_spills,
                                           target_labels_good_spills_local))

        processed_count = 0
        total_events_plotted = {label: 0 for label in datasets_to_plot_labels}
        for future in concurrent.futures.as_completed(futures):
            processed_count += 1
            try:
                result = future.result()
                current_label = result['label']
                mass_data_collections[current_label] = result['mass_data']
                total_events_plotted[current_label] = len(result['mass_data']) # Store count for legend
                
                if result['error']:
                    print(f"Error was reported while processing dataset {current_label}: {result['message']}", file=sys.stderr)
                
                print(f"Completed processing for {current_label} ({processed_count}/{len(futures)}). Events for mass plot: {total_events_plotted[current_label]}")
                sys.stdout.flush()

            except Exception as e:
                # This catches errors if future.result() itself raises an exception (e.g., worker process died)
                print(f"❌ Exception occurred while retrieving result for a processed task: {e}", file=sys.stderr)
                sys.stdout.flush()

    processing_time = time.time() - start_time_main
    print(f"\nData processing phase finished in {processing_time:.2f} seconds.")

    if not mass_data_collections or all(v is None or v.empty for v in mass_data_collections.values()):
        print("No mass data collected from any dataset after cuts. Plot will be empty or not generated.", file=sys.stderr)
        # Optionally, still generate an empty plot for consistency or just return
        # For now, let it proceed, plt.hist handles empty data by plotting nothing for that dataset.
    
    print("Generating mass distribution plot...")
    sys.stdout.flush()
    plt.style.use('seaborn-v0_8-whitegrid') # Using a seaborn style
    plt.figure(figsize=(12, 7))
    
    mass_range_plot = (4, 10)  # X-axis range: 4 to 10 GeV/c^2
    bin_width = 0.1            # Bin size: 0.1 GeV/c^2
    bins_plot = int(np.ceil((mass_range_plot[1] - mass_range_plot[0]) / bin_width))

    for label, mass_values in mass_data_collections.items():
        num_entries = total_events_plotted.get(label, 0) # Use count from processing phase
        if mass_values is not None and not mass_values.empty:
            plt.hist(mass_values, bins=bins_plot, range=mass_range_plot, 
                     histtype='step', linewidth=1.5, label=f"{label} (Entries: {num_entries})")
        else: # Plot empty histogram to keep label in legend if needed
            plt.hist([], bins=bins_plot, range=mass_range_plot, label=f"{label} (Entries: 0)")

    plt.xlabel("Dimuon Mass (GeV/$c^2$)", fontsize=12)
    plt.ylabel(f"Events / ({bin_width:.1f} GeV/$c^2$)", fontsize=12) # Format bin width in label
    plt.title(f"Dimuon Mass Distribution after cut: '{target_cut_name_local}'", fontsize=14)
    plt.xlim(mass_range_plot) # Explicitly set x-axis limits
    
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    
    # Handle y-scale (log or linear)
    # Only use log scale if there's substantial data to show.
    if any(n > 0 for n in total_events_plotted.values()):
        plt.yscale('log')
        current_ylim = plt.ylim()
        # Ensure bottom y-limit is reasonable for log scale (e.g., 0.5 or 1 if counts start low)
        min_val_for_log = 0.5
        if current_ylim[0] < min_val_for_log :
            plt.ylim(bottom=min_val_for_log, top=current_ylim[1] * 1.5 if current_ylim[1] > min_val_for_log else 100)
    else: # No data or all zeros, linear scale might be better or set fixed y-range
        plt.ylim(bottom=0, top=10) # Default for empty plot

    plt.tight_layout()
    
    output_filename_pdf = output_filename_base + ".pdf"
    try:
        plt.savefig(output_filename_pdf, format='pdf', dpi=150) # Save as PDF
        print(f"Plot successfully saved as '{output_filename_pdf}'")
    except Exception as e:
        print(f"Error saving plot: {e}", file=sys.stderr)
    
    # plt.show() # Uncomment to display plot interactively after saving

if __name__ == '__main__':
    print("Starting parallel mass distribution plotting script...")
    sys.stdout.flush() # Ensure print statements appear before multiprocessing may fork
    
    datasets_to_process_and_plot = {
        "Data(RS67)": FILE_PATHS_CONFIG["Data(RS67)"],
        "Mixed(RS67)": FILE_PATHS_CONFIG["Mixed(RS67)"],
        "empty flask (RS67)": FILE_PATHS_CONFIG["empty flask (RS67)"],
        "empty flask (RS67) mixed": FILE_PATHS_CONFIG["empty flask (RS67) mixed"],
    }
    
    output_plot_name_base = "mass_dist_after_Dcut_parallel_4to10GeV_bin01" 
    
    generate_mass_plots_after_cut_parallel(
        datasets_to_process_and_plot, 
        cuts_dict_full, 
        TARGET_CUT_NAME,
        variables_list, # Pass the module-level list
        GOOD_SPILLS_FILE_PATH,
        TARGET_LABELS_FOR_GOOD_SPILLS_CUT_CONFIG, # Pass the module-level config
        output_filename_base=output_plot_name_base
    )
    print("Script finished.")