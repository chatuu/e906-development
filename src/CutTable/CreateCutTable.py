# -*- coding: utf-8 -*-
"""
Process High Energy Physics data from ROOT files to generate a cut flow table.

This script reads particle physics data (dimuon events) from specified ROOT filabout:blank#blockedes,
including experimental data, Monte Carlo (MC) simulations, and mixed event samples.
It applies a series of predefined selection cuts sequentially to filter events.

Key Features:
- Reads data from ROOT files using the 'uproot' library.
- Handles multiple datasets (Data, MC types, Mixed Events).
- Applies run-independent beam offset corrections (currently fixed).
- Processes files in parallel using 'concurrent.futures.ProcessPoolExecutor' for speed.
- Applies sequential event selection cuts defined in a dictionary.
- Calculates weighted sums for MC events using 'ReWeight' column.
- Aggregates results into a pandas DataFrame.
- Calculates derived columns like 'Total MC', 'Purity (DY MC)', 'Efficiency (DY MC)'.
- Prints a formatted cut flow table to the console using 'tabulate'.
- Saves the resulting cut flow table to a CSV file.
- Includes error handling and warnings during file reading and processing.
"""

import uproot
import pandas as pd
from tabulate import tabulate
import numpy as np
import concurrent.futures
import time # Import time to potentially measure speedup
import re # Import regex for improved variable extraction from cuts
import sys # Import sys to flush stdout for better interleaved printing

# Reuse the existing helper functions, ensuring they are self-contained for multiprocessing
# (i.e., they don't rely on shared global state that isn't passed as arguments).
# read_tree, add_beam_offset, and apply_cuts are suitable.

def read_tree(file_path, tree_name, variables):
    """
    Reads specified variables from a ROOT tree into a pandas DataFrame.

    Handles missing files, trees, and variables gracefully. If the file or
    tree cannot be read, an empty DataFrame with expected columns is returned.
    If variables are missing within a tree, they are added as columns
    filled with NaN.

    It ensures essential columns ('runID', 'dy') are attempted and converts
    them to numeric types if present. The 'ReWeight' column, crucial for MC,
    is also handled: if present, NaNs are filled with 1.0; if absent, it's
    added with a default value of 1.0. Warnings are printed for missing
    trees or variables.

    Parameters
    ----------
    file_path : str
        Path to the ROOT file.
    tree_name : str
        Name of the TTree within the ROOT file.
    variables : list[str]
        List of variable names (branches) to read from the TTree. 'runID', 'dy',
        and 'ReWeight' are always added to this list internally if not present.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the requested variables.
        - Columns specified in `variables` (if found).
        - 'runID', 'dy' columns (if found in the tree).
        - 'ReWeight' column (using tree value or default 1.0).
        - Missing requested variables added as columns with NaN values.
        Returns an empty DataFrame (with expected columns but no rows)
        if the file or tree is inaccessible or reading fails critically.

    Raises
    ------
    Exception
        Catches and prints generic exceptions during file opening or reading,
        returning an empty DataFrame instead of raising further.
    """
    try:
        with uproot.open(file_path) as file:
            if tree_name not in file:
                # Print error if tree is missing (Note: May interleave in parallel execution)
                print(f"❌ Error: Tree '{tree_name}' not found in file '{file_path}'.")
                sys.stdout.flush() # Try to flush output buffer
                # Return empty DataFrame with relevant columns for downstream processing structure
                empty_df = pd.DataFrame()
                # Include potentially required columns in empty df
                for var in set(variables + ['runID', 'dy', 'ReWeight']):
                     empty_df[var] = pd.Series(dtype='float64') # Use float to accommodate NaN
                return empty_df


            tree = file[tree_name]
            all_keys = tree.keys()

            # Ensure variables needed for beamOffset calculation and the cuts are requested
            # Use a set to avoid duplicates, then convert back to list
            # Keep runID as it might still be used in cuts_dict
            vars_to_request = list(set(variables + ['runID', 'dy', 'ReWeight']))

            vars_to_read = [var for var in vars_to_request if var in all_keys]
            missing = [var for var in vars_to_request if var not in all_keys]
            if missing:
                # Filter out 'ReWeight' from missing warning if it's likely data
                is_mc_file = any(mc_part in file_path for mc_part in ['mc_drellyan', 'mc_jpsi', 'mc_psiprime'])
                # Create a copy to modify for warning message
                missing_for_warning = list(missing)
                if 'ReWeight' in missing_for_warning and not is_mc_file:
                     missing_for_warning.remove('ReWeight')

                # Print warning for missing variables if any remain after filtering ReWeight
                # (Note: May interleave in parallel execution)
                if missing_for_warning:
                    print(f"⚠️ Missing variables in {file_path}: {missing_for_warning}. These will be added as NaN.")
                    sys.stdout.flush() # Try to flush output buffer


            df = tree.arrays(vars_to_read, library="pd")

            # Add missing columns with NaN values to ensure DataFrame has expected structure
            for var in missing:
                 if var not in df.columns: # Avoid adding if it somehow got added during read
                     df[var] = np.nan

            # Ensure runID, dy, and ReWeight are numeric
            if 'runID' in df.columns:
                df['runID'] = pd.to_numeric(df['runID'], errors='coerce')
            if 'dy' in df.columns:
                df['dy'] = pd.to_numeric(df['dy'], errors='coerce')
            if 'ReWeight' in df.columns:
                # Fill NaN weights with 1.0 (no reweighting)
                df['ReWeight'] = pd.to_numeric(df['ReWeight'], errors='coerce').fillna(1.0)
            else:
                # If ReWeight wasn't even requested/found, add it with 1.0 for consistency
                 df['ReWeight'] = 1.0


            return df

    except Exception as e:
        # Print error if reading fails (Note: May interleave in parallel execution)
        print(f"❌ Error reading tree '{tree_name}' from '{file_path}': {e}")
        sys.stdout.flush() # Try to flush output buffer
        # Return empty DataFrame on major error for structure consistency
        empty_df = pd.DataFrame()
        # Keep runID as it might still be used in cuts_dict
        for var in set(variables + ['runID', 'dy', 'ReWeight']):
             empty_df[var] = pd.Series(dtype='float64')
        return empty_df


def add_beam_offset(df):
    """
    Adds the 'beamOffset' column to the DataFrame with a fixed value.

    This simplified version sets a constant beam offset for all events,
    currently fixed at 1.6, regardless of the run ID or other event properties.
    It modifies the DataFrame in place by adding the new column.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame, potentially containing event data.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with an added 'beamOffset' column set to 1.6.
    """
    # Set the beamOffset to a fixed value of 1.6 for all rows
    df['beamOffset'] = 1.6
    return df


def apply_cuts(df, cuts, use_weights=True):
    """
    Applies a dictionary of sequential cuts to a pandas DataFrame.

    Iterates through the provided dictionary of cuts. For each cut, it
    evaluates the cut string expression on the *currently filtered* DataFrame.
    The results (number of events or sum of weights passing) are stored.
    Cuts are applied cumulatively; events failing one cut are excluded from
    evaluation in subsequent cuts.

    Handles potential errors during cut evaluation by printing a warning and
    assigning a result of 0 for the problematic cut, allowing subsequent
    cuts to proceed on the data that passed previous *successful* cuts.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing event data. Must include columns referenced
        in the `cuts` expressions, plus 'ReWeight' if `use_weights` is True.
    cuts : dict[str, str]
        A dictionary where keys are descriptive cut names (str) and values are
        string expressions (str) that can be evaluated by `pandas.eval()`.
        The order of cuts in the dictionary determines the sequence of application.
    use_weights : bool, optional
        If True (default), uses the 'ReWeight' column in `df` to calculate the
        sum of weights passing each cut. If False, counts the number of events
        (effectively weight=1).

    Returns
    -------
    pd.Series
        A pandas Series containing the results.
        - Index: Cut names from the `cuts` dictionary keys, plus "Total Events".
        - Values: Sum of weights (or event count if `use_weights=False`)
          remaining after applying the corresponding cut and all preceding cuts.
    """
    cut_results = {}
    # Use 'ReWeight' column if it exists and use_weights is True, otherwise use count of 1
    weight = df['ReWeight'] if use_weights and 'ReWeight' in df.columns else pd.Series(1.0, index=df.index)

    # Ensure weight is numeric
    weight = pd.to_numeric(weight, errors='coerce').fillna(1.0)

    cut_results["Total Events"] = weight.sum()

    # Ensure beamOffset column exists before cuts that might use it
    # This should be guaranteed by calling add_beam_offset before apply_cuts
    if 'beamOffset' not in df.columns:
         df['beamOffset'] = 1.6 # Default to 1.6 if somehow missing

    # Ensure dy is numeric before cuts that use it
    if 'dy' in df.columns:
         df['dy'] = pd.to_numeric(df['dy'], errors='coerce').fillna(np.nan)

    # Ensure runID is numeric before cuts that use it (still needed for some existing cuts)
    if 'runID' in df.columns:
         df['runID'] = pd.to_numeric(df['runID'], errors='coerce').fillna(-1) # Fill NaN runIDs for safe evaluation

    # Keep track of the mask from successful cuts
    cumulative_mask = pd.Series(True, index=df.index)


    # It's more robust to apply cuts to the *original* DataFrame using the cumulative mask
    original_df = df.copy()
    original_weight = weight.copy()


    for cut_name, cut_string in cuts.items():
        if cut_name == "Total Events":
             continue
        try:
            # Evaluate the cut string on the *current* filtered DataFrame 'df'
            mask = df.eval(cut_string, engine='python')

            # Reindex the mask to the original_df's index
            reindexed_mask = pd.Series(False, index=original_df.index)
            reindexed_mask[mask.index] = mask

            # Update the cumulative mask using the reindexed mask
            cumulative_mask = cumulative_mask & reindexed_mask

            # Calculate the sum using the cumulative mask applied to the original weights
            cut_results[cut_name] = original_weight[cumulative_mask].sum()

            # Filter the DataFrame 'df' for the *next* iteration.
            df = original_df[cumulative_mask].copy()

        except Exception as e:
            # Print warning if a cut fails (Note: May interleave in parallel execution)
            print(f"⚠️ Error applying cut '{cut_name}': {e}. Skipping cut for this dataset.")
            sys.stdout.flush() # Try to flush output buffer
            cut_results[cut_name] = 0 # Assign 0 count for the failed cut


    return pd.Series(cut_results)

def process_single_file(file_path, tree_name, label, variables, cuts, use_weights):
    """
    Processes a single ROOT file/tree: reads data, adds offset, applies cuts.

    This function serves as the target for parallel execution. It encapsulates
    the logic for handling one dataset specified by its file path and tree name.
    It calls `read_tree`, `add_beam_offset`, and `apply_cuts` sequentially.
    Robust error handling ensures that failure in one file doesn't stop others.

    Parameters
    ----------
    file_path : str
        Path to the input ROOT file.
    tree_name : str
        Name of the TTree within the ROOT file.
    label : str
        A descriptive label for this dataset (e.g., "Data", "DY MC"). Used
        as a key in the results dictionary and potentially as a column name.
    variables : list[str]
        List of variable names required from the TTree.
    cuts : dict[str, str]
        Dictionary defining the cuts to apply (see `apply_cuts`).
    use_weights : bool
        Flag passed to `apply_cuts` to indicate whether to use event weights.

    Returns
    -------
    dict
        A dictionary containing the processing results or error information:
        - 'label' (str): The input label for identification.
        - 'results' (pd.Series): The cut flow results from `apply_cuts`.
          Contains results with 0 counts if an error occurred during processing.
        - 'error' (bool): True if a significant error occurred during processing,
          False otherwise.
        - 'message' (str, optional): An error message if `error` is True.
    """
    try:
        # Add process start message here for better tracking in parallel output
        print(f"Starting processing for {label} from {file_path}...")
        sys.stdout.flush()
        # Make sure runID is requested if cuts use it, even if add_beam_offset doesn't
        vars_to_read = list(set(variables + ['runID']))
        df = read_tree(file_path, tree_name, vars_to_read)

        if df.empty and not ('runID' in df.columns and 'dy' in df.columns and 'ReWeight' in df.columns):
             print(f"INFO: read_tree returned empty DataFrame for {label}. Skipping further processing.")
             sys.stdout.flush()
             expected_cut_names = ["Total Events"] + list(cuts.keys())
             # Return structure indicating success but empty data
             return {'label': label, 'results': pd.Series(0.0, index=expected_cut_names), 'error': False}


        df = add_beam_offset(df) # This now sets beamOffset to 1.6

        # Ensure all variables needed for cuts plus beamOffset and dy exist before applying cuts
        cut_variables = set()
        for cut_string in cuts.values():
            found_vars = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', cut_string)
            valid_vars = [var for var in found_vars if var not in ['and', 'or', 'not', 'in', 'is', 'True', 'False', 'None', 'abs', 'np', 'sqrt', 'min', 'max']]
            cut_variables.update(valid_vars)

        # Combine requested variables with variables found in cut strings
        # Make sure runID is included if the cuts explicitly use it
        required_for_cuts = list(set(variables + list(cut_variables) + ['beamOffset', 'dy', 'runID']))


        for var in required_for_cuts:
             if var not in df.columns:
                  # Print warning if adding missing variable (Note: May interleave in parallel execution)
                  print(f"⚠️ Adding missing required cut variable '{var}' as NaN to {label} DataFrame.")
                  sys.stdout.flush()
                  df[var] = np.nan # Add any missing required variables


        cut_series = apply_cuts(df, cuts, use_weights=use_weights)
        print(f"Finished processing {label}.")
        sys.stdout.flush()
        return {'label': label, 'results': cut_series, 'error': False}

    except Exception as e:
        # Print error from the overall processing function (Note: May interleave in parallel execution)
        print(f"❌ Error processing {label} from {file_path}: {e}")
        sys.stdout.flush()
        expected_cut_names = ["Total Events"] + list(cuts.keys())
        # Return structure indicating error
        return {'label': label, 'results': pd.Series(0.0, index=expected_cut_names), 'error': True, 'message': str(e)}


def create_cut_table_parallel(data_file_path, mc_file_paths, mc_labels, tree_name, cuts, variables, mixed_file_path, mixed_tree_name):
    """
    Creates a cut flow table by processing multiple datasets in parallel.

    Orchestrates the parallel processing of data, Monte Carlo (MC), and
    mixed event samples using `concurrent.futures.ProcessPoolExecutor`.
    It submits each dataset to the `process_single_file` function.
    After all tasks complete, it aggregates the results into a single
    pandas DataFrame. Finally, it calculates summary columns like
    'Total MC', 'Purity (DY MC)', and 'Efficiency (DY MC)'.

    Parameters
    ----------
    data_file_path : str
        Path to the main data ROOT file.
    mc_file_paths : list[str]
        List of paths to the MC simulation ROOT files.
    mc_labels : list[str]
        List of labels corresponding to each MC file in `mc_file_paths`.
        Used for identifying results and as column names.
    tree_name : str
        Name of the TTree containing event data within the data and MC files.
    cuts : dict[str, str]
        Dictionary defining the cuts to apply (see `apply_cuts`).
    variables : list[str]
        List of essential variable names to read from the ROOT trees.
    mixed_file_path : str
        Path to the ROOT file containing mixed event samples (and potentially
        a corresponding data sample, e.g., RS67 Data).
    mixed_tree_name : str
        Name of the TTree for the mixed event sample within `mixed_file_path`.
        Assumes the corresponding data sample (if present) is in a tree named "result".

    Returns
    -------
    pd.DataFrame or None
        A pandas DataFrame representing the complete cut flow table.
        Rows correspond to cuts (index="Cut Name"), and columns correspond to
        the different datasets and calculated values (Data, MCs, Total MC, etc.).
        Returns None or an empty DataFrame if critical errors prevent table
        construction (e.g., no tasks complete successfully).

    Notes
    -----
    - Assumes the data sample corresponding to the mixed events is in a tree
      named "result" within the `mixed_file_path` file.
    - Parallel execution uses a number of worker processes based on available CPU cores.
    - Handles errors from individual tasks gracefully, printing warnings and
      filling results with zeros for failed tasks.
    """
    start_time = time.time()
    print("Starting parallel processing...")
    sys.stdout.flush()

    results_dict = {} # To store results Series keyed by label
    futures = [] # To store future objects

    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        # Submit Data task
        futures.append(executor.submit(process_single_file, data_file_path, tree_name, "Data", variables, cuts, False))

        # Submit Data(RS67) task
        futures.append(executor.submit(process_single_file, mixed_file_path, "result", "Data(RS67)", variables, cuts, False))

        # Submit Mixed(RS67) task
        futures.append(executor.submit(process_single_file, mixed_file_path, mixed_tree_name, "Mixed(RS67)", variables, cuts, False))

        # Submit MC tasks
        for i, mc_file_path in enumerate(mc_file_paths):
            label = mc_labels[i]
            mc_vars = list(set(variables + ["ReWeight"]))
            futures.append(executor.submit(process_single_file, mc_file_path, tree_name, label, mc_vars, cuts, True))

        # Collect results as they complete
        print("Waiting for tasks to complete...")
        sys.stdout.flush()
        processed_tasks = 0
        total_tasks = len(futures)
        for future in concurrent.futures.as_completed(futures):
            processed_tasks += 1
            print(f"Processing task {processed_tasks}/{total_tasks}...")
            sys.stdout.flush()
            try:
                result = future.result() # Get the result from the worker process
                label = result['label']
                cut_series = result['results']
                is_error = result['error']
                error_message = result.get('message', '')

                if is_error:
                    # Print error message returned from the worker process
                    print(f"❗ Task for '{label}' completed with error: {error_message}. Using zero counts.")
                    sys.stdout.flush()
                    # Ensure the Series index matches expected cuts even on error
                    expected_cut_names = ["Total Events"] + list(cuts.keys())
                    results_dict[label] = pd.Series(0.0, index=expected_cut_names)
                else:
                     results_dict[label] = cut_series
                     print(f"Successfully processed results for '{label}'.")
                     sys.stdout.flush()

            except Exception as exc:
                # Catch exceptions raised *during* future handling/result retrieval
                print(f'❌ Exception occurred while processing future result (label unknown): {exc}')
                sys.stdout.flush()
                # We don't know which label failed here, so we can't easily add a placeholder column.
                # This indicates a more fundamental problem potentially.


    # --- Aggregate results into a DataFrame ---
    print("Aggregating results into DataFrame...")
    sys.stdout.flush()
    cut_table = pd.DataFrame(results_dict)

    # Ensure all expected cut names are in the index
    expected_indices = ["Total Events"] + list(cuts.keys())
    cut_table = cut_table.reindex(expected_indices, fill_value=0.0)


    # --- Perform final sequential calculations (Total MC, Purity, Efficiency) ---
    print("Aggregating results and calculating final columns...")
    sys.stdout.flush()

    # Calculate Total MC
    mc_cols = [label for label in mc_labels if label in cut_table.columns]
    if mc_cols:
         cut_table["Total MC"] = cut_table[mc_cols].sum(axis=1, min_count=1).fillna(0)
    else:
         print("⚠️ No MC columns found in results to calculate 'Total MC'. Setting to 0.")
         sys.stdout.flush()
         cut_table["Total MC"] = pd.Series(0.0, index=cut_table.index)


    # Calculate Purity & Efficiency for DY MC
    try:
        dy_col = "DY MC"
        # Check if necessary columns and index exist before calculation
        if dy_col in cut_table.columns and "Total Events" in cut_table.index and "Total MC" in cut_table.columns:
             total_dy = cut_table.loc["Total Events", dy_col]

             # Check if Total DY events is zero
             if total_dy == 0:
                 print(f"⚠️ Total Events for '{dy_col}' is zero. Efficiency will be zero.")
                 sys.stdout.flush()

             total_mc_safe = cut_table["Total MC"].replace(0, np.nan) # Avoid division by zero
             purity = (cut_table[dy_col] / total_mc_safe).replace([np.inf, -np.inf], np.nan).fillna(0).copy()

             total_dy_safe = total_dy if total_dy != 0 else np.nan # Avoid division by zero
             efficiency = (cut_table[dy_col] / total_dy_safe).replace([np.inf, -np.inf], np.nan).fillna(0).copy()

             cut_table["Purity (DY MC)"] = purity * 100
             cut_table["Efficiency (DY MC)"] = efficiency * 100
        else:
            # Print warning if calculation cannot be performed
            print(f"⚠️ Cannot calculate Purity/Efficiency: Check if '{dy_col}' column, 'Total Events' index, and 'Total MC' column exist.")
            sys.stdout.flush()
            cut_table["Purity (DY MC)"] = 0.0
            cut_table["Efficiency (DY MC)"] = 0.0

    except Exception as e:
        # Print error if calculation fails unexpectedly
        print(f"❌ Error calculating Purity/Efficiency: {e}")
        sys.stdout.flush()
        cut_table["Purity (DY MC)"] = 0.0
        cut_table["Efficiency (DY MC)"] = 0.0


    cut_table.index.name = "Cut Name"

    end_time = time.time()
    print(f"\nParallel processing finished in {end_time - start_time:.2f} seconds.")
    sys.stdout.flush()

    return cut_table

# ------------------- Configuration -------------------
# Script configuration section: Defines cuts, variables, file paths, etc.

cuts_dict = {
    #"Trigger Level Cut": "matrix1 == 1",
    "z1_v and z2_v within -320cm to -5cm": "z1_v > -320 and z1_v < -5 and z2_v > -320 and z2_v < -5",
    "xt^2 + (yt - beamOffset)^2 < 320cm^2": "((x1_t**2 + (y1_t - beamOffset)**2 < 320)) and ((x2_t**2 + (y2_t - beamOffset)**2 < 320))",
    "xd^2 + (yd - beamOffset)^2 within 16cm^2 to 1100cm^2": "((x1_d**2 + (y1_d - beamOffset)**2 > 16) and (x1_d**2 + (y1_d - beamOffset)**2 < 1100)) and ((x2_d**2 + (y2_d - beamOffset)**2 > 16) and (x2_d**2 + (y2_d - beamOffset)**2 < 1100))",
    "abs(abs(px_st1 - px_st3) - 0.416) < 0.008GeV": "abs(abs(px1_st1 - px1_st3) - 0.416) < 0.008 and abs(abs(px2_st1 - px2_st3) - 0.416) < 0.008",
    "abs(py_st1 - py_st3) < 0.008GeV": "abs(py1_st1 - py1_st3) < 0.008 and abs(py2_st1 - py2_st3) < 0.008",
    "abs(pz_st1 - pz_st3) < 0.08GeV": "abs(pz1_st1 - pz1_st3) < 0.08 and abs(pz2_st1 - pz2_st3) < 0.08",
    "chisq1_target and chisq2_target < 15": "chisq1_target < 15 and chisq2_target < 15",
    "nHits Track1 and Track2 > 13": "nHits1 > 13 and nHits2 > 13",
    "chisq/NDF > 12": "(chisq1/(nHits1-5)) < 12 and (chisq2/(nHits2-5)) < 12",
    "pz1_st1 and pz2_st1 within 9 to 75": "pz1_st1 > 9 and pz1_st1 < 75 and pz2_st1 > 9 and pz2_st1 < 75",
    "y1_st1/y1_st3 and y2_st1/y2_st3 < 1": "(y1_st1/y1_st3) < 1 and (y2_st1/y2_st3) < 1",
    "y1_st1 x y1_st3 and y2_st1 x y2_st3 > 0": "(y1_st1 * y1_st3) > 0 and (y2_st1 * y2_st3) > 0",
    "abs(py1_st1) and abs(py2_st1) > 0.02": "abs(py1_st1) > 0.02 and abs(py2_st1) > 0.02",

    "dz within -280 to -5": "dz > -280 and dz < -5",
    "dx within -0.25 to 0.25": "dx > -0.25 and dx < 0.25",
    "dy - beamOffset within -0.22 to 0.22": "abs(dy - beamOffset) < 0.22",
    "dx^2 + (dy - beamOffset)^2 < 0.06":"(dx**2 + (dy - beamOffset)**2) < 0.06",
    "abs(dpx) < 1.8": "abs(dpx) < 1.8",
    "abs(dpy) < 2": "abs(dpy) < 2",
    "dpz within 38 to 116": "dpz > 38 and dpz < 116",
    "dpx^2 + dpy^2 < 5": "(dpx**2 + dpy**2) < 5",
    "mass within 4.2 to 8.8": "mass > 4.2 and mass < 8.8",
    "xF within -0.1 to 0.95": "xF > -0.1 and xF < 0.95",
    "xT within -0.1 to 0.58": "xT > -0.1 and xT < 0.58",
    "cosTheta within -0.5 to 0.5": "costh > -0.5 and costh < 0.5",
    "abs(trackSeparation) < 270": "abs(trackSeparation) < 270",
    "nhits1 + nhits2 > 29": "(nHits1 + nHits2) > 29",
    "nhits1st1 + nhits2st1 > 8": "(nHits1St1 + nHits2St1) > 8",
    "y1_st3 x y2_st3 < 0": "y1_st3 * y2_st3 < 0",
    "chisq_dimuon < 18": "chisq_dimuon < 18",
    "abs(x1_st1 + x2_st1) < 42cm": "abs(x1_st1 + x2_st1) < 42",
    "chisq Target within 2": "abs(chisq1_target + chisq2_target - chisq_dimuon) < 2",
    
    "D1 < 400": "D1 < 400",
    "D2 < 400": "D2 < 400",
    "D3 < 400": "D3 < 400",
    "D1 + D2 + D3 < 1000": "D1 + D2 + D3 < 1000",
    "intensity within 0 to 80000": "intensityP > 0 and intensityP < 80000"
} #: Dictionary defining the selection cuts. Keys are descriptive names, values are pandas evaluation strings.

variables_list = [
    "matrix1",
    "nHits1", 
    "nHits2", 
    "nHits1St1", 
    "nHits2St1",
    "chisq1", 
    "chisq2", 
    "chisq_dimuon", 
    "chisq1_target", 
    "chisq2_target",
    "dx", 
    "dy", 
    "dz", 
    "dpx", 
    "dpy", 
    "dpz", 
    "mass",
    "D1", 
    "D2", 
    "D3", 
    "xF", 
    "xT", 
    "xB", 
    "costh", 
    "intensityP",
    "runID", # Keep runID in case any cuts *still* use it
    "trackSeparation",
    "x1_t", 
    "x2_t", 
    "y1_t", 
    "y2_t", 
    "x1_d", 
    "x2_d", 
    "y1_d", 
    "y2_d", 
    "z1_v", 
    "z2_v",
    "x1_st1", 
    "x2_st1", 
    "y1_st1", 
    "y2_st1", 
    "y1_st3", 
    "y2_st3",
    "py1_st1", 
    "py2_st1", 
    "py1_st3", 
    "py2_st3", 
    "pz1_st1", 
    "pz2_st1",
    "pz1_st3", 
    "pz2_st3",
    "px1_st1", 
    "px2_st1", 
    "px1_st3", 
    "px2_st3",
    # 'beamOffset' is calculated, not read.
    # 'ReWeight' is handled automatically.
]#: List of essential variables to read from the ROOT trees.

# --- File Paths and Labels ---
# NOTE: Please ensure these paths are correct for your system
try:
    # Attempt to define paths assuming a certain directory structure
    base_path_hugo = "../../ROOTFiles/Hugo/" #: Base path for Hugo's ROOT files.
    base_path_mixed = "../../ROOTFiles/MixedEvents/" #: Base path for Mixed Event ROOT files.

    data_file = base_path_hugo + "roadset57_70_R008_2111v42_tmp_noPhys.root" #: Path to the main data file.

    mc_files = [
        base_path_hugo + "mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        base_path_hugo + "mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        base_path_hugo + "mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    ]#: List of paths to Monte Carlo simulation files.

    mixed_file = base_path_mixed + "merged_RS67_3089LH2.root" #: Path to the mixed event file.

except FileNotFoundError as fnf_error:
    print(f"❌ Configuration Error: File not found at expected path - {fnf_error}")
    print("Please check the base_path variables and file names.")
    sys.exit(1) # Exit if configuration paths are wrong
except NameError: # Handle if base paths aren't defined
    print("⚠️ Warning: Using placeholder file paths. Ensure paths are correct.")
    sys.stdout.flush()
    data_file = "path/to/your/data.root"
    mc_files = ["path/to/drellyan.root", "path/to/jpsi.root", "path/to/psiprime.root"]
    mixed_file = "path/to/mixed.root"


mc_labels_list = ["DY MC", "J/Psi MC", "Psi Prime MC"] #: Labels for the MC datasets.
tree = "Tree" #: Default TTree name for data and MC files.
mixed_tree = "result_mix" #: TTree name for mixed events in the mixed file.

# ------------------- Main Execution Logic -------------------

if __name__ == "__main__":
    # Generate the cut table by running the parallel processing function
    cut_table_df = create_cut_table_parallel(
        data_file,
        mc_files,
        mc_labels_list,
        tree,
        cuts_dict,
        variables_list,
        mixed_file,
        mixed_tree
    )

    # --- Display and Save ---
    # Check if the DataFrame was successfully created and is not empty
    if cut_table_df is not None and not cut_table_df.empty :
        display_df = cut_table_df.copy()

        # Define the desired column order for the final table
        desired_order = [
            "Data",
            "Data(RS67)",
            "Mixed(RS67)",
            "DY MC",
            "J/Psi MC",
            "Psi Prime MC",
            "Total MC",
            "Purity (DY MC)",
            "Efficiency (DY MC)"
        ] #: List defining the desired column order in the output table.

        # Filter the desired order list to only include columns that actually exist
        columns_to_display_ordered = [col for col in desired_order if col in display_df.columns]

        if not columns_to_display_ordered:
            print("❌ Error: No columns found to display after filtering. Check processing results.")
            sys.stdout.flush()
        else:
            # Reindex the DataFrame columns to match the desired order
            display_df = display_df[columns_to_display_ordered]

            # --- Formatting steps for display ---
            cut_names_order = ["Total Events"] + list(cuts_dict.keys())
            # Reindex the rows (index) to ensure correct cut order
            display_df = display_df.reindex(cut_names_order, fill_value=0.0)

            # Insert the "Cut" column at the beginning for display purposes
            display_df.insert(0, "Cut", display_df.index)

            # Define columns that should be formatted as integers (event counts)
            count_cols_base = ["Data", "Data(RS67)", "Mixed(RS67)", "Total MC", "DY MC", "J/Psi MC", "Psi Prime MC"]
            count_cols = [col for col in count_cols_base if col in display_df.columns]

            # Apply integer formatting to count columns
            for col in count_cols:
                 display_df[col] = pd.to_numeric(display_df[col], errors="coerce").fillna(0).round(0).astype(int)

            def format_value(value, col_name):
                """Formats a value for table display based on its column name.

                Parameters
                ----------
                value : number or str
                    The value from the DataFrame cell.
                col_name : str
                    The name of the column the value belongs to.

                Returns
                -------
                str
                    The formatted string representation of the value.
                """
                if col_name in count_cols:
                    return "{:.0f}".format(value) # Format as integer
                elif col_name in ["Purity (DY MC)", "Efficiency (DY MC)"]:
                     try:
                         # Format as percentage with 2 decimal places
                         return "{:.2f}".format(float(value))
                     except (ValueError, TypeError):
                          # Fallback if value is not numeric
                          return str(value)
                else:
                     # Default string conversion for any other columns
                     return str(value)

            # Prepare data in list-of-lists format for tabulate
            data_to_display = []
            for index, row in display_df.iterrows():
                # Create a formatted row, starting with the cut name
                formatted_row = [row["Cut"]] + [format_value(row[col], col) for col in display_df.columns if col != "Cut"]
                data_to_display.append(formatted_row)

            # Get headers from the (now ordered) display DataFrame columns
            headers = list(display_df.columns)

            # Print the final formatted table
            print("\n--- Final Cut Table ---")
            sys.stdout.flush()
            print(tabulate(data_to_display, headers=headers, tablefmt="fancy_grid"))
            sys.stdout.flush()

            # --- Save Results ---
            try:
                output_csv_file = "cut_table_with_mixed_and_reweighting_parallel.csv"
                # Save the original cut_table_df (column order might not match display)
                cut_table_df.to_csv(output_csv_file)
                print(f"\n✅ Cut table saved as '{output_csv_file}'")
                sys.stdout.flush()
                # Optional: Save the re-ordered and potentially formatted DataFrame
                # display_df.to_csv("cut_table_display_ordered.csv", index=False) # Save without index if 'Cut' column is used
            except Exception as e:
                print(f"\n❌ Error saving CSV file '{output_csv_file}': {e}")
                sys.stdout.flush()

    else:
        # Message if the DataFrame wasn't generated or is empty
        print("❌ Cut table DataFrame was not generated or is empty. Cannot display or save.")
        sys.stdout.flush()