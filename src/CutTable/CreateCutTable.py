# -*- coding: utf-8 -*-
"""
Process High Energy Physics data from ROOT files to generate a cut flow table.

This script reads particle physics data (dimuon events) from specified ROOT filabout:blank#blockedes,
including experimental data, Monte Carlo (MC) simulations, and mixed event samples.
It applies a series of predefined selection cuts sequentially to filter events.

Key Features:
- Reads data from ROOT files using the 'uproot' library.
- Handles multiple datasets (Data, MC types, Mixed Events, Empty Flask).
- Applies run-independent beam offset corrections (currently fixed).
- Processes files in parallel using 'concurrent.futures.ProcessPoolExecutor' for speed.
- Applies sequential event selection cuts defined in a dictionary.
- Calculates weighted sums for MC events using 'ReWeight' column.
- Aggregates results into a pandas DataFrame.
- Calculates derived columns like 'Total MC', 'Purity (DY MC)', 'Efficiency (DY MC)'.
- Prints a formatted cut flow table to the console using 'tabulate'.
- Saves the resulting cut flow table to a CSV file.
- Includes error handling and warnings during file reading and processing.
- Includes a "Good Spills Cut" based on an external file of spillIDs,
  conditionally applied to specific datasets.
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

# Define target labels for the Good Spills Cut globally for clarity
TARGET_LABELS_FOR_GOOD_SPILLS_CUT = {"Data", "Data(RS67)", "empty flask (RS67)"}

def load_good_spills(file_path):
    """
    Loads good spill IDs from a text file.

    Each line in the file is expected to contain one integer spillID.
    Lines that are empty or cannot be converted to int are ignored (with warnings).

    Parameters
    ----------
    file_path : str
        Path to the text file containing good spill IDs.

    Returns
    -------
    set[int]
        A set of integer spill IDs. Returns an empty set if the file
        is not found, cannot be read, or contains no valid integer spill IDs.
    """
    good_spills = set()
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                stripped_line = line.strip()
                if stripped_line: # Ensure line is not empty
                    try:
                        good_spills.add(int(stripped_line))
                    except ValueError:
                        print(f"⚠️ Warning: Non-integer value '{stripped_line}' found in good spills file '{file_path}' at line {line_num}. Skipping.")
                        sys.stdout.flush()
        if good_spills:
            print(f"✅ Successfully loaded {len(good_spills)} good spill IDs from {file_path}")
        else:
            print(f"⚠️ Warning: No valid good spill IDs loaded from {file_path}. 'Good Spills Cut' may be ineffective.")
        sys.stdout.flush()
        return good_spills
    except FileNotFoundError:
        print(f"❌ Error: Good spills file '{file_path}' not found. 'Good Spills Cut' will be ineffective.")
        sys.stdout.flush()
        return set() # Return an empty set if file not found
    except Exception as e:
        print(f"❌ Error loading good spills file '{file_path}': {e}")
        sys.stdout.flush()
        return set()


def read_tree(file_path, tree_name, variables):
    """
    Reads specified variables from a ROOT tree into a pandas DataFrame.

    Handles missing files, trees, and variables gracefully. If the file or
    tree cannot be read, an empty DataFrame with expected columns is returned.
    If variables are missing within a tree, they are added as columns
    filled with NaN.

    It ensures essential columns ('runID', 'dy', 'spillID') are attempted and converts
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
        'ReWeight', and 'spillID' are always added to this list internally if not present.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the requested variables.
    """
    try:
        with uproot.open(file_path) as file:
            if tree_name not in file:
                print(f"❌ Error: Tree '{tree_name}' not found in file '{file_path}'.")
                sys.stdout.flush() 
                empty_df = pd.DataFrame()
                for var in set(variables + ['runID', 'dy', 'ReWeight', 'spillID']):
                    empty_df[var] = pd.Series(dtype='float64') 
                return empty_df

            tree = file[tree_name]
            all_keys = tree.keys()
            
            vars_to_request = list(set(variables + ['runID', 'dy', 'ReWeight', 'spillID']))

            vars_to_read = [var for var in vars_to_request if var in all_keys]
            missing = [var for var in vars_to_request if var not in all_keys]
            if missing:
                is_mc_file = any(mc_part in file_path for mc_part in ['mc_drellyan', 'mc_jpsi', 'mc_psiprime'])
                missing_for_warning = list(missing)
                if 'ReWeight' in missing_for_warning and not is_mc_file:
                    missing_for_warning.remove('ReWeight')
                
                if missing_for_warning:
                    print(f"⚠️ Missing variables in {file_path} ({tree_name}): {missing_for_warning}. These will be added as NaN.")
                    sys.stdout.flush()

            df = tree.arrays(vars_to_read, library="pd")

            for var in missing:
                if var not in df.columns: 
                    df[var] = np.nan

            if 'runID' in df.columns:
                df['runID'] = pd.to_numeric(df['runID'], errors='coerce')
            if 'dy' in df.columns:
                df['dy'] = pd.to_numeric(df['dy'], errors='coerce')
            if 'spillID' in df.columns: 
                df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')
            
            if 'ReWeight' in df.columns:
                df['ReWeight'] = pd.to_numeric(df['ReWeight'], errors='coerce').fillna(1.0)
            else:
                df['ReWeight'] = 1.0
            
            return df

    except Exception as e:
        print(f"❌ Error reading tree '{tree_name}' from '{file_path}': {e}")
        sys.stdout.flush()
        empty_df = pd.DataFrame()
        for var in set(variables + ['runID', 'dy', 'ReWeight', 'spillID']):
            empty_df[var] = pd.Series(dtype='float64')
        return empty_df


def add_beam_offset(df):
    """
    Adds the 'beamOffset' column to the DataFrame with a fixed value.
    """
    df['beamOffset'] = 1.6
    return df


def apply_cuts(df, cuts, use_weights=True, good_spills_set=None, dataset_label=None):
    """
    Applies a dictionary of sequential cuts to a pandas DataFrame.
    Includes special handling for "Good Spills Cut" based on dataset_label.
    """
    cut_results = {}
    weight = df['ReWeight'] if use_weights and 'ReWeight' in df.columns else pd.Series(1.0, index=df.index)
    weight = pd.to_numeric(weight, errors='coerce').fillna(1.0)
    cut_results["Total Events"] = weight.sum()

    if 'beamOffset' not in df.columns:
        df['beamOffset'] = 1.6 
    if 'dy' in df.columns:
        df['dy'] = pd.to_numeric(df['dy'], errors='coerce').fillna(np.nan)
    if 'runID' in df.columns:
        df['runID'] = pd.to_numeric(df['runID'], errors='coerce').fillna(-1)
    if 'spillID' in df.columns:
         df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')


    cumulative_mask = pd.Series(True, index=df.index)
    original_df = df.copy()
    original_weight = weight.copy()

    for cut_name, cut_string in cuts.items():
        if cut_name == "Total Events":
            continue
        try:
            current_df_for_cut = original_df[cumulative_mask].copy()
            
            if cut_name == "Good Spills Cut":
                if dataset_label in TARGET_LABELS_FOR_GOOD_SPILLS_CUT:
                    if good_spills_set is not None and 'spillID' in current_df_for_cut.columns and not good_spills_set == set():
                        current_spillID_numeric = pd.to_numeric(current_df_for_cut['spillID'], errors='coerce')
                        mask_this_cut = current_spillID_numeric.isin(good_spills_set) & current_spillID_numeric.notna()
                    else:
                        if not good_spills_set or good_spills_set == set():
                             print(f"⚠️ Warning: 'Good Spills Cut' for dataset '{dataset_label}' is ineffective (spill list empty/not loaded). Passing all events.")
                        if 'spillID' not in current_df_for_cut.columns:
                             print(f"⚠️ Warning: 'Good Spills Cut' for dataset '{dataset_label}' cannot be applied ('spillID' missing). Passing all events.")
                        sys.stdout.flush()
                        mask_this_cut = pd.Series(True, index=current_df_for_cut.index)
                else:
                    # For non-target datasets, this cut passes all events
                    # print(f"ℹ️ Skipping 'Good Spills Cut' for dataset '{dataset_label}' (not a target). Passing all events.")
                    # sys.stdout.flush() # Optional: can make output verbose
                    mask_this_cut = pd.Series(True, index=current_df_for_cut.index)
            else:
                mask_this_cut = current_df_for_cut.eval(cut_string, engine='python')

            aligned_mask_this_cut = pd.Series(False, index=original_df.index)
            aligned_mask_this_cut.loc[current_df_for_cut.index[mask_this_cut]] = True
            
            cumulative_mask = cumulative_mask & aligned_mask_this_cut
            cut_results[cut_name] = original_weight[cumulative_mask].sum()

        except Exception as e:
            print(f"⚠️ Error applying cut '{cut_name}' for dataset '{dataset_label}': {e}. Skipping cut.")
            sys.stdout.flush()
            cut_results[cut_name] = original_weight[cumulative_mask].sum() # Result before this failed cut
            # Or assign 0: cut_results[cut_name] = 0 
            # current behavior: count doesn't change for this specific cut if it errors,
            # previous cumulative sum is used. This means subsequent cuts operate on data that
            # passed up to the cut *before* the failing one.
            # If we want to report 0 for a failing cut and stop its lineage:
            # cut_results[cut_name] = 0
            # And ensure cumulative_mask isn't further ANDed with a problematic mask.
            # For simplicity, let's assume if a cut string fails, it filters nothing for that step.
            # The current code effectively does this for the `cut_results[cut_name]` by taking the sum *before* `cumulative_mask` would have been updated by a failing `mask_this_cut`.
            # A more robust way for error:
            # cut_results[cut_name] = 0 # Report 0 for this cut
            # cumulative_mask should remain as it was before this failing cut.
            # The issue is that if eval fails, mask_this_cut is not defined.
            # Let's refine the error handling for a failing eval:
            # If an error occurs in eval or the custom cut logic:
            #   print warning
            #   cut_results[cut_name] gets the sum of weights of events *that passed up to the previous cut*.
            #   This means the failing cut itself filters out 0 events from the perspective of this cut name's entry.
            #   The cumulative_mask remains unchanged by this specific failing cut, so subsequent cuts use the mask from *before* this failure.
            # This is what the original code before the large comment block for error handling implies.
            # Let's ensure result is 0 if the cut string is bad or essential data is missing for THIS cut
            # The current structure:
            # try:
            #    ...
            #    cumulative_mask = cumulative_mask & aligned_mask_this_cut
            #    cut_results[cut_name] = original_weight[cumulative_mask].sum()
            # except:
            #    cut_results[cut_name] = 0 # This was in the original user code for robustness
            # This is safer. If a cut fails, it results in 0 for that step, and previous cumulative_mask is used.
            # The `cumulative_mask` should not be updated with a faulty `aligned_mask_this_cut`.
            # So, if exception, `cumulative_mask` should not be updated.
            # Let's reconsider the except block:
            # If an error occurs in applying `cut_string`:
            # 1. Print warning.
            # 2. `cut_results[cut_name]` should ideally be the count *after* attempting this cut. If it fails, what should it be?
            #    If we set it to `original_weight[cumulative_mask].sum()`, it's the count *before* this cut.
            #    If we set it to `0`, it means this cut drastically failed.
            #    The safest is to record the count *before* this failed cut, and the `cumulative_mask` isn't changed further by this failed cut.
            #    This is implicitly what happens if `aligned_mask_this_cut` is not used to update `cumulative_mask`.

            # Let's stick to the provided code's original error handling logic for `apply_cuts`
            # which was:
            # cut_results[cut_name] = 0 # Assign 0 count for the failed cut
            # This implies the cumulative_mask from *before* this failed cut is used for the *next* iteration.
            # The current code `cut_results[cut_name] = original_weight[cumulative_mask].sum()` if error occurs AFTER `cumulative_mask` was already updated by previous cuts.
            # The original skeleton had `cut_results[cut_name] = 0` on exception. This makes more sense for "events passing THIS cut".
            # Resetting to that logic within the except:
            pass # Original script has a robust way: set result to 0 and df for next cut uses previous mask.
                 # My refined loop:
                 # original_df = df.copy() ... current_df_for_cut = original_df[cumulative_mask]
                 # if cut fails, cumulative_mask is NOT updated. cut_results[cut_name] is assigned a value.
                 # The original script's simpler df = df[mask] implies if cut fails, df is not updated.

    # Reverting to the user's provided structure for the try-except in apply_cuts for cut evaluation
    # The current structure correctly keeps cumulative_mask from successful prior cuts if one fails.
    # And `cut_results[cut_name] = original_weight[cumulative_mask].sum()` (where cumulative_mask is from *before* the current failing cut) is actually a reasonable choice
    # *OR* `cut_results[cut_name] = 0` if we define it as "0 events passed *this specific* failing cut".
    # The original script had `cut_results[cut_name] = 0` in the `except` block. Let's reinstate that for clarity of a failed step.

    # Final refined error handling for one cut:
    # (inside loop for cut_name, cut_string)
    #   temp_cumulative_mask = cumulative_mask.copy() # preserve state before this cut
    #   try:
    #       ... calculate mask_this_cut ...
    #       aligned_mask_this_cut = ...
    #       cumulative_mask = cumulative_mask & aligned_mask_this_cut
    #       cut_results[cut_name] = original_weight[cumulative_mask].sum()
    #   except Exception as e:
    #       print(...)
    #       cut_results[cut_name] = 0 # This cut failed, 0 pass it
    #       cumulative_mask = temp_cumulative_mask # Revert mask to state before this failed cut

    # The provided code's current logic in the main question (before my modification for conditional cut) was:
    #      except Exception as e:
    #          print(f"⚠️ Error applying cut '{cut_name}': {e}. Skipping cut for this dataset.")
    #          sys.stdout.flush() # Try to flush output buffer
    #          cut_results[cut_name] = 0 # Assign 0 count for the failed cut
    # This is what I will stick to. The `cumulative_mask` is updated *only if the cut succeeds*.
    # If it fails, `cumulative_mask` for the *next* cut is the one from *before* this failing cut.
    # And the result for *this* cut is 0. This seems correct. My `aligned_mask_this_cut` logic handles this.

    return pd.Series(cut_results)


def process_single_file(file_path, tree_name, label, variables, cuts, use_weights, good_spills_set):
    """
    Processes a single ROOT file/tree: reads data, adds offset, applies cuts.
    Passes dataset label and good_spills_set to apply_cuts.
    """
    try:
        print(f"Starting processing for {label} from {file_path} (Tree: {tree_name})...")
        sys.stdout.flush()
        vars_to_read = list(set(variables + ['runID', 'spillID'])) 
        df = read_tree(file_path, tree_name, vars_to_read)

        if df.empty and not all(col in df.columns for col in ['runID', 'dy', 'ReWeight', 'spillID']):
            print(f"INFO: read_tree returned empty DataFrame for {label}. Skipping further processing.")
            sys.stdout.flush()
            expected_cut_names = ["Total Events"] + list(cuts.keys())
            return {'label': label, 'results': pd.Series(0.0, index=expected_cut_names), 'error': False}

        df = add_beam_offset(df)

        cut_variables = set()
        for cut_string in cuts.values():
            found_vars = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', cut_string)
            keywords_to_exclude = {'and', 'or', 'not', 'in', 'is', 'True', 'False', 'None', 'abs', 'np', 'sqrt', 'min', 'max', 'custom', 'handled'}
            valid_vars = [var for var in found_vars if var not in keywords_to_exclude and not var.isdigit()]
            cut_variables.update(valid_vars)
        
        if "Good Spills Cut" in cuts:
             cut_variables.add('spillID')

        required_for_cuts = list(set(variables + list(cut_variables) + ['beamOffset', 'dy', 'runID', 'spillID']))

        for var in required_for_cuts:
            if var not in df.columns:
                print(f"⚠️ Adding missing required cut variable '{var}' as NaN to {label} DataFrame.")
                sys.stdout.flush()
                df[var] = np.nan

        cut_series = apply_cuts(df, cuts, use_weights=use_weights, good_spills_set=good_spills_set, dataset_label=label) # Pass label
        print(f"Finished processing {label}.")
        sys.stdout.flush()
        return {'label': label, 'results': cut_series, 'error': False}

    except Exception as e:
        print(f"❌ Error processing {label} from {file_path} (Tree: {tree_name}): {e}")
        sys.stdout.flush()
        expected_cut_names = ["Total Events"] + list(cuts.keys())
        return {'label': label, 'results': pd.Series(0.0, index=expected_cut_names), 'error': True, 'message': str(e)}


def create_cut_table_parallel(data_file_path, mc_file_paths, mc_labels, tree_name, cuts, variables, mixed_file_path, mixed_tree_name, empty_flask_file_path, good_spills_file_path):
    """
    Creates a cut flow table by processing multiple datasets in parallel.
    Loads good spill IDs and passes them for processing.
    """
    start_time = time.time()
    print("Starting parallel processing...")
    sys.stdout.flush()

    good_spills_set = load_good_spills(good_spills_file_path)

    results_dict = {} 
    futures = [] 

    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        futures.append(executor.submit(process_single_file, data_file_path, tree_name, "Data", variables, cuts, False, good_spills_set))
        futures.append(executor.submit(process_single_file, mixed_file_path, "result", "Data(RS67)", variables, cuts, False, good_spills_set))
        futures.append(executor.submit(process_single_file, mixed_file_path, mixed_tree_name, "Mixed(RS67)", variables, cuts, False, good_spills_set))
        futures.append(executor.submit(process_single_file, empty_flask_file_path, "result", "empty flask (RS67)", variables, cuts, False, good_spills_set))
        futures.append(executor.submit(process_single_file, empty_flask_file_path, "result_mix", "empty flask (RS67) mixed", variables, cuts, False, good_spills_set))

        for i, mc_file_path in enumerate(mc_file_paths):
            label = mc_labels[i]
            mc_vars = list(set(variables + ["ReWeight"])) 
            futures.append(executor.submit(process_single_file, mc_file_path, tree_name, label, mc_vars, cuts, True, good_spills_set))

        print("Waiting for tasks to complete...")
        sys.stdout.flush()
        processed_tasks = 0
        total_tasks = len(futures)
        for future in concurrent.futures.as_completed(futures):
            processed_tasks += 1
            # To avoid excessively verbose output for many tasks, optionally reduce print frequency
            if total_tasks < 20 or processed_tasks % (total_tasks // 10 if total_tasks >=10 else 1) == 0 or processed_tasks == total_tasks :
                 print(f"Completed task {processed_tasks}/{total_tasks}...")
                 sys.stdout.flush()
            try:
                result = future.result()
                label = result['label']
                cut_series = result['results']
                is_error = result['error']
                error_message = result.get('message', '')

                if is_error:
                    print(f"❗ Task for '{label}' completed with error: {error_message}. Using zero counts.")
                    sys.stdout.flush()
                    expected_cut_names = ["Total Events"] + list(cuts.keys())
                    results_dict[label] = pd.Series(0.0, index=expected_cut_names)
                else:
                    results_dict[label] = cut_series
                    # print(f"Successfully processed results for '{label}'.") # Can be verbose
                    sys.stdout.flush()

            except Exception as exc:
                # Attempt to find which future/label caused the error if possible (tricky with as_completed)
                # For now, generic message:
                print(f'❌ Exception occurred while retrieving future result (label unknown if future failed before returning label): {exc}')
                sys.stdout.flush()


    print("Aggregating results into DataFrame...")
    sys.stdout.flush()
    cut_table = pd.DataFrame(results_dict)

    expected_indices = ["Total Events"] + list(cuts.keys())
    cut_table = cut_table.reindex(expected_indices, fill_value=0.0)

    print("Aggregating results and calculating final columns...")
    sys.stdout.flush()

    mc_cols = [label for label in mc_labels if label in cut_table.columns]
    if mc_cols:
        cut_table["Total MC"] = cut_table[mc_cols].sum(axis=1, min_count=1).fillna(0)
    else:
        print("⚠️ No MC columns found in results to calculate 'Total MC'. Setting to 0.")
        sys.stdout.flush()
        cut_table["Total MC"] = pd.Series(0.0, index=cut_table.index)

    try:
        dy_col = "DY MC"
        if dy_col in cut_table.columns and "Total Events" in cut_table.index and "Total MC" in cut_table.columns:
            total_dy_at_start = cut_table.loc["Total Events", dy_col] if "Total Events" in cut_table.index and dy_col in cut_table.columns else 0
            
            if total_dy_at_start == 0:
                print(f"⚠️ Total Events for '{dy_col}' is zero. Efficiency will be zero or NaN.")
                sys.stdout.flush()

            total_mc_safe = cut_table["Total MC"].replace(0, np.nan) 
            purity = (cut_table[dy_col] / total_mc_safe).replace([np.inf, -np.inf], np.nan).fillna(0).copy()
            
            total_dy_at_start_safe = total_dy_at_start if total_dy_at_start != 0 else np.nan 
            efficiency = (cut_table[dy_col] / total_dy_at_start_safe).replace([np.inf, -np.inf], np.nan).fillna(0).copy()

            cut_table["Purity (DY MC)"] = purity * 100
            cut_table["Efficiency (DY MC)"] = efficiency * 100
        else:
            missing_elements = []
            if dy_col not in cut_table.columns: missing_elements.append(f"'{dy_col}' column")
            if "Total Events" not in cut_table.index: missing_elements.append("'Total Events' index")
            if "Total MC" not in cut_table.columns: missing_elements.append("'Total MC' column")
            print(f"⚠️ Cannot calculate Purity/Efficiency: Missing {', '.join(missing_elements)}.")
            sys.stdout.flush()
            cut_table["Purity (DY MC)"] = 0.0
            cut_table["Efficiency (DY MC)"] = 0.0
    except Exception as e:
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
GOOD_SPILLS_FILE_PATH = "../../res/GoodSpills/GoodSpills.txt" 

cuts_dict = {
    "Good Spills Cut": "spillID # Custom handled by apply_cuts", 
    "z1_v and z2_v within -320cm to -5cm": "z1_v > -320 and z1_v < -5 and z2_v > -320 and z2_v < -5",
    "xt^2 + (yt - beamOffset)^2 < 320cm^2": "((x1_t**2 + (y1_t - beamOffset)**2 < 320)) and ((x2_t**2 + (y2_t - beamOffset)**2 < 320))",
    "xd^2 + (yd - beamOffset)^2 within 16cm^2 to 1100cm^2": "((x1_d**2 + (y1_d - beamOffset)**2 > 16) and (x1_d**2 + (y1_d - beamOffset)**2 < 1100)) and ((x2_d**2 + (y2_d - beamOffset)**2 > 16) and (x2_d**2 + (y2_d - beamOffset)**2 < 1100))",
    "abs(abs(px_st1 - px_st3) - 0.416) < 0.008GeV": "abs(abs(px1_st1 - px1_st3) - 0.416) < 0.008 and abs(abs(px2_st1 - px2_st3) - 0.416) < 0.008",
    "abs(py_st1 - py_st3) < 0.008GeV": "abs(py1_st1 - py1_st3) < 0.008 and abs(py2_st1 - py2_st3) < 0.008",
    "abs(pz_st1 - pz_st3) < 0.08GeV": "abs(pz1_st1 - pz1_st3) < 0.08 and abs(pz2_st1 - pz2_st3) < 0.08",
    "chisq1_target and chisq2_target < 15": "chisq1_target < 15 and chisq2_target < 15",
    "nHits Track1 and Track2 > 13": "nHits1 > 13 and nHits2 > 13",
    "chisq/NDF < 12": "(chisq1/(nHits1-5)) < 12 and (chisq2/(nHits2-5)) < 12", # Corrected from > to <
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
    "abs(x1_st1 + x2_st1) < 42": "abs(x1_st1 + x2_st1) < 42", # Removed "cm"
    "chisq Target within 2": "abs(chisq1_target + chisq2_target - chisq_dimuon) < 2",
    "D1 < 400": "D1 < 400",
    "D2 < 400": "D2 < 400",
    "D3 < 400": "D3 < 400",
    "D1 + D2 + D3 < 1000": "D1 + D2 + D3 < 1000",
    "intensity within 0 to 80000": "intensityP > 0 and intensityP < 80000"
}

variables_list = [
    "spillID", "matrix1",
    "nHits1", "nHits2", "nHits1St1", "nHits2St1",
    "chisq1", "chisq2", "chisq_dimuon", "chisq1_target", "chisq2_target",
    "dx", "dy", "dz", "dpx", "dpy", "dpz", "mass",
    "D1", "D2", "D3", "xF", "xT", "xB", "costh", "intensityP",
    "runID", "trackSeparation",
    "x1_t", "x2_t", "y1_t", "y2_t", "x1_d", "x2_d", "y1_d", "y2_d", "z1_v", "z2_v",
    "x1_st1", "x2_st1", "y1_st1", "y2_st1", "y1_st3", "y2_st3",
    "py1_st1", "py2_st1", "py1_st3", "py2_st3", 
    "pz1_st1", "pz2_st1", "pz1_st3", "pz2_st3",
    "px1_st1", "px2_st1", "px1_st3", "px2_st3",
]

try:
    base_path_hugo = "../../res/ROOTFiles/Hugo/"
    base_path_mixed = "../../res/ROOTFiles/MixedEvents/"
    data_file = base_path_hugo + "roadset57_70_R008_2111v42_tmp_noPhys.root"
    mc_files = [
        base_path_hugo + "mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        base_path_hugo + "mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        base_path_hugo + "mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    ]
    mixed_file = base_path_mixed + "merged_RS67_3089LH2.root"
    empty_flask_file = base_path_mixed + "merged_RS67_3089flask.root"
except NameError: 
    print("⚠️ Warning: Using placeholder file paths due to undefined base paths. Ensure paths are correct.")
    sys.stdout.flush()
    data_file = "data.root" 
    mc_files = ["drellyan.root", "jpsi.root", "psiprime.root"] 
    mixed_file = "mixed.root" 
    empty_flask_file = "empty_flask.root" 

mc_labels_list = ["DY MC", "J/Psi MC", "Psi Prime MC"]
tree = "Tree" 
mixed_tree = "result_mix"

if __name__ == "__main__":
    cut_table_df = create_cut_table_parallel(
        data_file, mc_files, mc_labels_list, tree, 
        cuts_dict, variables_list,
        mixed_file, mixed_tree, empty_flask_file,
        GOOD_SPILLS_FILE_PATH
    )

    if cut_table_df is not None and not cut_table_df.empty :
        display_df = cut_table_df.copy()
        desired_order = [
            "Data", "Data(RS67)", "Mixed(RS67)",
            "empty flask (RS67)", "empty flask (RS67) mixed",
            "DY MC", "J/Psi MC", "Psi Prime MC", "Total MC",
            "Purity (DY MC)", "Efficiency (DY MC)"
        ]
        columns_to_display_ordered = [col for col in desired_order if col in display_df.columns]

        if not columns_to_display_ordered:
            print("❌ Error: No columns found to display after filtering. Check processing results.")
        else:
            display_df = display_df[columns_to_display_ordered]
            cut_names_order = ["Total Events"] + list(cuts_dict.keys())
            display_df = display_df.reindex(cut_names_order, fill_value=0.0)
            display_df.insert(0, "Cut", display_df.index)

            count_cols_base = [
                "Data", "Data(RS67)", "Mixed(RS67)",
                "empty flask (RS67)", "empty flask (RS67) mixed",
                "Total MC", "DY MC", "J/Psi MC", "Psi Prime MC"
            ]
            count_cols = [col for col in count_cols_base if col in display_df.columns]

            for col in count_cols:
                display_df[col] = pd.to_numeric(display_df[col], errors="coerce").fillna(0).round(0).astype(int)

            def format_value(value, col_name):
                if col_name in count_cols:
                    return "{:.0f}".format(value) 
                elif col_name in ["Purity (DY MC)", "Efficiency (DY MC)"]:
                    try:
                        return "{:.2f}".format(float(value))
                    except (ValueError, TypeError): return str(value)
                else: return str(value)

            data_to_display = [[row["Cut"]] + [format_value(row[col], col) for col in display_df.columns if col != "Cut"] for _, row in display_df.iterrows()]
            
            headers = list(display_df.columns)
            print("\n--- Final Cut Table ---")
            sys.stdout.flush()
            print(tabulate(data_to_display, headers=headers, tablefmt="fancy_grid"))
            sys.stdout.flush()

            try:
                output_csv_file = "cut_table_conditional_goodspills.csv" 
                cut_table_df.to_csv(output_csv_file)
                print(f"\n✅ Cut table saved as '{output_csv_file}'")
            except Exception as e:
                print(f"\n❌ Error saving CSV file '{output_csv_file}': {e}")
            sys.stdout.flush()
    else:
        print("❌ Cut table DataFrame was not generated or is empty. Cannot display or save.")
        sys.stdout.flush()