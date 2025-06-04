# -*- coding: utf-8 -*-
"""
Process High Energy Physics data from ROOT files to generate a cut flow table
with statistical uncertainties and MC normalization.

This script reads particle physics data (dimuon events) from specified ROOT filabout:blank#blockedes,
including experimental data, Monte Carlo (MC) simulations, and mixed event samples.
It applies a series of predefined selection cuts sequentially to filter events.
MC samples are scaled by normalization constants, and uncertainties are propagated.

Key Features:
- Reads data from ROOT files using the 'uproot' library.
- Handles multiple datasets (Data, MC types, Mixed Events, Empty Flask).
- Applies run-independent beam offset corrections (currently fixed).
- Processes files in parallel using 'concurrent.futures.ProcessPoolExecutor' for speed.
- Applies sequential event selection cuts defined in a dictionary.
- Calculates weighted sums and sum of weights squared for MC events using 'ReWeight' column.
- Applies normalization constants to specified MC columns and propagates uncertainties.
- Calculates counts for data events.
- Aggregates results into pandas DataFrames.
- Calculates derived columns like 'Total MC', 'Purity (DY MC)', 'Efficiency (DY MC)'.
- Prints a formatted cut flow table to the console using 'tabulate', showing value ± uncertainty.
- Saves the resulting counts and sum-of-weights-squared tables to CSV files.
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
                df['ReWeight'] = 1.0 # Default to 1.0 if not present (important for use_weights logic)
            
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
    Returns two series: one for counts/sum_of_weights, and one for N/sum_of_weights_squared.
    Includes special handling for "Good Spills Cut" based on dataset_label.
    """
    cut_counts = {}
    cut_sum_w2 = {} 

    # Define base_weight_series based on use_weights and ReWeight column presence
    if use_weights and 'ReWeight' in df.columns:
        df['ReWeight'] = pd.to_numeric(df['ReWeight'], errors='coerce').fillna(1.0)
        base_weight_series = df['ReWeight']
    else:
        base_weight_series = pd.Series(1.0, index=df.index)

    cut_counts["Total Events"] = base_weight_series.sum()
    if use_weights: 
        cut_sum_w2["Total Events"] = (base_weight_series**2).sum()
    else: 
        cut_sum_w2["Total Events"] = base_weight_series.sum()

    if 'beamOffset' not in df.columns: df['beamOffset'] = 1.6 
    if 'dy' in df.columns: df['dy'] = pd.to_numeric(df['dy'], errors='coerce').fillna(np.nan)
    if 'runID' in df.columns: df['runID'] = pd.to_numeric(df['runID'], errors='coerce').fillna(-1)
    if 'spillID' in df.columns and not pd.api.types.is_numeric_dtype(df['spillID']):
         df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')


    cumulative_mask = pd.Series(True, index=df.index) 

    for cut_name, cut_string in cuts.items():
        if cut_name == "Total Events": 
            continue
        
        mask_for_this_specific_cut = pd.Series(False, index=df.index) 

        try:
            if cut_name == "Good Spills Cut":
                if dataset_label in TARGET_LABELS_FOR_GOOD_SPILLS_CUT:
                    if good_spills_set is not None and 'spillID' in df.columns and not good_spills_set == set():
                        if not pd.api.types.is_numeric_dtype(df['spillID']): 
                            df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')
                        mask_for_this_specific_cut = df['spillID'].isin(good_spills_set) & df['spillID'].notna()
                        mask_for_this_specific_cut = mask_for_this_specific_cut.reindex(df.index, fill_value=False)
                    else:
                        if not good_spills_set or good_spills_set == set():
                             print(f"⚠️ Warning: 'Good Spills Cut' for dataset '{dataset_label}' is ineffective (spill list empty/not loaded). Passing all events.")
                        if 'spillID' not in df.columns:
                             print(f"⚠️ Warning: 'Good Spills Cut' for dataset '{dataset_label}' cannot be applied ('spillID' missing). Passing all events.")
                        sys.stdout.flush()
                        mask_for_this_specific_cut = pd.Series(True, index=df.index)
                else: 
                    mask_for_this_specific_cut = pd.Series(True, index=df.index)
            else:
                mask_for_this_specific_cut = df.eval(cut_string, engine='python')

            if not isinstance(mask_for_this_specific_cut, pd.Series) or mask_for_this_specific_cut.dtype != bool:
                 raise TypeError(f"Cut '{cut_name}' for '{dataset_label}' did not produce a boolean Series. Got type: {type(mask_for_this_specific_cut)}")
            if not mask_for_this_specific_cut.index.equals(df.index):
                 print(f"⚠️ Warning: Index mismatch for cut '{cut_name}' on dataset '{dataset_label}'. Reindexing mask.")
                 sys.stdout.flush()
                 mask_for_this_specific_cut = mask_for_this_specific_cut.reindex(df.index, fill_value=False)

            new_cumulative_mask = cumulative_mask & mask_for_this_specific_cut
            
            current_weights_passing = base_weight_series[new_cumulative_mask]
            cut_counts[cut_name] = current_weights_passing.sum()
            if use_weights:
                cut_sum_w2[cut_name] = (current_weights_passing**2).sum()
            else: 
                cut_sum_w2[cut_name] = current_weights_passing.sum() 
            
            cumulative_mask = new_cumulative_mask 

        except Exception as e:
            print(f"⚠️ Error applying cut '{cut_name}' for dataset '{dataset_label}': {e}. Reporting 0 for this cut.")
            sys.stdout.flush()
            cut_counts[cut_name] = 0.0
            cut_sum_w2[cut_name] = 0.0
            
    return pd.Series(cut_counts), pd.Series(cut_sum_w2)


def process_single_file(file_path, tree_name, label, variables, cuts, use_weights, good_spills_set):
    """
    Processes a single ROOT file/tree: reads data, adds offset, applies cuts.
    Passes dataset label and good_spills_set to apply_cuts.
    Returns counts and sum_of_weights_squared (or N for data).
    """
    try:
        print(f"Starting processing for {label} from {file_path} (Tree: {tree_name})...")
        sys.stdout.flush()
        
        vars_to_read = list(set(variables + ['runID', 'spillID'])) 
        if use_weights:
            vars_to_read = list(set(vars_to_read + ['ReWeight']))
        
        df = read_tree(file_path, tree_name, vars_to_read)

        essential_cols_for_df_check = ['runID'] 
        if df.empty and not all(col in df.columns for col in essential_cols_for_df_check):
            print(f"INFO: read_tree returned a critically empty DataFrame for {label}. Skipping further processing and returning zeros.")
            sys.stdout.flush()
            expected_cut_names = ["Total Events"] + list(cuts.keys())
            return {
                'label': label, 
                'results_counts': pd.Series(0.0, index=expected_cut_names),
                'results_sum_w2': pd.Series(0.0, index=expected_cut_names),
                'error': False 
            }

        df = add_beam_offset(df) 

        cut_variables = set()
        keywords_to_exclude_for_var_extraction = {
            'and', 'or', 'not', 'in', 'is', 'True', 'False', 'None', 
            'abs', 'np', 'sqrt', 'min', 'max', 'pow', 'float', 'int' 
        }

        for cut_name_local, cut_string_original in cuts.items(): 
            cut_string_for_parsing = "" 
            if cut_name_local == "Good Spills Cut":
                cut_string_for_parsing = "spillID" 
            else:
                cut_string_for_parsing = cut_string_original.split('#')[0].strip()

            if not cut_string_for_parsing: 
                continue

            found_vars = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', cut_string_for_parsing)
            
            valid_vars = [var for var in found_vars 
                          if var not in keywords_to_exclude_for_var_extraction and \
                          not var.isdigit()]
            cut_variables.update(valid_vars)
        
        all_required_vars_for_df_check = list(set(variables + list(cut_variables) + 
                                                  ['beamOffset', 'dy', 'runID', 'spillID']))
        if use_weights:
            all_required_vars_for_df_check = list(set(all_required_vars_for_df_check + ['ReWeight']))

        for var in all_required_vars_for_df_check:
            if var not in df.columns:
                print(f"⚠️ Adding missing required variable '{var}' as NaN to '{label}' DataFrame before applying cuts. Check variable list and ROOT file contents if this is unexpected.")
                sys.stdout.flush()
                df[var] = np.nan
            elif df[var].dtype == 'object': 
                 df[var] = pd.to_numeric(df[var], errors='coerce')

        cut_series_counts, cut_series_sum_w2 = apply_cuts(df, cuts, use_weights=use_weights, good_spills_set=good_spills_set, dataset_label=label)
        
        print(f"Finished processing {label}.")
        sys.stdout.flush()
        return {
            'label': label, 
            'results_counts': cut_series_counts, 
            'results_sum_w2': cut_series_sum_w2, 
            'error': False
        }

    except Exception as e:
        print(f"❌ Error processing {label} from {file_path} (Tree: {tree_name}): {e}")
        import traceback
        traceback.print_exc() 
        sys.stdout.flush()
        expected_cut_names = ["Total Events"] + list(cuts.keys())
        return {
            'label': label, 
            'results_counts': pd.Series(0.0, index=expected_cut_names), 
            'results_sum_w2': pd.Series(0.0, index=expected_cut_names), 
            'error': True, 
            'message': str(e)
        }


def create_cut_table_parallel(data_file_path, mc_file_paths, mc_labels, tree_name, cuts, variables, mixed_file_path, mixed_tree_name, empty_flask_file_path, good_spills_file_path):
    """
    Creates a cut flow table by processing multiple datasets in parallel.
    Applies normalization to MC columns and propagates uncertainties.
    """
    start_time = time.time()
    print("Starting parallel processing...")
    sys.stdout.flush()

    good_spills_set = load_good_spills(good_spills_file_path)

    results_counts_dict = {} 
    results_sum_w2_dict = {}
    futures = [] 

    expected_cut_names = ["Total Events"] + list(cuts.keys())

    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        futures.append(executor.submit(process_single_file, data_file_path, tree_name, "Data", variables, cuts, False, good_spills_set))
        futures.append(executor.submit(process_single_file, mixed_file_path, "result", "Data(RS67)", variables, cuts, False, good_spills_set)) 
        futures.append(executor.submit(process_single_file, mixed_file_path, mixed_tree_name, "Mixed(RS67)", variables, cuts, False, good_spills_set))
        futures.append(executor.submit(process_single_file, empty_flask_file_path, "result", "empty flask (RS67)", variables, cuts, False, good_spills_set))
        futures.append(executor.submit(process_single_file, empty_flask_file_path, "result_mix", "empty flask (RS67) mixed", variables, cuts, False, good_spills_set))

        for i, mc_file_path in enumerate(mc_file_paths):
            label = mc_labels[i]
            futures.append(executor.submit(process_single_file, mc_file_path, tree_name, label, variables, cuts, True, good_spills_set))

        print("Waiting for tasks to complete...")
        sys.stdout.flush()
        processed_tasks = 0
        total_tasks = len(futures)
        
        for future in concurrent.futures.as_completed(futures):
            processed_tasks += 1
            if total_tasks < 20 or processed_tasks % (total_tasks // 10 if total_tasks >=10 else 1) == 0 or processed_tasks == total_tasks :
                print(f"Completed task {processed_tasks}/{total_tasks}...")
                sys.stdout.flush()
            
            try:
                result = future.result() 
                label = result['label']
                is_error = result['error']
                
                if is_error:
                    print(f"❗ Task for '{label}' completed with error. Using zero counts.")
                    sys.stdout.flush()
                    results_counts_dict[label] = pd.Series(0.0, index=expected_cut_names)
                    results_sum_w2_dict[label] = pd.Series(0.0, index=expected_cut_names)
                else:
                    results_counts_dict[label] = result['results_counts']
                    results_sum_w2_dict[label] = result['results_sum_w2']
            
            except Exception as exc: 
                print(f'❌ Exception occurred while retrieving future result: {exc}')
                import traceback
                traceback.print_exc()
                sys.stdout.flush()

    print("Aggregating results into DataFrames...")
    sys.stdout.flush()
    cut_table_counts = pd.DataFrame(results_counts_dict)
    cut_table_sum_w2 = pd.DataFrame(results_sum_w2_dict)

    cut_table_counts = cut_table_counts.reindex(expected_cut_names, fill_value=0.0)
    cut_table_sum_w2 = cut_table_sum_w2.reindex(expected_cut_names, fill_value=0.0)

    # --- Apply Normalization to specified MC columns ---
    norm_factors = {
        "DY MC":      {"value": 2.5679e-03, "error": 2.73e-05},
        "J/Psi MC":   {"value": 1.8529e-03, "error": 2.02e-05},
        "Psi Prime MC": {"value": 2.2016e-03, "error": 7.41e-05}
    }
    mc_cols_to_normalize = list(norm_factors.keys())

    print("Applying normalization to MC columns and propagating uncertainties...")
    sys.stdout.flush()
    for mc_col in mc_cols_to_normalize:
        if mc_col in cut_table_counts.columns and mc_col in cut_table_sum_w2.columns:
            c_val = norm_factors[mc_col]["value"]
            sigma_c = norm_factors[mc_col]["error"]

            N_mc_unscaled = cut_table_counts[mc_col].copy()  # Unscaled counts (sum of weights)
            Sw_unscaled = cut_table_sum_w2[mc_col].copy()    # Variance of unscaled counts (sum of weights^2)

            # Scale the counts
            cut_table_counts[mc_col] = c_val * N_mc_unscaled

            # Propagate uncertainties for sum_w2 (variance)
            # Var(c*N) = c^2 * Var(N) + N^2 * Var(c)  (if c and N are uncorrelated)
            # Var(N) is Sw_unscaled
            # Var(c) is sigma_c^2
            term1_variance = (c_val**2) * Sw_unscaled
            term2_variance = (N_mc_unscaled**2) * (sigma_c**2)
            cut_table_sum_w2[mc_col] = term1_variance + term2_variance
            
            print(f"  Normalized '{mc_col}' with c = {c_val:.3e} ± {sigma_c:.3e}")
            sys.stdout.flush()
        else:
            print(f"⚠️ Warning: MC column '{mc_col}' for normalization not found in tables.")
            sys.stdout.flush()

    print("Recalculating 'Total MC' and derived columns (Purity, Efficiency) with normalized MC...")
    sys.stdout.flush()

    mc_cols_present = [label for label in mc_labels if label in cut_table_counts.columns]
    if mc_cols_present:
        cut_table_counts["Total MC"] = cut_table_counts[mc_cols_present].sum(axis=1, min_count=1).fillna(0)
        cut_table_sum_w2["Total MC"] = cut_table_sum_w2[mc_cols_present].sum(axis=1, min_count=1).fillna(0)
    else:
        print("⚠️ No MC columns found in results to calculate 'Total MC'. Setting to 0.")
        sys.stdout.flush()
        cut_table_counts["Total MC"] = pd.Series(0.0, index=cut_table_counts.index)
        cut_table_sum_w2["Total MC"] = pd.Series(0.0, index=cut_table_sum_w2.index) 

    try:
        dy_col = "DY MC" # This is now the normalized DY MC
        if dy_col in cut_table_counts.columns and \
           "Total Events" in cut_table_counts.index and \
           "Total MC" in cut_table_counts.columns:
            
            # total_dy_at_start is now the normalized count at "Total Events"
            total_dy_at_start = cut_table_counts.loc["Total Events", dy_col]
            if pd.isna(total_dy_at_start) or total_dy_at_start == 0:
                # This check needs to be robust if total_dy_at_start can be legitimately zero after normalization
                print(f"⚠️ Normalized total events for '{dy_col}' is zero or NaN. Efficiency will be zero or NaN.")
                sys.stdout.flush()
                total_dy_at_start_safe = np.nan 
            else:
                total_dy_at_start_safe = total_dy_at_start

            total_mc_safe = cut_table_counts["Total MC"].replace(0, np.nan) 
            purity = (cut_table_counts[dy_col] / total_mc_safe).fillna(0) * 100
            efficiency = (cut_table_counts[dy_col] / total_dy_at_start_safe).fillna(0) * 100 # Efficiency w.r.t normalized start
            
            cut_table_counts["Purity (DY MC)"] = purity
            cut_table_counts["Efficiency (DY MC)"] = efficiency
        else:
            missing_elements = []
            if dy_col not in cut_table_counts.columns: missing_elements.append(f"'{dy_col}' column")
            if "Total Events" not in cut_table_counts.index: missing_elements.append("'Total Events' index")
            if "Total MC" not in cut_table_counts.columns: missing_elements.append("'Total MC' column")
            print(f"⚠️ Cannot calculate Purity/Efficiency: Missing {', '.join(missing_elements)}.")
            sys.stdout.flush()
            cut_table_counts["Purity (DY MC)"] = 0.0
            cut_table_counts["Efficiency (DY MC)"] = 0.0
    except Exception as e:
        print(f"❌ Error calculating Purity/Efficiency: {e}")
        sys.stdout.flush()
        cut_table_counts["Purity (DY MC)"] = 0.0
        cut_table_counts["Efficiency (DY MC)"] = 0.0

    cut_table_counts.index.name = "Cut Name"
    cut_table_sum_w2.index.name = "Cut Name" 

    end_time = time.time()
    print(f"\nParallel processing finished in {end_time - start_time:.2f} seconds.")
    sys.stdout.flush()
    return cut_table_counts, cut_table_sum_w2

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
    "chisq/NDF < 12": "(chisq1/(nHits1-5)) < 12 and (chisq2/(nHits2-5)) < 12", 
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
    "abs(x1_st1 + x2_st1) < 42": "abs(x1_st1 + x2_st1) < 42", 
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
    "runID", "trackSeparation", "ReWeight", 
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
    cut_table_counts_df, cut_table_sum_w2_df = create_cut_table_parallel(
        data_file, mc_files, mc_labels_list, tree, 
        cuts_dict, variables_list,
        mixed_file, mixed_tree, empty_flask_file,
        GOOD_SPILLS_FILE_PATH
    )

    if cut_table_counts_df is not None and not cut_table_counts_df.empty and \
       cut_table_sum_w2_df is not None and not cut_table_sum_w2_df.empty:
        
        desired_order = [
            "Data", "Data(RS67)", "Mixed(RS67)",
            "empty flask (RS67)", "empty flask (RS67) mixed",
            "DY MC", "J/Psi MC", "Psi Prime MC", "Total MC",
            "Purity (DY MC)", "Efficiency (DY MC)"
        ]
        columns_to_display_ordered = [col for col in desired_order if col in cut_table_counts_df.columns]

        if not columns_to_display_ordered:
            print("❌ Error: No columns found to display after filtering. Check processing results.")
        else:
            display_df_counts = cut_table_counts_df[columns_to_display_ordered]
            cut_names_order = ["Total Events"] + list(cuts_dict.keys())
            display_df_counts = display_df_counts.reindex(cut_names_order, fill_value=0.0)

            mc_like_cols_for_error = set(mc_labels_list + ["Total MC"])
            data_like_cols_for_error = {
                "Data", "Data(RS67)", "Mixed(RS67)", 
                "empty flask (RS67)", "empty flask (RS67) mixed"
            }
            cols_for_val_pm_err_display = mc_like_cols_for_error.union(data_like_cols_for_error)

            def format_value_with_uncertainty(value_count, col_name, cut_name_idx):
                if col_name in cols_for_val_pm_err_display:
                    sum_w2_val = 0.0 
                    if cut_name_idx in cut_table_sum_w2_df.index and col_name in cut_table_sum_w2_df.columns:
                         sum_w2_val = cut_table_sum_w2_df.loc[cut_name_idx, col_name]
                    else:
                        print(f"⚠️ Warning: Value for uncertainty calculation not found for '{col_name}' at cut '{cut_name_idx}'. Using 0 uncertainty.")
                        sys.stdout.flush()


                    if pd.isna(value_count): value_count = 0.0 
                    if pd.isna(sum_w2_val): sum_w2_val = 0.0
                    
                    uncertainty = np.sqrt(sum_w2_val) if sum_w2_val >= 0 else 0.0
                    
                    # For MC columns (now normalized), and Total MC, typically show some decimals
                    # For Data columns, counts are integers usually, but error makes it float display
                    if col_name in mc_like_cols_for_error: 
                        return "{:.1f} \u00B1 {:.1f}".format(value_count, uncertainty)
                    else: # Data-like columns (counts are integers, error can be float)
                        return "{:.0f} \u00B1 {:.1f}".format(value_count, uncertainty) # error with 1 decimal for data too
                
                elif col_name in ["Purity (DY MC)", "Efficiency (DY MC)"]:
                    try:
                        return "{:.2f}%".format(float(value_count))
                    except (ValueError, TypeError): return str(value_count) 
                else: 
                    return str(value_count)

            headers_for_tabulate = ["Cut"] + columns_to_display_ordered
            formatted_data_for_tabulate = []

            for cut_idx_name in display_df_counts.index: 
                row_values = [cut_idx_name] 
                for col_name in columns_to_display_ordered: 
                    raw_count_value = display_df_counts.loc[cut_idx_name, col_name]
                    row_values.append(format_value_with_uncertainty(raw_count_value, col_name, cut_idx_name))
                formatted_data_for_tabulate.append(row_values)
            
            print("\n--- Final Cut Table (with Uncertainties) ---")
            sys.stdout.flush()
            print(tabulate(formatted_data_for_tabulate, headers=headers_for_tabulate, tablefmt="fancy_grid"))
            sys.stdout.flush()

            try:
                output_csv_counts = "cut_table_normalized_counts_conditional_goodspills.csv"
                cut_table_counts_df.to_csv(output_csv_counts) 
                print(f"\n✅ Normalized counts table saved as '{output_csv_counts}'")

                output_csv_sum_w2 = "cut_table_propagated_sum_w2_conditional_goodspills.csv"
                cut_table_sum_w2_df.to_csv(output_csv_sum_w2) 
                print(f"✅ Propagated sum of weights squared (variances) table saved as '{output_csv_sum_w2}'")
            except Exception as e:
                print(f"\n❌ Error saving CSV files: {e}")
            sys.stdout.flush()
    else:
        print("❌ Cut table DataFrames were not generated or are empty. Cannot display or save.")
        sys.stdout.flush()