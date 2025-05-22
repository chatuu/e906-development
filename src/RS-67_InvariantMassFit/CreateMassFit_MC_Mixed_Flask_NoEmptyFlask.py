# -*- coding: utf-8 -*-
"""
Process High Energy Physics data from ROOT files to generate a cut flow table
and perform fits using SciPy.
"""

import uproot
import pandas as pd
# from tabulate import tabulate # No longer needed
import numpy as np
import concurrent.futures
import time
import re
import sys
import csv # Added for CSV output, though pandas can also do it

# SciPy and Matplotlib for fitting and plotting
from scipy.optimize import minimize, approx_fprime
import matplotlib.pyplot as plt


# Define target labels for the Good Spills Cut globally for clarity
TARGET_LABELS_FOR_GOOD_SPILLS_CUT = {"Data", "Data(RS67)", "empty flask (RS67)"}

# Define the constant factor for flask subtraction/scaling (No longer used in the new fit model, but kept for now if other parts rely on it)
# FLASK_FACTOR = 4.39933159


def load_good_spills(file_path):
    good_spills = set()
    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                stripped_line = line.strip()
                if stripped_line:
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
        return set()
    except Exception as e:
        print(f"❌ Error loading good spills file '{file_path}': {e}")
        sys.stdout.flush()
        return set()

def read_tree(file_path, tree_name, variables):
    try:
        with uproot.open(file_path) as file:
            if tree_name not in file:
                print(f"❌ Error: Tree '{tree_name}' not found in file '{file_path}'.")
                sys.stdout.flush()
                empty_df = pd.DataFrame()
                all_expected_vars = set(variables + ['runID', 'dy', 'ReWeight', 'spillID', 'mass'])
                for var in all_expected_vars:
                    empty_df[var] = pd.Series(dtype='float64')
                return empty_df

            tree = file[tree_name]
            all_keys = tree.keys()
            vars_to_request = list(set(variables + ['runID', 'dy', 'ReWeight', 'spillID', 'mass']))
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
                if var not in df.columns: df[var] = np.nan
            if 'runID' in df.columns: df['runID'] = pd.to_numeric(df['runID'], errors='coerce')
            if 'dy' in df.columns: df['dy'] = pd.to_numeric(df['dy'], errors='coerce')
            if 'spillID' in df.columns: df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')
            if 'mass' in df.columns: df['mass'] = pd.to_numeric(df['mass'], errors='coerce')
            if 'ReWeight' in df.columns: df['ReWeight'] = pd.to_numeric(df['ReWeight'], errors='coerce').fillna(1.0)
            else: df['ReWeight'] = 1.0
            return df
    except Exception as e:
        print(f"❌ Error reading tree '{tree_name}' from '{file_path}': {e}")
        sys.stdout.flush()
        empty_df = pd.DataFrame()
        all_expected_vars = set(variables + ['runID', 'dy', 'ReWeight', 'spillID', 'mass'])
        for var in all_expected_vars: empty_df[var] = pd.Series(dtype='float64')
        return empty_df

def add_beam_offset(df):
    df['beamOffset'] = 1.6
    return df

def apply_cuts(df, cuts, use_weights=True, good_spills_set=None, dataset_label=None, return_final_state=False):
    cut_results = {}
    weight_series = df['ReWeight'] if use_weights and 'ReWeight' in df.columns else pd.Series(1.0, index=df.index)
    weight_series = pd.to_numeric(weight_series, errors='coerce').fillna(1.0)
    cut_results["Total Events"] = weight_series.sum()

    if 'beamOffset' not in df.columns: df['beamOffset'] = 1.6
    if 'dy' in df.columns: df['dy'] = pd.to_numeric(df['dy'], errors='coerce').fillna(np.nan)
    else: df['dy'] = np.nan
    if 'runID' in df.columns: df['runID'] = pd.to_numeric(df['runID'], errors='coerce').fillna(-1)
    else: df['runID'] = -1
    if 'spillID' in df.columns: df['spillID'] = pd.to_numeric(df['spillID'], errors='coerce')
    else: df['spillID'] = np.nan

    cumulative_mask = pd.Series(True, index=df.index)
    original_df_for_cuts = df
    original_weights_for_cuts = weight_series.copy()

    for cut_name, cut_string in cuts.items():
        if cut_name == "Total Events": continue
        temp_cumulative_mask_before_this_cut = cumulative_mask.copy()
        try:
            current_df_subset = original_df_for_cuts[cumulative_mask]
            if current_df_subset.empty:
                mask_this_cut_on_subset = pd.Series(False, index=current_df_subset.index)
            elif cut_name == "Good Spills Cut":
                if dataset_label in TARGET_LABELS_FOR_GOOD_SPILLS_CUT:
                    if good_spills_set is not None and 'spillID' in current_df_subset.columns and not good_spills_set == set() and not current_df_subset['spillID'].isna().all():
                        current_spillID_numeric = pd.to_numeric(current_df_subset['spillID'], errors='coerce')
                        mask_this_cut_on_subset = current_spillID_numeric.isin(good_spills_set) & current_spillID_numeric.notna()
                    else:
                        if not good_spills_set or good_spills_set == set(): print(f"⚠️ Warning: 'Good Spills Cut' for '{dataset_label}' ineffective (spill list empty/not loaded). Passing subset.")
                        if 'spillID' not in current_df_subset.columns or current_df_subset['spillID'].isna().all(): print(f"⚠️ Warning: 'Good Spills Cut' for '{dataset_label}' cannot be applied ('spillID' missing or all NaN). Passing subset.")
                        sys.stdout.flush()
                        mask_this_cut_on_subset = pd.Series(True, index=current_df_subset.index)
                else:
                    mask_this_cut_on_subset = pd.Series(True, index=current_df_subset.index)
            else:
                mask_this_cut_on_subset = current_df_subset.eval(cut_string, engine='python')

            aligned_mask_for_this_cut_step = pd.Series(False, index=original_df_for_cuts.index)
            if not current_df_subset.empty:
                aligned_mask_for_this_cut_step.loc[current_df_subset[mask_this_cut_on_subset].index] = True

            cumulative_mask = cumulative_mask & aligned_mask_for_this_cut_step
            cut_results[cut_name] = original_weights_for_cuts[cumulative_mask].sum()
        except Exception as e:
            print(f"⚠️ Error applying cut '{cut_name}' for dataset '{dataset_label}': {e}. This cut results in 0 for this step.")
            sys.stdout.flush()
            cut_results[cut_name] = 0.0
            cumulative_mask = temp_cumulative_mask_before_this_cut
    if return_final_state:
        df_passed_all_cuts = original_df_for_cuts[cumulative_mask]
        weights_passed_all_cuts = original_weights_for_cuts[cumulative_mask]
        return pd.Series(cut_results), df_passed_all_cuts, weights_passed_all_cuts
    else:
        return pd.Series(cut_results)

def process_single_file(file_path, tree_name, label, variables, cuts, use_weights, good_spills_set, is_fit_component=False):
    try:
        print(f"Starting processing for {label} from {file_path} (Tree: {tree_name})...")
        sys.stdout.flush()
        vars_to_read_list = list(set(variables + ['runID', 'spillID', 'mass']))
        if use_weights: vars_to_read_list.append('ReWeight')

        df = read_tree(file_path, tree_name, vars_to_read_list)

        if df.empty or ('mass' not in df.columns and is_fit_component):
            print(f"INFO: read_tree returned unusable DataFrame for {label} (empty or missing 'mass' for fit component). Skipping.")
            sys.stdout.flush()
            expected_cut_names = ["Total Events"] + list(cuts.keys())
            cut_series_dummy = pd.Series(0.0, index=expected_cut_names)
            return {'label': label, 'results': cut_series_dummy, 'error': False,
                    'mass_data': np.array([]), 'weights_data': np.array([])}

        df = add_beam_offset(df)

        cut_variables_in_strings = set()
        for cut_string in cuts.values():
            if "#" in cut_string:
                custom_part = cut_string.split("#")[0].strip()
                if not custom_part: continue
                current_eval_string = custom_part
            else:
                current_eval_string = cut_string
            found_vars = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', current_eval_string)
            keywords_to_exclude = {'and', 'or', 'not', 'in', 'is', 'True', 'False', 'None', 'abs', 'np', 'sqrt', 'min', 'max'}
            valid_vars = [var for var in found_vars if var not in keywords_to_exclude and not re.match(r'^[0-9]+(\.[0-9]*)?$', var)]
            cut_variables_in_strings.update(valid_vars)
        if "Good Spills Cut" in cuts: cut_variables_in_strings.add('spillID')

        all_required_vars_for_df = list(set(variables_list + list(cut_variables_in_strings) + ['beamOffset', 'dy', 'runID', 'spillID', 'mass', 'ReWeight']))
        for var in all_required_vars_for_df:
            if var not in df.columns: df[var] = np.nan

        cols_to_numeric = ['mass', 'dy', 'runID', 'spillID', 'ReWeight']
        for col_num in cols_to_numeric:
            if col_num in df.columns: df[col_num] = pd.to_numeric(df[col_num], errors='coerce')

        if 'ReWeight' not in df.columns:
            df['ReWeight'] = 1.0
        df['ReWeight'] = df['ReWeight'].fillna(1.0)

        mass_data_final, weights_data_final = np.array([]), np.array([])
        cut_series, df_passed_all, weights_series_passed_all = apply_cuts(df, cuts, use_weights=use_weights, good_spills_set=good_spills_set, dataset_label=label, return_final_state=True)

        if is_fit_component:
            if 'mass' in df_passed_all.columns and not df_passed_all.empty:
                mass_values_passed = df_passed_all['mass']
                valid_mass_mask = mass_values_passed.notna()
                mass_data_final = mass_values_passed[valid_mass_mask].to_numpy()

                if weights_series_passed_all is not None:
                    weights_data_final = weights_series_passed_all[valid_mass_mask].to_numpy()
                else:
                    weights_data_final = np.ones_like(mass_data_final)

        print(f"Finished processing {label}.")
        sys.stdout.flush()
        return {'label': label, 'results': cut_series, 'error': False, 'mass_data': mass_data_final, 'weights_data': weights_data_final}

    except Exception as e:
        print(f"❌ Error processing {label} from {file_path} (Tree: {tree_name}): {e}")
        sys.stdout.flush()
        expected_cut_names = ["Total Events"] + list(cuts.keys())
        return {'label': label, 'results': pd.Series(0.0, index=expected_cut_names), 'error': True, 'message': str(e), 'mass_data': None, 'weights_data': None}

def prepare_fit_data_parallel(data_file_path, mc_file_paths, mc_labels, tree_name, current_cuts_dict, variables_param,
                              mixed_file_path, mixed_tree_name, empty_flask_file_path, good_spills_file_path):
    start_time = time.time()
    print("Starting parallel processing for FIT DATA...")
    sys.stdout.flush()
    good_spills_set = load_good_spills(good_spills_file_path)
    fit_data_arrays, futures = {}, []

    all_datasets_for_fit_construction = {
        "Data(RS67)", "Mixed(RS67)", "empty flask (RS67)",
        "empty flask (RS67) mixed",
        "DY MC", "J/Psi MC", "Psi Prime MC"
    }
    tasks = [
        (data_file_path, tree_name, "Data", variables_param, current_cuts_dict, False, good_spills_set, False),
        (mixed_file_path, "result", "Data(RS67)", variables_param, current_cuts_dict, False, good_spills_set, True),
        (mixed_file_path, mixed_tree_name, "Mixed(RS67)", variables_param, current_cuts_dict, False, good_spills_set, True),
        (empty_flask_file_path, "result", "empty flask (RS67)", variables_param, current_cuts_dict, False, good_spills_set, True),
        (empty_flask_file_path, "result_mix", "empty flask (RS67) mixed", variables_param, current_cuts_dict, False, good_spills_set, True)
    ]
    for i, mc_file_path in enumerate(mc_file_paths):
        label = mc_labels[i]
        mc_vars = list(set(variables_param + ["ReWeight"]))
        tasks.append((mc_file_path, tree_name, label, mc_vars, current_cuts_dict, True, good_spills_set, True))

    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        for task_args in tasks: futures.append(executor.submit(process_single_file, *task_args))
        print(f"Submitted {len(futures)} tasks to executor...")
        sys.stdout.flush()
        processed_tasks, total_tasks = 0, len(futures)
        for future in concurrent.futures.as_completed(futures):
            processed_tasks += 1
            if total_tasks < 20 or processed_tasks % max(1, total_tasks // 10) == 0 or processed_tasks == total_tasks:
                print(f"Completed task {processed_tasks}/{total_tasks}...")
                sys.stdout.flush()
            try:
                result = future.result()
                label, is_error = result['label'], result['error']
                if is_error:
                    error_message = result.get('message', 'Unknown error during processing')
                    print(f"❗ Task for '{label}' completed with error: {error_message}.")
                    sys.stdout.flush()
                    if label in all_datasets_for_fit_construction:
                        fit_data_arrays[label] = {'mass': np.array([]), 'weights': np.array([])}
                else:
                    mass_data = result['mass_data']
                    weights_data = result['weights_data']
                    fit_data_arrays[label] = {
                        'mass': mass_data if mass_data is not None else np.array([]),
                        'weights': weights_data if weights_data is not None else np.array([])
                    }
            except Exception as exc:
                print(f'❌ Exception occurred while retrieving future result: {exc}. Affected task details might be lost.')
                sys.stdout.flush()

    print(f"\nParallel processing for fit data finished in {time.time() - start_time:.2f} seconds.")
    sys.stdout.flush()
    return fit_data_arrays

def perform_scipy_fit(fit_data_arrays, mass_min=3.0, mass_max=7.0, n_bins=100,
                      fit_range_gev_min=None, fit_range_gev_max=None,
                      parameter_constraints=None, use_logy_plot=False,
                      initial_param_guesses=None):
    """
    Performs a binned template fit to Data(RS67) histogram using SciPy.
    Target: Data(RS67)
    Model: Total = c_mixed*Mixed(RS67) + c_empty*EF(RS67)
                 + c_dy*DY_MC + c_jpsi*JPsi_MC + c_psip*PsiP_MC
    Attempts to numerically calculate Hessian for error estimation.
    """
    print("\n--- Starting SciPy Template Fit (Direct Data Model Modified) ---")
    print(f"Histograms: {n_bins} bins from {mass_min:.2f} to {mass_max:.2f} GeV.")

    bin_edges = np.linspace(mass_min, mass_max, n_bins + 1) # Full range bin edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2     # Full range bin centers

    all_component_labels_needed = ["Data(RS67)", "Mixed(RS67)", "empty flask (RS67)",
                                   "empty flask (RS67) mixed",
                                   "DY MC", "J/Psi MC", "Psi Prime MC"]

    component_contents = {} # To store full-range histograms
    component_raw_counts_for_error_calc = {} # Full-range raw counts for target

    template_definitions_ordered_new_model = [
        ("Mixed(RS67)", "c_mixed"),
        ("empty flask (RS67)", "c_empty"),
        ("DY MC", "c_dy"),
        ("J/Psi MC", "c_jpsi"),
        ("Psi Prime MC", "c_psip")
    ]

    for label in all_component_labels_needed:
        if label not in fit_data_arrays or fit_data_arrays[label]['mass'] is None:
            print(f"⚠️ Missing mass data for component '{label}'.")
            is_essential_for_this_fit = label == "Data(RS67)" or label in [tpl[0] for tpl in template_definitions_ordered_new_model]
            if is_essential_for_this_fit:
                print(f"❌ Component '{label}' is essential for the current fit and missing. Cannot perform fit.")
                return None
            component_contents[label] = np.zeros(n_bins, dtype=float)
            if label == "Data(RS67)": component_raw_counts_for_error_calc[label] = np.zeros(n_bins, dtype=float)
            continue

        mass_data = fit_data_arrays[label]['mass']
        weights_data = fit_data_arrays[label]['weights']
        current_hist_contents = np.zeros(n_bins, dtype=float)
        if isinstance(mass_data, np.ndarray) and mass_data.size > 0:
            current_hist_contents, _ = np.histogram(mass_data, bins=bin_edges, weights=weights_data) # Histogram over full range
            if label == "Data(RS67)":
                raw_counts_target, _ = np.histogram(mass_data, bins=bin_edges) # Raw counts over full range
                component_raw_counts_for_error_calc[label] = raw_counts_target.astype(float)
        elif isinstance(mass_data, np.ndarray) and mass_data.size == 0:
            print(f"ℹ️ No mass data entries for component '{label}'. Assuming zero content.")
        else:
            print(f"ℹ️ Mass data for '{label}' is not a numpy array or None. Assuming zero content.")
        component_contents[label] = current_hist_contents.astype(float) # Store full-range histogram
        if label == "Data(RS67)" and label not in component_raw_counts_for_error_calc:
             component_raw_counts_for_error_calc[label] = np.zeros(n_bins, dtype=float)
        print(f"Prepared component '{label}': {np.sum(current_hist_contents):.2f} effective entries (sum of weights/counts).")

    target_data_label = "Data(RS67)"
    if target_data_label not in component_contents or target_data_label not in component_raw_counts_for_error_calc:
        print(f"❌ Target data '{target_data_label}' not prepared correctly. Aborting fit.")
        return None

    target_data_contents = component_contents[target_data_label] # Full-range target data
    target_data_errors = np.sqrt(np.maximum(1.0, component_raw_counts_for_error_calc[target_data_label])) # Full-range errors
    target_data_errors[target_data_errors < 1e-9] = 1.0
    mask_zeroed_content_and_small_error = (target_data_contents == 0) & (target_data_errors < 1.0)
    target_data_errors[mask_zeroed_content_and_small_error] = 1.0

    print(f"Constructed Target Data ('{target_data_label}') histogram with integral: {np.sum(target_data_contents):.2f}")
    if np.sum(target_data_contents) <= 1e-9 and np.all(target_data_contents <=1e-9):
        print("❌ Target Data histogram has zero or negligible integral. Aborting fit.")
        return None

    templates_hist_content = {} # To store full-range template histograms selected for the fit model
    active_templates_info = []
    min_sum_weights_for_template = 1.0

    for source_label, param_name_for_fit in template_definitions_ordered_new_model:
        if source_label not in component_contents:
            print(f"⚠️ Missing data for template component '{source_label}'. Excluding from fit.")
            # templates_hist_content[source_label] = np.zeros(n_bins, dtype=float) # Not strictly needed if excluded
            continue
        current_template_hist = component_contents[source_label] # This is the full-range histogram
        current_template_hist[current_template_hist < 0] = 0.0
        sum_of_hist_content = np.sum(current_template_hist)
        is_mc_component = "MC" in source_label
        print_label_type = "MC TEMPLATE" if is_mc_component else "BACKGROUND TEMPLATE"
        print(f"Prepared {print_label_type} '{source_label}': Sum of (weighted) entries = {sum_of_hist_content:.2f}.")
        if sum_of_hist_content >= min_sum_weights_for_template:
            templates_hist_content[source_label] = current_template_hist # Store the full-range histogram
            active_templates_info.append((source_label, param_name_for_fit))
        else:
            # templates_hist_content[source_label] = np.zeros(n_bins, dtype=float) # Not strictly needed if excluded
            print(f"⚠️ EXCLUDING template '{source_label}' (Sum {sum_of_hist_content:.2f} < {min_sum_weights_for_template}).")

    if not active_templates_info:
        print("❌ No active templates with sufficient statistics. Cannot perform SciPy fit.")
        return None

    active_param_names_ordered = [info[1] for info in active_templates_info]
    print(f"ℹ️ Active parameters for SciPy fit: {active_param_names_ordered}")

    fit_idx_start, fit_idx_end = 0, n_bins
    if fit_range_gev_min is not None and fit_range_gev_max is not None and fit_range_gev_max > fit_range_gev_min:
        actual_bin_idx_start = np.searchsorted(bin_edges, fit_range_gev_min, side='left')
        actual_bin_idx_end = np.searchsorted(bin_edges, fit_range_gev_max, side='right')
        if actual_bin_idx_start < actual_bin_idx_end and actual_bin_idx_start < n_bins and actual_bin_idx_end > 0 :
            fit_idx_start = actual_bin_idx_start
            fit_idx_end = min(actual_bin_idx_end, n_bins)
            print(f"ℹ️ SciPy fit range: bins from index {fit_idx_start} to {fit_idx_end-1} "
                  f"(approx {bin_edges[fit_idx_start]:.2f} - {bin_edges[fit_idx_end]:.2f} GeV).")
        else:
            print(f"⚠️ Invalid SciPy fit range definition. Using full histogram range for fit.")
            fit_idx_start, fit_idx_end = 0, n_bins
    else:
        print("ℹ️ SciPy fit using full histogram range (no specific sub-range defined).")

    # Data for the FIT (sliced to fit region)
    target_data_fit = target_data_contents[fit_idx_start:fit_idx_end]
    target_errors_fit = target_data_errors[fit_idx_start:fit_idx_end]
    # Templates for the FIT (sliced to fit region)
    templates_fit = {label_key: templates_hist_content[label_key][fit_idx_start:fit_idx_end]
                     for label_key, _ in active_templates_info}
    target_errors_fit[target_errors_fit <= 1e-9] = 1.0

    # Model function for the FIT (operates on fit-region data)
    def model_total_func_fit_region(params_array):
        prediction = np.zeros_like(target_data_fit, dtype=float)
        for i, (active_label, _) in enumerate(active_templates_info):
            prediction += params_array[i] * templates_fit[active_label] # Uses sliced templates_fit
        return prediction

    def chi_squared_total_func(params_array):
        model = model_total_func_fit_region(params_array) # Uses fit-region model
        residuals = target_data_fit - model
        chi2_terms = (residuals**2) / (target_errors_fit**2)
        return np.sum(chi2_terms)

    initial_params_list = []
    bounds_list = []
    default_guess = 0.1
    default_bounds = (0.0, None)

    for p_name_key in active_param_names_ordered:
        guess = (initial_param_guesses.get(p_name_key, default_guess) if initial_param_guesses else default_guess)
        low_b, high_b = (parameter_constraints.get(p_name_key, default_bounds) if parameter_constraints else default_bounds)
        initial_params_list.append(guess)
        bounds_list.append((low_b, high_b))

    print(f"ℹ️ SciPy initial parameter values: {initial_params_list}")
    print(f"ℹ️ SciPy parameter bounds: {bounds_list}")

    if not initial_params_list:
        print("❌ No active parameters to fit. Aborting SciPy fit.")
        return None

    print("Performing SciPy minimization for Total Model fit...")
    minimizer_options = {'maxiter': 20000, 'ftol': 1e-9, 'disp': False}
    result = minimize(chi_squared_total_func, initial_params_list, method='L-BFGS-B', bounds=bounds_list, options=minimizer_options)

    fit_params = result.x
    fit_success = result.success
    chi2_val = result.fun
    num_bins_in_fit = len(target_data_fit)
    ndf_reduction = len(fit_params)
    ndf = num_bins_in_fit - ndf_reduction if num_bins_in_fit > ndf_reduction else 0
    param_errors = np.full_like(fit_params, np.nan)
    err_chi2_ndf = np.nan

    params_at_bounds_info = []
    for i, (_, p_name_key) in enumerate(active_templates_info):
        val = fit_params[i]
        low_b, high_b = bounds_list[i]
        if low_b is not None and np.isclose(val, low_b, rtol=1e-5, atol=1e-8):
            params_at_bounds_info.append(f"{p_name_key} at lower bound ({low_b})")
        if high_b is not None and np.isclose(val, high_b, rtol=1e-5, atol=1e-8):
            params_at_bounds_info.append(f"{p_name_key} at upper bound ({high_b})")
    if params_at_bounds_info:
        print("⚠️ WARNING: Some fit parameters are at their bounds, error estimation might be unreliable:")
        for info in params_at_bounds_info: print(f"   - {info}")

    if fit_success:
        print("✅ SciPy Fit (Total Model): Minimization successful.")
        # Gradient/Hessian checks remain for the fit-region chi2 function
        current_gradient_at_solution = approx_fprime(fit_params, chi_squared_total_func, np.sqrt(np.finfo(float).eps) * np.maximum(1.0, np.abs(fit_params)))
        print(f"ℹ️ Gradient at solution (approx): {current_gradient_at_solution}")
        print(f"ℹ️ Fit parameters at solution: {fit_params}")
        print("Attempting numerical calculation of Hessian for error estimation...")
        epsilon_for_grad_func = np.sqrt(np.finfo(float).eps) * np.maximum(1.0, np.abs(fit_params))
        def grad_chi2_at_minimum(p_array_for_grad): # This is for the fit-region chi2
            return approx_fprime(p_array_for_grad, chi_squared_total_func, epsilon_for_grad_func)
        default_epsilon_for_hess_step = np.sqrt(np.finfo(float).eps)
        numerical_hessian = np.zeros((len(fit_params), len(fit_params)))
        try:
            for i in range(len(fit_params)):
                def single_grad_component_func(p_array_for_hess_row, component_idx_i):
                    grad_vector = grad_chi2_at_minimum(p_array_for_hess_row)
                    return grad_vector[component_idx_i]
                hessian_row_i = approx_fprime(fit_params, single_grad_component_func,
                                              default_epsilon_for_hess_step * np.maximum(1.0, np.abs(fit_params[i])), i)
                numerical_hessian[i, :] = hessian_row_i
            numerical_hessian = (numerical_hessian + numerical_hessian.T) / 2.0
            eigenvalues, eigenvectors = np.linalg.eigh(numerical_hessian)
            print(f"ℹ️ Original eigenvalues of numerical Hessian: {eigenvalues}")
            min_eigenvalue_allowed = 1e-6
            if np.any(eigenvalues <= min_eigenvalue_allowed / 10):
                print(f"⚠️ Numerical Hessian is not robustly positive definite. Original eigenvalues: {eigenvalues}")
                print(f"   Attempting to correct Hessian by flooring eigenvalues to >= {min_eigenvalue_allowed}.")
                corrected_eigenvalues = np.maximum(eigenvalues, min_eigenvalue_allowed)
                if not np.array_equal(eigenvalues, corrected_eigenvalues): print(f"   Corrected eigenvalues: {corrected_eigenvalues}")
                numerical_hessian_corrected = eigenvectors @ np.diag(corrected_eigenvalues) @ eigenvectors.T
                try:
                    covariance_matrix = np.linalg.inv(numerical_hessian_corrected)
                    diag_cov = np.diag(covariance_matrix)
                    param_errors = np.sqrt(np.abs(diag_cov))
                    if np.any(diag_cov < 0): print("⚠️ Warning: Some diagonal elements of the 'corrected' covariance matrix were negative. abs() taken for sqrt.")
                    print(f"ℹ️ SciPy parameter errors estimated from 'corrected' numerical Hessian: {param_errors}")
                except np.linalg.LinAlgError as lae_corr:
                    print(f"⚠️ Linear algebra error during 'corrected' Hessian inversion: {lae_corr}. Errors remain NaN.")
                    param_errors = np.full_like(fit_params, np.nan)
            else:
                covariance_matrix = np.linalg.inv(numerical_hessian)
                param_errors = np.sqrt(np.diag(covariance_matrix))
                print(f"ℹ️ SciPy parameter errors estimated from numerical Hessian: {param_errors}")
        except np.linalg.LinAlgError as lae:
            print(f"⚠️ Linear algebra error during Hessian processing: {lae}. Errors set to NaN.")
            param_errors = np.full_like(fit_params, np.nan)
        except Exception as e:
            print(f"⚠️ Error during numerical Hessian calculation: {e}. Errors set to NaN.")
            param_errors = np.full_like(fit_params, np.nan)
    else:
        print(f"❌ SciPy Fit (Total Model) failed. Status: {result.status}, Message: {result.message}")

    if ndf > 0: err_chi2_ndf = np.sqrt(2.0 / ndf)

    fit_results_output = {"chi2": chi2_val, "ndf": ndf, "status": result.status,
                          "message": result.message, "success": fit_success, "err_chi2_ndf": err_chi2_ndf}
    for i, (_, param_name_key) in enumerate(active_templates_info):
        fit_results_output[param_name_key] = fit_params[i]
        fit_results_output[f"err_{param_name_key}"] = param_errors[i]

    # --- Plotting ---
    plt.figure(figsize=(12, 8))
    # Plot target data points over the FULL histogram range
    plt.errorbar(bin_centers, target_data_contents, yerr=target_data_errors, fmt='o', label=f'Target: {target_data_label}', color='black', capsize=3, elinewidth=1, markeredgewidth=1, ms=5, zorder=5)

    if fit_success:
        # Construct Total Model for PLOTTING using FULL-RANGE templates and fitted parameters
        total_model_plot_full_range = np.zeros_like(bin_centers, dtype=float) # bin_centers is full range
        for i, (active_label_for_param, _) in enumerate(active_templates_info):
            if active_label_for_param in templates_hist_content:
                full_range_template = templates_hist_content[active_label_for_param]
                total_model_plot_full_range += fit_params[i] * full_range_template
            else:
                print(f"⚠️ Plotting Warning: Active template '{active_label_for_param}' not found in templates_hist_content.")

        # Plot the full-range total model line
        try:
            plt.stairs(total_model_plot_full_range, bin_edges, label='Total Fit Model', color='red', linewidth=2, zorder=4, baseline=None)
        except AttributeError:
            plt.step(bin_edges[:-1], total_model_plot_full_range, where='post', label='Total Fit Model', color='red', linewidth=2, zorder=4)

        # Plot Individual Scaled Components for PLOTTING using FULL-RANGE templates
        component_colors = ['green', 'blue', 'purple', 'orange', 'cyan', 'magenta', 'brown', 'pink']
        for i, (label, param_name_key) in enumerate(active_templates_info):
            if label in templates_hist_content and i < len(fit_params):
                full_range_template_component = templates_hist_content[label] # Use full-range template
                scaled_component_plot_full_range = fit_params[i] * full_range_template_component

                if np.sum(np.abs(scaled_component_plot_full_range)) > 1e-6:
                    try:
                        plt.stairs(scaled_component_plot_full_range, bin_edges, label=f"{label} ({param_name_key}={fit_params[i]:.2e})",
                                   linestyle='-', color=component_colors[i % len(component_colors)], linewidth=1.5, baseline=None, zorder=3)
                    except AttributeError:
                        plt.step(bin_edges[:-1], scaled_component_plot_full_range, where='post', label=f"{label} ({param_name_key}={fit_params[i]:.2e})",
                                   linestyle='-', color=component_colors[i % len(component_colors)], linewidth=1.5, zorder=3)
            else:
                 print(f"⚠️ Plotting Warning: Template '{label}' for param '{param_name_key}' not found or param index issue.")


    plt.xlabel("Mass (GeV)", fontsize=12)
    bin_width_gev = (mass_max - mass_min) / n_bins # Uses full range mass_min, mass_max
    plt.ylabel(f"Events / ({bin_width_gev:.3f} GeV)", fontsize=12)

    ndf_val = fit_results_output['ndf']
    chi2_val_title = fit_results_output['chi2']
    chi2_ndf_val_str = f"{chi2_val_title/ndf_val:.2f}" if ndf_val > 0 else "N/A"
    err_chi2_ndf_val = fit_results_output.get('err_chi2_ndf', np.nan)
    chi2_ndf_err_str = f"{err_chi2_ndf_val:.2f}" if ndf_val > 0 and not np.isnan(err_chi2_ndf_val) else "N/A"

    if ndf_val > 0 and chi2_ndf_val_str != "N/A" and chi2_ndf_err_str != "N/A":
        title_chi2_part = rf"$\chi^2$/NDF = {chi2_val_title:.2f}/{ndf_val} = {chi2_ndf_val_str} $\pm$ {chi2_ndf_err_str}"
    else:
        title_chi2_part = rf"$\chi^2$/NDF = {chi2_val_title:.2f}/{ndf_val} = {chi2_ndf_val_str}"

    plt.title(f"SciPy Fit to {target_data_label} (Status: {'Success' if fit_success else 'Failed'})\n{title_chi2_part}", fontsize=14)

    if use_logy_plot:
        plt.yscale('log'); plt.ylim(bottom=max(0.1, plt.ylim()[0] if plt.ylim()[0] > 0 else 0.1 ))
    else:
        current_ylim_min = plt.ylim()[0]
        if current_ylim_min < 0 and np.any(target_data_contents >= 0):
            positive_data = target_data_contents[target_data_contents > 0]
            if len(positive_data) > 0 and np.mean(positive_data) > -2*current_ylim_min:
                 plt.ylim(bottom=0)

    plt.legend(loc='upper right', fontsize='small'); plt.grid(True, linestyle=':', alpha=0.6); plt.tight_layout()
    plot_filename_base = "scipy_fit_direct_data_model_modified"
    try:
        plt.savefig(plot_filename_base + ".png"); plt.savefig(plot_filename_base + ".pdf")
        print(f"ℹ️ SciPy fit plot saved as {plot_filename_base}.png/pdf")
    except Exception as e: print(f"⚠️ Error saving SciPy fit plot: {e}")
    plt.close()
    return fit_results_output

# ------------------- Configuration -------------------
GOOD_SPILLS_FILE_PATH = "../../res/GoodSpills/GoodSpills.txt"
full_cuts_dict = {
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
    "D1 < 400": "D1 < 400", "D2 < 400": "D2 < 400", "D3 < 400": "D3 < 400",
    "D1 + D2 + D3 < 1000": "D1 + D2 + D3 < 1000",
}
variables_list = ["mass", "spillID", "matrix1", "nHits1", "nHits2", "nHits1St1", "nHits2St1",
                  "chisq1", "chisq2", "chisq_dimuon", "chisq1_target", "chisq2_target",
                  "dx", "dy", "dz", "dpx", "dpy", "dpz", "D1", "D2", "D3",
                  "xF", "xT", "xB", "costh", "intensityP", "runID", "trackSeparation",
                  "x1_t", "x2_t", "y1_t", "y2_t", "x1_d", "x2_d", "y1_d", "y2_d",
                  "z1_v", "z2_v", "x1_st1", "x2_st1", "y1_st1", "y2_st1", "y1_st3", "y2_st3",
                  "py1_st1", "py2_st1", "py1_st3", "py2_st3", "pz1_st1", "pz2_st1", "pz1_st3", "pz2_st3",
                  "px1_st1", "px2_st1", "px1_st3", "px2_st3",
                  "ReWeight", "beamOffset"
                 ]
variables_list = list(set(variables_list))

try:
    base_path_hugo = "../../res/ROOTFiles/Hugo/"
    base_path_mixed = "../../res/ROOTFiles/MixedEvents/"
    data_file = base_path_hugo + "roadset57_70_R008_2111v42_tmp_noPhys.root"
    mc_files = [ base_path_hugo + "mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
                 base_path_hugo + "mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
                 base_path_hugo + "mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root" ]
    mixed_file = base_path_mixed + "merged_RS67_3089LH2.root"
    empty_flask_file = base_path_mixed + "merged_RS67_3089flask.root"
except NameError:
    print("⚠️ Warning: Using placeholder file paths due to undefined base paths. Ensure paths are correct.")
    sys.stdout.flush()
    data_file, mc_files, mixed_file, empty_flask_file = "data.root", ["drellyan.root", "jpsi.root", "psiprime.root"], "mixed.root", "empty_flask.root"

mc_labels_list, tree, mixed_tree = ["DY MC", "J/Psi MC", "Psi Prime MC"], "Tree", "result_mix"

def format_param_err_str(value, error, val_fmt=".4e", err_fmt=".2e"):
    if error is None or np.isnan(error) or not np.isfinite(error):
        return f"{value:{val_fmt}} +/- nan"
    else:
        return f"{value:{val_fmt}} +/- {error:{err_fmt}}"

if __name__ == "__main__":
    active_cuts_dict, target_cut_name_to_stop_after, found_target_cut = {}, "D1 + D2 + D3 < 1000", False
    for cut_name, cut_string in full_cuts_dict.items():
        active_cuts_dict[cut_name] = cut_string
        if cut_name == target_cut_name_to_stop_after:
            found_target_cut = True; break
    if not found_target_cut:
        print(f"⚠️ FATAL: Target cut '{target_cut_name_to_stop_after}' for truncation was not found in full_cuts_dict. Aborting."); sys.exit(1)
    print(f"ℹ️ Using a reduced set of cuts. Applying {len(active_cuts_dict)} cuts up to and including: '{target_cut_name_to_stop_after}'.")

    collected_fit_data = prepare_fit_data_parallel(
        data_file, mc_files, mc_labels_list, tree, active_cuts_dict, variables_list,
        mixed_file, mixed_tree, empty_flask_file, GOOD_SPILLS_FILE_PATH
    )

    print("ℹ️ Cut table generation and display has been removed for speedup.")

    if collected_fit_data and isinstance(collected_fit_data, dict) and len(collected_fit_data) > 0:
        scipy_param_bounds_dict = {
            "c_mixed": (0.0, 10.0), "c_empty": (0.0, 10.0),
            "c_dy": (0.0, 2.0), "c_jpsi": (0.0, 2.0), "c_psip": (0.0, 2.0)
        }
        scipy_initial_guesses_dict = {
            "c_mixed": 1.0, "c_empty": 1.0, "c_dy": 0.1, "c_jpsi": 0.1, "c_psip": 0.1
        }
        print(f"ℹ️ SciPy FIT (Direct Data Model Modified): Using parameter bounds: {scipy_param_bounds_dict}")
        print(f"ℹ️ SciPy FIT (Direct Data Model Modified): Using initial guesses: {scipy_initial_guesses_dict}")

        n_bins_for_fit = 140
        print(f"ℹ️ SciPy FIT: Using n_bins = {n_bins_for_fit} for fit histograms.")
        hist_mass_min_GeV = 2.0 # Plotting range min
        hist_mass_max_GeV = 9.0 # Plotting range max
        fit_sub_range_min_GeV = 3.0 # Fit sub-range min
        fit_sub_range_max_GeV = 9.0 # Fit sub-range max
        print(f"ℹ️ Histogram mass range: {hist_mass_min_GeV:.2f}-{hist_mass_max_GeV:.2f} GeV. Fit sub-range: {fit_sub_range_min_GeV:.2f}-{fit_sub_range_max_GeV:.2f} GeV.")
        use_log_y_on_plot = True

        fit_results = perform_scipy_fit(
            collected_fit_data,
            mass_min=hist_mass_min_GeV, mass_max=hist_mass_max_GeV, n_bins=n_bins_for_fit,
            fit_range_gev_min=fit_sub_range_min_GeV, fit_range_gev_max=fit_sub_range_max_GeV,
            parameter_constraints=scipy_param_bounds_dict,
            initial_param_guesses=scipy_initial_guesses_dict,
            use_logy_plot=use_log_y_on_plot
        )

        if fit_results:
            print("\n--- SciPy Fit Results (Direct Data Model Modified) ---")
            fit_success = fit_results.get('success', False)
            fit_status = fit_results.get('status', -999)
            fit_message = fit_results.get('message', "N/A")
            print(f"   Fit Success (SciPy): {fit_success}")
            print(f"   Fit Status Code (SciPy): {fit_status}")
            print(f"   Fit Message (SciPy): {fit_message}")

            print("\n   --- Fitted Scaling Factors (c_values) ---")
            param_names_to_report = ["c_mixed", "c_empty", "c_dy", "c_jpsi", "c_psip"]
            fitted_params_values = {}
            for p_name in param_names_to_report:
                val = fit_results.get(p_name, np.nan)
                err_val = fit_results.get(f"err_{p_name}", np.nan)
                fitted_params_values[p_name] = (val, err_val)
                if not np.isnan(val): print(f"   {p_name:<15} : {format_param_err_str(val, err_val)}")

            chi2 = fit_results.get('chi2', np.nan)
            ndf = fit_results.get('ndf', -1)
            err_chi2_ndf = fit_results.get('err_chi2_ndf', np.nan)
            chi2_ndf_val = np.nan

            if ndf > 0:
                chi2_ndf_val = chi2 / ndf
                if not np.isnan(err_chi2_ndf): print(f"   Chi2/NDF: {chi2_ndf_val:.2f} +/- {err_chi2_ndf:.2f} ({chi2:.2f} / {ndf})")
                else: print(f"   Chi2/NDF: {chi2_ndf_val:.2f} ({chi2:.2f} / {ndf}) (error N/A)")
            else: print(f"   Chi2: {chi2:.2f} (NDF: {ndf}, Chi2/NDF not meaningful)")

            print(f"\nℹ️ To view plots, check 'scipy_fit_direct_data_model_modified.png/pdf'.")
            csv_file_name = "fit_results_direct_model_modified_output.csv"
            csv_data = []
            for p_name_key in param_names_to_report:
                if p_name_key in fit_results:
                    val, err_val = fitted_params_values.get(p_name_key, (np.nan, np.nan))
                    csv_data.append({
                        'Component': p_name_key, 'ScalingFactor': val, 'ErrorScalingFactor': err_val,
                        'Chisq_NDF': chi2_ndf_val if ndf > 0 else 'N/A',
                        'Error_Chisq_NDF': err_chi2_ndf if ndf > 0 and not np.isnan(err_chi2_ndf) else 'N/A',
                        'MinMassFit': fit_sub_range_min_GeV, 'MaxMassFit': fit_sub_range_max_GeV
                    })
            if csv_data:
                df_results = pd.DataFrame(csv_data)
                column_order = ['Component', 'ScalingFactor', 'ErrorScalingFactor', 'Chisq_NDF', 'Error_Chisq_NDF', 'MinMassFit', 'MaxMassFit']
                df_results = df_results[column_order]
                try:
                    df_results.to_csv(csv_file_name, index=False, float_format='%.6e')
                    print(f"\n✅ Fit results successfully saved to: {csv_file_name}")
                except Exception as e: print(f"\n❌ Error saving fit results to CSV: {e}")
            else: print("\nℹ️ No active parameters were fitted (or values were NaN), CSV file not generated.")

            if not fit_success : print("\n⚠️ WARNING: SciPy fit FAILED (check status and message).")
            elif ndf > 0 and (chi2 / ndf > 10.0) : print(f"\n⚠️ WARNING: SciPy fit converged, but Chi2/NDF ({chi2/ndf:.2f}) is very high.")
            elif ndf <= 0 and fit_success: print(f"\n⚠️ WARNING: SciPy fit converged, but NDF ({ndf}) is not positive (Chi2/NDF not meaningful).")
        else: print("❌ SciPy fit analysis did not return a results dictionary.")
    else: print("⚠️ No data collected for SciPy fit analysis, or data was empty/invalid.")
    print("\nScript finished.")