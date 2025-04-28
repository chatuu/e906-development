import uproot
import pandas as pd
from tabulate import tabulate
import numpy as np

def create_cut_table(data_file_path, mc_file_paths, mc_labels, tree_name, cuts, variables, mixed_file_path, mixed_tree_name):
    """
    Creates a cut flow table by applying sequential cuts to data, MC, and mixed event samples.

    Reads data from ROOT files, calculates a run-dependent beam offset, applies a
    series of cuts, and reports the number of events remaining after each cut.
    Includes calculations for total MC, purity, and efficiency for the DY MC sample.

    Parameters
    ----------
    data_file_path : str
        Path to the data ROOT file.
    mc_file_paths : list of str
        List of paths to the MC ROOT files.
    mc_labels : list of str
        List of labels corresponding to the MC files (used for column names).
    tree_name : str
        Name of the tree within the ROOT files containing the event data
        (used for data and MC files).
    cuts : dict
        A dictionary where keys are descriptive cut names and values are
        string expressions parsable by pandas.DataFrame.eval(). Cuts are applied
        sequentially in the order they appear in the dictionary.
    variables : list of str
        List of variable names to read from the ROOT trees. This list should
        include all variables used in the cut expressions, plus 'runID' and 'dy'.
    mixed_file_path : str
        Path to the file containing both Data(RS67) and Mixed(RS67) samples.
    mixed_tree_name : str
        Name of the tree within the mixed file containing the Mixed(RS67) data.
        Assumes the Data(RS67) is in a tree named "result".

    Returns
    -------
    pd.DataFrame or None
        A pandas DataFrame representing the cut flow table with cut names as index
        and samples (Data, MCs, Total MC, Mixed, Purity, Efficiency) as columns.
        Returns None if a critical error occurs during data processing.
    """
    cut_table = pd.DataFrame()

    def apply_cuts(df, cuts, use_weights=True):
        """
        Applies a dictionary of sequential cuts to a pandas DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame containing event data.
        cuts : dict
            A dictionary of cut names and their corresponding string expressions.
            Cuts are applied sequentially.
        use_weights : bool, optional
            Whether to use the 'ReWeight' column for summing events. Defaults to True.

        Returns
        -------
        pd.Series
            A pandas Series where the index is the cut name (including "Total Events")
            and the value is the sum of event weights (or counts if no weights)
            remaining after applying the cut and all previous cuts.
        """
        cut_results = {}
        # Use 'ReWeight' column if it exists and use_weights is True, otherwise use count of 1
        weight = df['ReWeight'] if use_weights and 'ReWeight' in df.columns else pd.Series(1.0, index=df.index)

        # Ensure weight is numeric
        weight = pd.to_numeric(weight, errors='coerce').fillna(1.0)

        cut_results["Total Events"] = weight.sum()

        # Ensure beamOffset column exists before cuts that might use it
        if 'beamOffset' not in df.columns:
             df['beamOffset'] = 0.0 # Default to 0 if calculation failed or column is missing

        # Ensure dy is numeric before cuts that use it
        if 'dy' in df.columns:
             df['dy'] = pd.to_numeric(df['dy'], errors='coerce').fillna(np.nan)

        # Keep track of the mask from successful cuts
        cumulative_mask = pd.Series(True, index=df.index)


        for cut_name, cut_string in cuts.items():
            # Skip the total events entry if it somehow ended up in cuts
            if cut_name == "Total Events":
                 continue
            try:
                # Evaluate the cut string on the *current* filtered DataFrame
                mask = df.eval(cut_string, engine='python')

                # Apply the mask to the cumulative mask and the original weight series
                cumulative_mask = cumulative_mask & mask

                # Calculate the number of events remaining after this cut (using original weights)
                cut_results[cut_name] = weight[cumulative_mask].sum()

                # Filter the DataFrame and weight series for the *next* iteration
                # We filter df here so subsequent cuts are applied to the reduced dataset
                df = df[mask].copy() # Use .copy() to avoid SettingWithCopyWarning
                # Note: weight is *not* filtered here; we apply the cumulative_mask to the *original* weight series later
                # This is important for correct sequential summing. The filtered 'df' is used for the *next* df.eval()

            except Exception as e:
                # This catch block now handles actual errors from df.eval(),
                # e.g., if a column name used in the string is missing
                print(f"⚠️ Error applying cut '{cut_name}': {e}. Skipping cut for this dataset.")
                cut_results[cut_name] = 0
                # If a cut fails, the cumulative_mask is NOT updated by this cut's mask.
                # The dataframe 'df' is also NOT filtered by the failed mask.
                # The next cut will be applied to the dataframe filtered only by *previous successful* cuts.


        return pd.Series(cut_results)


    def read_tree(file_path, tree_name, variables):
        """
        Reads specified variables from a ROOT tree into a pandas DataFrame.

        Handles missing files, trees, and variables gracefully by returning an empty
        DataFrame on major errors or adding missing columns with NaN. Ensures
        essential columns ('runID', 'dy', 'ReWeight') are always attempted and
        converted to numeric types.

        Parameters
        ----------
        file_path : str
            Path to the ROOT file.
        tree_name : str
            Name of the tree within the file.
        variables : list of str
            List of variable names to read.

        Returns
        -------
        pd.DataFrame
            A pandas DataFrame containing the requested variables, plus 'runID', 'dy',
            and 'ReWeight' if available. Missing columns are added as NaN. Returns
            an empty DataFrame if the file or tree cannot be read.
        """
        try:
            with uproot.open(file_path) as file:
                if tree_name not in file:
                    print(f"❌ Error: Tree '{tree_name}' not found in file '{file_path}'.")
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
                vars_to_request = list(set(variables + ['runID', 'dy', 'ReWeight']))

                vars_to_read = [var for var in vars_to_request if var in all_keys]
                missing = [var for var in vars_to_request if var not in all_keys]
                if missing:
                    # Filter out 'ReWeight' from missing warning if it's likely data
                    is_mc_file = any(mc_f in file_path for mc_f in mc_files) # Basic check
                    if 'ReWeight' in missing and not is_mc_file:
                         missing.remove('ReWeight')
                         if missing:
                             print(f"⚠️ Missing variables in {file_path}: {missing}. These will be added as NaN.")
                         # else: print(f"✅ All expected variables except ReWeight found in {file_path} (assuming data).") # Optional verbose
                    elif missing:
                         print(f"⚠️ Missing variables in {file_path}: {missing}. These will be added as NaN.")


                df = tree.arrays(vars_to_read, library="pd")

                # Add missing columns with NaN values to ensure DataFrame has expected structure
                for var in missing:
                     if var not in df.columns: # Avoid adding if it somehow got added during read
                         df[var] = np.nan

                # Ensure runID, dy, and ReWeight are numeric (important for beamOffset, cuts, and weights)
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
            print(f"❌ Error reading tree '{tree_name}' from '{file_path}': {e}")
            # Return empty DataFrame on major error for structure consistency
            empty_df = pd.DataFrame()
            for var in set(variables + ['runID', 'dy', 'ReWeight']):
                 empty_df[var] = pd.Series(dtype='float64')
            return empty_df


    def add_beam_offset(df):
        """
        Calculates and adds the 'beamOffset' column to the DataFrame based on 'runID'.

        Uses a fixed set of run ID ranges to determine the beam offset value for
        each event. Events outside the specified ranges will have a beamOffset of 0.0.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame which must contain a 'runID' column.

        Returns
        -------
        pd.DataFrame
            The DataFrame with an added 'beamOffset' column. Returns the original
            DataFrame with beamOffset set to 0.0 if 'runID' column is missing or all NaN.
        """
        if 'runID' not in df.columns or df['runID'].isnull().all():
            print("⚠️ Cannot calculate beamOffset: 'runID' column not found or is all NaN. Setting beamOffset to 0.0.")
            df['beamOffset'] = 0.0
            return df

        # Fill NaN in runID with a value outside all ranges (e.g., -1) before calculation
        # This ensures np.select doesn't raise errors on NaN input and maps NaNs to the default value
        runids_for_calc = df['runID'].fillna(-1).astype(int)

        conditions = [
            (runids_for_calc >= 8912) & (runids_for_calc <= 10420),
            (runids_for_calc >= 10421) & (runids_for_calc <= 10912),
            (runids_for_calc >= 11075) & (runids_for_calc <= 12435),
            (runids_for_calc >= 12525) & (runids_for_calc <= 15789),
            (runids_for_calc >= 15793) & (runids_for_calc <= 16076)
        ]
        values = [0.4, 0.4, 1.6, 1.6, 1.6]

        # Calculate beamOffset, default to 0.0 if no condition is met (including original NaN runIDs)
        df['beamOffset'] = np.select(conditions, values, default=0.0)

        return df


    # --- Data ---
    try:
        data = read_tree(data_file_path, tree_name, variables)
        data = add_beam_offset(data) # Calculate beamOffset
        # Ensure all variables needed for cuts plus beamOffset and dy exist before applying cuts
        # This helps apply_cuts avoid errors for datasets missing specific variables used in cuts
        required_for_cuts = list(set(variables + ['beamOffset', 'dy']))
        for var in required_for_cuts:
             if var not in data.columns:
                  data[var] = np.nan # Add any missing required variables
        cut_table["Data"] = apply_cuts(data.copy(), cuts, use_weights=False)
    except Exception as e:
        print(f"❌ Error processing Data file: {e}")
        # Ensure Data column exists even on error for table structure consistency
        # Need to know the expected cut names to create a Series with zeros
        expected_cut_names = ["Total Events"] + list(cuts.keys())
        cut_table["Data"] = pd.Series(0.0, index=expected_cut_names)


    # --- Data(RS67) ---
    # Assuming 'result' tree also has runID and dy
    try:
        rs67_data = read_tree(mixed_file_path, "result", variables)
        rs67_data = add_beam_offset(rs67_data) # Calculate beamOffset
        # Ensure all cut variables plus beamOffset and dy exist
        required_for_cuts = list(set(variables + ['beamOffset', 'dy']))
        for var in required_for_cuts:
             if var not in rs67_data.columns:
                  rs67_data[var] = np.nan # Add any missing required variables
        cut_table["Data(RS67)"] = apply_cuts(rs67_data.copy(), cuts, use_weights=False)
    except Exception as e:
        print(f"❌ Error processing RS67 data: {e}")
        expected_cut_names = ["Total Events"] + list(cuts.keys())
        cut_table["Data(RS67)"] = pd.Series(0.0, index=expected_cut_names)


    # --- Mixed(RS67) ---
    # Assuming mixed_tree also has runID and dy
    try:
        mixed_data = read_tree(mixed_file_path, mixed_tree_name, variables)
        mixed_data = add_beam_offset(mixed_data) # Calculate beamOffset
        # Ensure all cut variables plus beamOffset and dy exist
        required_for_cuts = list(set(variables + ['beamOffset', 'dy']))
        for var in required_for_cuts:
             if var not in mixed_data.columns:
                  mixed_data[var] = np.nan # Add any missing required variables
        cut_table["Mixed(RS67)"] = apply_cuts(mixed_data.copy(), cuts, use_weights=False)
    except Exception as e:
        print(f"❌ Error processing Mixed file: {e}")
        expected_cut_names = ["Total Events"] + list(cuts.keys())
        cut_table["Mixed(RS67)"] = pd.Series(0.0, index=expected_cut_names)


    # --- MC files ---
    mc_results = {} # Store individual MC results Series by label
    for i, mc_file_path in enumerate(mc_file_paths):
        label = mc_labels[i]
        try:
            # Read variables + ReWeight + runID + dy for MC
            mc_data = read_tree(mc_file_path, tree_name, variables + ["ReWeight"])
            mc_data = add_beam_offset(mc_data) # Calculate beamOffset
             # Ensure all cut variables plus beamOffset and dy exist
            required_for_cuts = list(set(variables + ['beamOffset', 'dy']))
            for var in required_for_cuts:
                if var not in mc_data.columns:
                     mc_data[var] = np.nan # Add any missing required variables

            mc_results[label] = apply_cuts(mc_data.copy(), cuts, use_weights=True)
            cut_table[label] = mc_results[label] # Add to cut_table immediately

        except Exception as e:
            print(f"❌ Error processing MC file '{label}': {e}")
            # Ensure MC column exists even on error
            expected_cut_names = ["Total Events"] + list(cuts.keys())
            cut_table[label] = pd.Series(0.0, index=expected_cut_names)
            mc_results[label] = pd.Series(0.0, index=expected_cut_names) # Add zero series for summation


    # --- Total MC ---
    try:
         # Sum the results from the mc_results dictionary
         if mc_results:
             # Ensure all Series in mc_results have the same index before summing
             # This handles cases where one dataset failed and returned a Series with zeros
             all_mc_series = [s for s in mc_results.values()]
             if all_mc_series:
                 # Align indices - fill missing cuts with 0 for summation
                 cut_table["Total MC"] = pd.concat(all_mc_series, axis=1).sum(axis=1, min_count=1).fillna(0)
             else:
                  cut_table["Total MC"] = pd.Series(0.0, index=["Total Events"] + list(cuts.keys()))

         else:
             cut_table["Total MC"] = pd.Series(0.0, index=["Total Events"] + list(cuts.keys()))

    except Exception as e:
        print(f"⚠️ Error calculating Total MC: {e}")
        cut_table["Total MC"] = pd.Series(0.0, index=["Total Events"] + list(cuts.keys()))


    # --- Purity & Efficiency for DY MC ---
    try:
        dy_col = "DY MC"
        # Ensure columns and index exist before calculation
        if dy_col in cut_table.columns and "Total Events" in cut_table.index and "Total MC" in cut_table.columns:
             total_dy = cut_table.loc["Total Events", dy_col]
             # Avoid division by zero, handle inf/nan results
             # Use .copy() here to prevent SettingWithCopyWarning if modifying in place later
             # Ensure divisor is not zero for purity
             total_mc_safe = cut_table["Total MC"].replace(0, np.nan) # Replace 0 with NaN for division
             purity = (cut_table[dy_col] / total_mc_safe).replace([np.inf, -np.inf], np.nan).fillna(0).copy()

             # Ensure total_dy is not zero for efficiency
             total_dy_safe = total_dy if total_dy != 0 else np.nan
             efficiency = (cut_table[dy_col] / total_dy_safe).replace([np.inf, -np.inf], np.nan).fillna(0).copy()


             cut_table["Purity (DY MC)"] = purity * 100
             cut_table["Efficiency (DY MC)"] = efficiency * 100
        else:
            print(f"⚠️ Cannot calculate Purity/Efficiency: Missing '{dy_col}' column, 'Total Events' index, or 'Total MC' column. Or Total Events for DY MC is zero.")
            cut_table["Purity (DY MC)"] = 0.0
            cut_table["Efficiency (DY MC)"] = 0.0

    except Exception as e:
        print(f"⚠️ Error calculating Purity/Efficiency: {e}")
        cut_table["Purity (DY MC)"] = 0.0
        cut_table["Efficiency (DY MC)"] = 0.0


    cut_table.index.name = "Cut Name"

    # Ensure "Total Events" row is the first index and all cuts from cuts_dict are present
    expected_indices = ["Total Events"] + list(cuts.keys())
    # Reindex the DataFrame to ensure the correct order and presence of all cuts (filling missing rows with 0.0)
    cut_table = cut_table.reindex(expected_indices, fill_value=0.0)


    return cut_table


# ------------------- Configuration -------------------

cuts_dict = {
    "nhits1 + nhits2 > 29": "(nHits1 + nHits2) > 29",
    "nhits1st1 + nhits2st1 > 8": "(nHits1St1 + nHits2St1) > 8",
    "chisq dimuon < 18": "chisq_dimuon < 18",
    "chisq Target within 2": "abs(chisq1_target + chisq2_target - chisq_dimuon) < 2",
    "dx within -0.25 to 0.25": "dx > -0.25 and dx < 0.25",
    "dz within -280 to -5": "dz > -280 and dz < -5",
    # Re-adding the dy - beamOffset cut
    #"dy - beamOffset within -0.22 to 0.22": "(dy - beamOffset) > -0.22 and (dy - beamOffset) < 0.22",
    "dy - beamOffset within -0.22 to 0.22": "((abs(dy-1.6) < 0.22 and runID > 11000) or (abs(dy-0.4) < 0.22 and runID < 11000))", # Example alternative
    "dx^2 + (dy - beamOffset)^2 < 0.06":"((dx*dx+(dy-1.6)*(dy-1.6)<0.06 and runID > 11000) or (dx*dx+(dy-0.4)*(dy-0.4)<0.06 and runID < 11000))", # Example alternative
    "abs(x1_st1 + x2_st1) < 42cm": "abs(x1_st1 + x2_st1) < 42",
    "abs(trackSeparation) < 270": "abs(trackSeparation) < 270",
    "abs(dpx) < 1.8": "abs(dpx) < 1.8",
    "abs(dpy) < 2": "abs(dpy) < 2",
    "dpz within 38 to 116": "dpz > 38 and dpz < 116",
    "dpx^2 + dpy^2 < 5": "(dpx**2 + dpy**2) < 5",
    "mass within 4.2 to 8.8": "mass > 4.2 and mass < 8.8",
    "xF within -0.1 to 0.95": "xF > -0.1 and xF < 0.95",
    "xT within -0.1 to 0.58": "xT > -0.1 and xT < 0.58",
    "cosTheta within -0.5 to 0.5": "costh > -0.5 and costh < 0.5",
    "D1 < 400": "D1 < 400",
    "D2 < 400": "D2 < 400",
    "D3 < 400": "D3 < 400",
    "D1 + D2 + D3 < 1000": "D1 + D2 + D3 < 1000",
    "intensity within 0 to 80000": "intensityP > 0 and intensityP < 80000"
}

variables_list = [
    "nHits1", "nHits2", "nHits1St1", "nHits2St1", "chisq_dimuon", "chisq1_target", "chisq2_target",
    "dx", "dy", "dz", "dpx", "dpy", "dpz", "mass", "D1", "D2", "D3", "xF", "xT", "xB", "costh", "intensityP",
    "runID", "trackSeparation", "x1_st1", "x2_st1" # <--- Need runID for beamOffset calculation, trackSeparation added from cuts_dict
    # 'beamOffset' is not in ROOT files, it's calculated, so don't add it here.
    # 'ReWeight' is handled specifically in read_tree
]

data_file = "../../ROOTFiles/Hugo/roadset57_70_R008_2111v42_tmp_noPhys.root"

mc_files = [
    "../../ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    "../../ROOTFiles/Hugo/mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    "../../ROOTFiles/Hugo/mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
]

mc_labels_list = ["DY MC", "J/Psi MC", "Psi Prime MC"]

tree = "Tree"
mixed_file = "../../ROOTFiles/MixedEvents/merged_RS67_3089LH2.root"
mixed_tree = "result_mix"

# ------------------- Generate Cut Table -------------------

cut_table_df = create_cut_table(data_file, mc_files, mc_labels_list, tree, cuts_dict, variables_list, mixed_file, mixed_tree)

if cut_table_df is not None:
    display_df = cut_table_df.copy()

    # Ensure "Total Events" is the first cut name and all cuts from cuts_dict are included
    cut_names_order = ["Total Events"] + list(cuts_dict.keys())
    # Reindex display_df to match the desired order and include all cut names
    display_df = display_df.reindex(cut_names_order)


    display_df.insert(0, "Cut", display_df.index) # Use the index for the 'Cut' column

    # Define columns that should be formatted as integers (event counts)
    count_cols = ["Data", "Data(RS67)", "Mixed(RS67)", "Total MC"] + mc_labels_list
    for col in count_cols:
        if col in display_df.columns:
            # Use .round(0) and .astype(int) to ensure they are integers
            display_df[col] = pd.to_numeric(display_df[col], errors="coerce").fillna(0).round(0).astype(int)


    # Define a function to format based on the column
    def format_value(value, col_name):
        if col_name in count_cols:
            return "{:.0f}".format(value)
        elif col_name in ["Purity (DY MC)", "Efficiency (DY MC)"]: # Percentage columns
             # Handle potential non-numeric after fillna(0) if source was bad
             try:
                 return "{:.2f}".format(float(value))
             except (ValueError, TypeError):
                  return str(value) # Fallback if conversion fails
        else:
             # Default formatting for unexpected columns
             return str(value)

    # Apply the formatting row by row
    # Use display_df.columns to get columns excluding the 'Cut' column we just inserted
    data_to_display = []
    for index, row in display_df.iterrows():
        formatted_row = [row["Cut"]] + [format_value(row[col], col) for col in display_df.columns[1:]]
        data_to_display.append(formatted_row)


    print(tabulate(data_to_display, headers=["Cut"] + list(display_df.columns[1:]), tablefmt="fancy_grid"))
    cut_table_df.to_csv("cut_table_with_mixed_and_reweighting.csv")
    print("\n✅ Cut table saved as 'cut_table_with_mixed_and_reweighting.csv'")