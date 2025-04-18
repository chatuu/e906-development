import uproot
import pandas as pd
from tabulate import tabulate
import numpy as np

def create_cut_table(data_file_path, mc_file_paths, mc_labels, tree_name, cuts, variables, mixed_file_path, mixed_tree_name):
    """
    Reads data, Monte Carlo, and mixed event ROOT files, applies cuts, and creates a combined cut table
    with 'Total MC', 'Purity', 'Efficiency', and 'Mixed' columns.

    Args:
        data_file_path (str): Path to the data ROOT file.
        mc_file_paths (list): List of paths to the Monte Carlo ROOT files.
        mc_labels (list): List of strings for the column headings of the MC files.
        tree_name (str): Name of the TTree in the data and MC files.
        cuts (dict): Dictionary of cuts, where keys are cut names and values are cut strings.
        variables (list): List of variables to be read from the TTrees.
        mixed_file_path (str): Path to the mixed events ROOT file.
        mixed_tree_name (str): Name of the TTree in the mixed events file.

    Returns:
        pandas.DataFrame: DataFrame containing the combined cut table with 'Total MC', 'Purity', 'Efficiency', and 'Mixed' columns.
    """

    all_tables = {}

    # Process Data file
    try:
        file = uproot.open(data_file_path)
        tree = file[tree_name]
        if tree is None:
            raise ValueError(f"Tree '{tree_name}' not found in data file: {data_file_path}")
        available_data_vars = set(tree.keys())
        data_vars_to_read = [var for var in variables if var in available_data_vars]
        data = tree.arrays(data_vars_to_read, library="pd")
        file.close()

        cut_table_data = pd.DataFrame(columns=["Data"])
        cut_table_data.loc["Total Events"] = len(data)
        data_cut_applied = data.copy() # Keep track of data after cuts

        for cut_name, cut_string in cuts.items():
            try:
                cut_mask = data_cut_applied.eval(cut_string)
                data_cut_applied = data_cut_applied[cut_mask]
                cut_table_data.loc[cut_name] = len(data_cut_applied)
            except Exception as e:
                print(f"Error applying cut '{cut_name}' to data: {e}")
                cut_table_data.loc[cut_name] = None
        all_tables["Data"] = cut_table_data

    except FileNotFoundError:
        print(f"Error: Data ROOT file '{data_file_path}' not found.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred with data file: {e}")
        return None

    # Process Mixed Events file
    try:
        file = uproot.open(mixed_file_path)
        mixed_tree = file[mixed_tree_name]
        if mixed_tree is None:
            raise ValueError(f"Tree '{mixed_tree_name}' not found in mixed events file: {mixed_file_path}")
        available_mixed_vars = set(mixed_tree.keys())
        mixed_vars_to_read = [var for var in variables if var in available_mixed_vars]
        mixed_data = mixed_tree.arrays(mixed_vars_to_read, library="pd")
        file.close()

        cut_table_mixed = pd.DataFrame(columns=["Mixed"])
        cut_table_mixed.loc["Total Events"] = len(mixed_data)
        mixed_data_cut_applied = mixed_data.copy() # Keep track of mixed data after cuts

        for cut_name, cut_string in cuts.items():
            try:
                required_vars_in_cut = set()
                for var in variables:
                    if var in cut_string:
                        required_vars_in_cut.add(var)

                if required_vars_in_cut.issubset(set(mixed_data_cut_applied.columns)):
                    cut_mask = mixed_data_cut_applied.eval(cut_string)
                    mixed_data_cut_applied = mixed_data_cut_applied[cut_mask]
                    cut_table_mixed.loc[cut_name] = len(mixed_data_cut_applied)
                else:
                    print(f"Warning: Not all required variables for cut '{cut_name}' found in mixed events file. Skipping this cut for mixed events.")
                    cut_table_mixed.loc[cut_name] = None
            except Exception as e:
                print(f"Error applying cut '{cut_name}' to mixed events: {e}")
                cut_table_mixed.loc[cut_name] = None
        all_tables["Mixed"] = cut_table_mixed

    except FileNotFoundError:
        print(f"Error: Mixed events ROOT file '{mixed_file_path}' not found.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred with mixed events file: {e}")
        return None

    # Process Monte Carlo files
    mc_tables = []
    for i, mc_file_path in enumerate(mc_file_paths):
        mc_label = mc_labels[i]
        try:
            file = uproot.open(mc_file_path)
            tree = file[tree_name]
            if tree is None:
                raise ValueError(f"Tree '{tree_name}' not found in MC file: {mc_file_path}")
            available_mc_vars = set(tree.keys())
            mc_vars_to_read = [var for var in variables if var in available_mc_vars]
            mc_data = tree.arrays(mc_vars_to_read, library="pd")
            file.close()

            mc_cut_table = pd.DataFrame(columns=[mc_label])
            mc_cut_table.loc["Total Events"] = len(mc_data)
            mc_data_cut_applied = mc_data.copy() # Keep track of MC data after cuts

            for cut_name, cut_string in cuts.items():
                try:
                    required_vars_in_cut = set()
                    for var in variables:
                        if var in cut_string:
                            required_vars_in_cut.add(var)

                    if required_vars_in_cut.issubset(set(mc_data_cut_applied.columns)):
                        cut_mask = mc_data_cut_applied.eval(cut_string)
                        mc_data_cut_applied = mc_data_cut_applied[cut_mask]
                        mc_cut_table.loc[cut_name] = len(mc_data_cut_applied)
                    else:
                        print(f"Warning: Not all required variables for cut '{cut_name}' found in MC file for {mc_label}. Skipping this cut for this MC file.")
                        mc_cut_table.loc[cut_name] = None

                except Exception as e:
                    print(f"Error applying cut '{cut_name}' to MC file for {mc_label}: {e}")
                    mc_cut_table.loc[cut_name] = None

            mc_tables.append(mc_cut_table)
            all_tables[mc_label] = mc_cut_table

        except FileNotFoundError:
            print(f"Error: MC ROOT file '{mc_file_path}' not found.")
            pass
        except Exception as e:
            print(f"An unexpected error occurred with MC file for {mc_label}: {e}")
            pass

    # Combine all tables
    combined_cut_table = pd.concat(all_tables, axis=1)

    # Calculate Total MC and insert after 'Data' column
    data_col_loc = combined_cut_table.columns.get_loc('Data')
    data_col_index = int(data_col_loc) if not isinstance(data_col_loc, slice) else int(data_col_loc.start)

    # Get the 'Mixed' column Series
    mixed_series = combined_cut_table[('Mixed', 'Mixed')]
    # Drop the original 'Mixed' column
    combined_cut_table = combined_cut_table.drop(('Mixed', 'Mixed'), axis=1)
    # Insert the 'Mixed' column after 'Data'
    combined_cut_table.insert(data_col_index + 1, 'Mixed', mixed_series)
    combined_cut_table.insert(data_col_index + 2, 'Total MC', combined_cut_table[mc_labels].sum(axis=1))

    # Calculate Purity (DY MC / Total MC * 100) and insert after 'Psi Prime MC'
    if 'DY MC' in combined_cut_table.columns and 'Psi Prime MC' in combined_cut_table.columns and 'Total MC' in combined_cut_table.columns:
        try:
            dy_mc_series = combined_cut_table['DY MC'].iloc[:, 0].astype(float)
            total_mc_series = combined_cut_table['Total MC'].astype(float)

            purity = (dy_mc_series / total_mc_series) * 100
            psi_prime_mc_loc = combined_cut_table.columns.get_loc(('Psi Prime MC', 'Psi Prime MC'))
            if isinstance(psi_prime_mc_loc, (int, slice)):
                psi_prime_mc_index = int(psi_prime_mc_loc) if not isinstance(psi_prime_mc_loc, slice) else int(psi_prime_mc_loc.start)
            elif isinstance(psi_prime_mc_loc, (list, tuple, pd.Index, np.ndarray)):
                psi_prime_mc_index = int(np.asarray(psi_prime_mc_loc).flatten()[0])
            else:
                raise TypeError(f"Unexpected type for column location: {type(psi_prime_mc_loc)}")

            # Insert Purity
            if 0 <= psi_prime_mc_index + 1 <= len(combined_cut_table.columns):
                combined_cut_table.insert(psi_prime_mc_index + 1, 'Purity (DY MC)', purity)
                purity_index = psi_prime_mc_index + 1
            else:
                combined_cut_table['Purity (DY MC)'] = purity # Append if index is out of bounds
                purity_index = len(combined_cut_table.columns) - 1 # Index of the last column
                print("Warning: Could not insert 'Purity (DY MC)' after 'Psi Prime MC', appending to the end.")

            # Calculate Efficiency and insert after 'Purity (DY MC)'
            total_dy_mc_events = combined_cut_table.loc['Total Events', ('DY MC', 'DY MC')]
            efficiency = (dy_mc_series / total_dy_mc_events) * 100

            # Ensure the insertion index for Efficiency is valid
            if 0 <= purity_index + 1 <= len(combined_cut_table.columns):
                combined_cut_table.insert(purity_index + 1, 'Efficiency (DY MC)', efficiency)
            else:
                combined_cut_table['Efficiency (DY MC)'] = efficiency # Append if index is out of bounds
                print("Warning: Could not insert 'Efficiency (DY MC)' after 'Purity (DY MC)', appending to the end.")

        except Exception as e:
            print(f"Error during purity/efficiency calculation: {e}")
            raise

    else:
        print("Warning: Could not calculate 'Purity (DY MC)' and 'Efficiency (DY MC)' columns. Ensure 'DY MC', 'Psi Prime MC', and 'Total MC' columns exist.")

    combined_cut_table.index.name = "Cut Name"
    return combined_cut_table

# Example usage
cuts_dict = {
    "nhits1 + nhits2 > 29": "(nHits1 + nHits2) > 29",
    "nhits1st1 + nhits2st1 > 8": "(nHits1St1 + nHits2St1) > 8",
    "chisq dimuon < 18": "chisq_dimuon < 18",
    "chisq Target within 2": "abs(chisq1_target + chisq2_target - chisq_dimuon) < 2",
    "dx within -0.25 to 0.25": "dx > -0.25 and dx < 0.25",
    "dz within -280 to -5": "dz > -280 and dz < -5",
    "abs(dpx) < 1.8": "abs(dpx) < 1.8",
    "abs(dpy) < 2": "abs(dpy) < 2",
    "dpz within 38 to 116": "dpz > 38 and dpz < 116",
    "dpx^2 + dpy^2 < 5": "(dpx**2 + dpy**2) < 5",
    "mass within 4.2 to 8.8": "mass > 4.2 and mass < 8.8",
    "xF within -0.1 to 0.95": "xF > -0.1 and xF < 0.95",
    "xT within 0.1 to 0.58": "xT > -0.1 and xT < 0.58",
    "cosTheta within -0.5 to 0.5": "costh > -0.5 and costh < 0.5",
    "D1 < 400": "D1 < 400",
    "D2 < 400": "D2 < 400",
    "D3 < 400": "D3 < 400",
    "D1 + D2 + D3 < 1000": "D1 + D2 + D3 < 1000",
    "intensity within 0 to 80000": "intensityP > 0 and intensityP < 80000"
}
variables_list = [
    "nHits1",
    "nHits2",
    "nHits1St1",
    "nHits2St1",
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
]

data_file = "/Users/ckuruppu/Documents/NMSU-Physics/e906-development/ROOTFiles/Hugo/roadset57_70_R008_2111v42_tmp_noPhys.root"
mc_files = [
    "/Users/ckuruppu/Documents/NMSU-Physics/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    "/Users/ckuruppu/Documents/NMSU-Physics/e906-development/ROOTFiles/Hugo/mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    "/Users/ckuruppu/Documents/NMSU-Physics/e906-development/ROOTFiles/Hugo/mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
]
mc_labels_list = ["DY MC", "J/Psi MC", "Psi Prime MC"]
tree = "Tree"
mixed_file = "/Users/ckuruppu/Documents/NMSU-Physics/e906-development/ROOTFiles/Hugo/merged_RS67_3089LH2.root"
mixed_tree = "result_mix"

cut_table_df = create_cut_table(data_file, mc_files, mc_labels_list, tree, cuts_dict, variables_list, mixed_file, mixed_tree)

if cut_table_df is not None:
    print(cut_table_df)

output_csv_file = "cut_table_with_mixed.csv"

if cut_table_df is not None:
    print(tabulate(cut_table_df, headers="keys", tablefmt="fancy_grid"))
    cut_table_df.to_csv(output_csv_file)
    print(f"\nCut table with mixed events saved to: {output_csv_file}")