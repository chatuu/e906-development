import uproot
import pandas as pd
from tabulate import tabulate
import numpy as np

def create_cut_table(data_file_path, mc_file_paths, mc_labels, tree_name, cuts, variables, mixed_file_path, mixed_tree_name):
    cut_table = pd.DataFrame()

    def apply_cuts(df, cuts, use_weights=True):
        cut_results = {}
        weight = df['ReWeight'] if use_weights and 'ReWeight' in df.columns else pd.Series(1, index=df.index)

        cut_results["Total Events"] = weight.sum()

        for cut_name, cut_string in cuts.items():
            try:
                mask = df.eval(cut_string)
                cut_results[cut_name] = weight[mask].sum()
                df = df[mask]
                weight = weight[mask]
            except Exception as e:
                print(f"⚠️ Error applying cut '{cut_name}': {e}")
                cut_results[cut_name] = 0

        return pd.Series(cut_results)

    def read_tree(file_path, tree_name, variables):
        with uproot.open(file_path) as file:
            tree = file[tree_name]
            all_keys = tree.keys()
            vars_to_read = [var for var in variables if var in all_keys]
            missing = [var for var in variables if var not in all_keys]
            if missing:
                print(f"⚠️ Missing variables in {file_path}: {missing}")
            df = tree.arrays(vars_to_read, library="pd")
            return df

    # --- Data ---
    try:
        data = read_tree(data_file_path, tree_name, variables)
        cut_table["Data"] = apply_cuts(data.copy(), cuts, use_weights=False)
    except Exception as e:
        print(f"❌ Error reading Data file: {e}")
        return None

    # --- Data(RS67) ---
    try:
        rs67_data = read_tree(mixed_file_path, "result", variables)
        cut_table["Data(RS67)"] = apply_cuts(rs67_data.copy(), cuts, use_weights=False)
    except Exception as e:
        print(f"❌ Error reading RS67 data: {e}")
        return None

    # --- Mixed(RS67) ---
    try:
        mixed_data = read_tree(mixed_file_path, mixed_tree_name, variables)
        cut_table["Mixed(RS67)"] = apply_cuts(mixed_data.copy(), cuts, use_weights=False)
    except Exception as e:
        print(f"❌ Error reading Mixed file: {e}")
        return None

    # --- MC files ---
    for i, mc_file_path in enumerate(mc_file_paths):
        label = mc_labels[i]
        try:
            mc_data = read_tree(mc_file_path, tree_name, variables + ["ReWeight"])
            cut_table[label] = apply_cuts(mc_data.copy(), cuts, use_weights=True)
        except Exception as e:
            print(f"❌ Error reading MC file '{label}': {e}")
            cut_table[label] = None

    # --- Total MC ---
    try:
        cut_table["Total MC"] = cut_table[mc_labels].sum(axis=1)
    except Exception as e:
        print(f"⚠️ Error calculating Total MC: {e}")

    # --- Purity & Efficiency for DY MC ---
    try:
        dy_col = "DY MC"
        total_dy = cut_table.loc["Total Events", dy_col]
        cut_table["Purity (DY MC)"] = (cut_table[dy_col] / cut_table["Total MC"]) * 100
        cut_table["Efficiency (DY MC)"] = (cut_table[dy_col] / total_dy) * 100
    except Exception as e:
        print(f"⚠️ Error calculating Purity/Efficiency: {e}")

    cut_table.index.name = "Cut Name"
    return cut_table


# ------------------- Configuration -------------------

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
    "nHits1", "nHits2", "nHits1St1", "nHits2St1", "chisq_dimuon", "chisq1_target", "chisq2_target",
    "dx", "dy", "dz", "dpx", "dpy", "dpz", "mass", "D1", "D2", "D3", "xF", "xT", "xB", "costh", "intensityP"
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

    for col in ["Data", "Data(RS67)", "Mixed(RS67)"]:
        if col in display_df.columns:
            display_df[col] = pd.to_numeric(display_df[col], errors="coerce").fillna(0).astype(int)

    # Define a function to format based on the column
    def format_value(value, col_name):
        if col_name in ["Data", "Data(RS67)", "Mixed(RS67)"]:
            return "{:.0f}".format(value)
        else:
            return "{:.2f}".format(value)

    # Apply the formatting row by row
    formatted_data = []
    for index, row in display_df.iterrows():
        formatted_row = [format_value(value, col) for col, value in row.items()]
        formatted_data.append(formatted_row)

    print(tabulate(formatted_data, headers=display_df.columns, tablefmt="fancy_grid"))
    cut_table_df.to_csv("cut_table_with_mixed_and_reweighting.csv")
    print("\n✅ Cut table saved as 'cut_table_with_mixed_and_reweighting.csv'")