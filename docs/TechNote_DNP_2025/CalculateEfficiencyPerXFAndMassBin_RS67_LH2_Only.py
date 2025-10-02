"""
Calculates Average Efficiency from NPZ Files.

This script reads a dataset of dimuon events from a Parquet file and a collection
of efficiency curves from .npz files. It then applies a series of physics and
kinematic cuts (dataCut) to the events, using a dynamic beam offset that depends
on the runID of each event. For the surviving events, it calculates the average
efficiency within a 2D grid of kinematic bins (xF vs. mass).

The final output is a single CSV file containing the average efficiency and
associated errors for each bin.
"""

from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import os

def get_plus_minus(val, x_array):
    """Finds the nearest lower and upper values in an array for a given value.

    Args:
        val (float): The value to search for.
        x_array (np.ndarray): The sorted array to search within.

    Returns:
        tuple: A tuple containing two tuples:
               ((lower_val, lower_idx), (upper_val, upper_idx)).
               Returns np.nan for values or indices that are not found.
    """
    x_array = np.asarray(x_array)

    lower_val, upper_val = None, None
    lower_idx, upper_idx = None, None

    # Find the closest value less than or equal to val
    le_val_mask = x_array <= val
    if np.any(le_val_mask):
        lower_idx = np.where(le_val_mask)[0][-1]
        lower_val = x_array[lower_idx]

    # Find the closest value greater than or equal to val
    ge_val_mask = x_array >= val
    if np.any(ge_val_mask):
        upper_idx = np.where(ge_val_mask)[0][0]
        upper_val = x_array[upper_idx]

    if lower_val is None:
        lower_val, lower_idx = np.nan, np.nan
    if upper_val is None:
        upper_val, upper_idx = np.nan, np.nan

    return (lower_val, lower_idx), (upper_val, upper_idx)


def effi_model(d2, e_d2, npz_filepath, model="interpolate"):
    """Calculates efficiency and its uncertainty using an interpolation model.

    This function loads an efficiency curve from a .npz file and uses linear
    interpolation to find the efficiency for a given series of input values (d2).

    Args:
        d2 (pd.Series): A series of values (e.g., D2) for which to find the efficiency.
        e_d2 (pd.Series): The errors associated with the d2 values (currently unused).
        npz_filepath (str): The full path to the .npz file containing the
            efficiency curve data ('x', 'y', 'y_error_low', 'y_error_high').
        model (str, optional): The type of model to use. Only 'interpolate' is
            supported. Defaults to "interpolate".

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
            - The first DataFrame holds the estimated efficiency values.
            - The second DataFrame holds the estimated errors on the efficiency.
    """
    data = np.load(npz_filepath)
    x, y = data["x"], data["y"]
    y_error_low, y_error_high = data["y_error_low"], data["y_error_high"]

    f_linear = interp1d(x, y, kind="linear", fill_value="extrapolate")
    d2 = pd.Series(d2)

    if model != "interpolate":
        raise ValueError("This script is configured for the 'interpolate' model.")

    effi_values = f_linear(d2)
    effi = pd.DataFrame(effi_values, columns=["estimate"])

    # A simplified error interpolation for calculating the average propagated error
    e_effi_vals = np.interp(d2, x, (y_error_low + y_error_high) / 2)
    e_effi = pd.DataFrame({"e_estimate": e_effi_vals})

    return effi, e_effi


def apply_datacut(df):
    """
    Applies the specified event selection criteria (dataCut) to the DataFrame.
    The beam offset is calculated dynamically for each row based on runID.

    Args:
        df (pd.DataFrame): The input DataFrame with event data.

    Returns:
        pd.DataFrame: A new DataFrame containing only the events that pass the cuts.
    """
    print("Applying dataCut with dynamic beam_offset based on runID...")
    initial_events = len(df)

    # Create a dynamic beam_offset Series based on the runID for each event
    # if runID > 11000, offset is 1.6, otherwise it's 0.4
    beam_offset = np.where(df['runID'] > 11000, 1.6, 0.4)

    # chuckCutsPositive_2111v42_tmp
    cut_pos = (
        (df['chisq1_target'] < 15) &
        (df['pz1_st1'] > 9) & (df['pz1_st1'] < 75) &
        (df['nHits1'] > 13) &
        (df['x1_t']**2 + (df['y1_t'] - beam_offset)**2 < 320) &
        (df['x1_d']**2 + (df['y1_d'] - beam_offset)**2 < 1100) &
        (df['x1_d']**2 + (df['y1_d'] - beam_offset)**2 > 8) &
        (df['chisq1_target'] < 1.5 * df['chisq1_upstream']) &
        (df['chisq1_target'] < 1.5 * df['chisq1_dump']) &
        (df['z1_v'] < -5) & (df['z1_v'] > -320) &
        (df['chisq1'] / (df['nHits1'] - 5) < 12) &
        (df['y1_st1'] / df['y1_st3'] < 1) &
        (((df['px1_st1'] - df['px1_st3']).abs() - 0.416).abs() < 0.008) &
        ((df['py1_st1'] - df['py1_st3']).abs() < 0.008) &
        ((df['pz1_st1'] - df['pz1_st3']).abs() < 0.08) &
        (df['y1_st1'] * df['y1_st3'] > 0) &
        (df['py1_st1'].abs() > 0.02)
    )

    # chuckCutsNegative_2111v42_tmp
    cut_neg = (
        (df['chisq2_target'] < 15) &
        (df['pz2_st1'] > 9) & (df['pz2_st1'] < 75) &
        (df['nHits2'] > 13) &
        (df['x2_t']**2 + (df['y2_t'] - beam_offset)**2 < 320) &
        (df['x2_d']**2 + (df['y2_d'] - beam_offset)**2 < 1100) &
        (df['x2_d']**2 + (df['y2_d'] - beam_offset)**2 > 8) &
        (df['chisq2_target'] < 1.5 * df['chisq2_upstream']) &
        (df['chisq2_target'] < 1.5 * df['chisq2_dump']) &
        (df['z2_v'] < -5) & (df['z2_v'] > -320) &
        (df['chisq2'] / (df['nHits2'] - 5) < 12) &
        (df['y2_st1'] / df['y2_st3'] < 1) &
        (((df['px2_st1'] - df['px2_st3']).abs() - 0.416).abs() < 0.008) &
        ((df['py2_st1'] - df['py2_st3']).abs() < 0.008) &
        ((df['pz2_st1'] - df['pz2_st3']).abs() < 0.08) &
        (df['y2_st1'] * df['y2_st3'] > 0) &
        (df['py2_st1'].abs() > 0.02)
    )

    # chuckCutsDimuon_2111v42
    cut_dimuon = (
        (df['dx'].abs() < 0.25) &
        ((df['dy'] - beam_offset).abs() < 0.22) &
        (df['dz'] > -280) & (df['dz'] < -5) &
        (df['dpx'].abs() < 1.8) &
        (df['dpy'].abs() < 2) &
        (df['dpx']**2 + df['dpy']**2 < 5) &
        (df['dpz'] > 38) & (df['dpz'] < 116) &
        (df['dx']**2 + (df['dy'] - beam_offset)**2 < 0.06) &
        (df['trackSeparation'].abs() < 270) &
        (df['chisq_dimuon'] < 18) &
        ((df['chisq1_target'] + df['chisq2_target'] - df['chisq_dimuon']).abs() < 2) &
        (df['y1_st3'] * df['y2_st3'] < 0) &
        (df['nHits1'] + df['nHits2'] > 29) &
        (df['nHits1St1'] + df['nHits2St1'] > 8) &
        ((df['x1_st1'] + df['x2_st1']).abs() < 42)
    )

    # physicsCuts_2111v42
    cut_phys = (
        (df['mass'] > 4.2) & (df['mass'] < 8.8) &
        (df['xF'] < 0.95) & (df['xF'] > -0.1) &
        (df['xT'] > 0.05) & (df['xT'] < 0.55) &
        (df['costh'].abs() < 0.5)
    )

    # occCuts_2111v42
    cut_occ = (
        (df['D1'] < 400) &
        (df['D2'] < 400) &
        (df['D3'] < 400) &
        (df['D1'] + df['D2'] + df['D3'] < 1000)
    )

    # Combine all cuts
    final_cut = cut_pos & cut_neg & cut_dimuon & cut_phys & cut_occ
    
    df_filtered = df[final_cut].copy() # Use .copy() to avoid SettingWithCopyWarning
    
    final_events = len(df_filtered)
    print(f"Initial events: {initial_events}. Events after dataCut: {final_events} ({final_events/initial_events:.2%} remaining).")
    
    return df_filtered


def calculate_average_efficiencies(df, npz_dir=".", var="D2"):
    """Calculates average efficiency for each (xF, Mass) bin.

    This function iterates through a predefined 2D grid of xF and mass bins.
    For each bin, it finds the corresponding .npz efficiency file, filters the
    input DataFrame for relevant events, and computes the average efficiency
    and associated errors.

    Args:
        df (pd.DataFrame): The input DataFrame containing dimuon events. Must
            include 'xF', 'mass', 'runID', and 'eventID' columns.
        npz_dir (str, optional): The directory where the .npz files are stored.
            Defaults to ".".
        var (str, optional): The name of the variable to use for efficiency
            calculation (e.g., 'D2'). Defaults to "D2".

    Returns:
        pd.DataFrame: A DataFrame containing the calculated average efficiencies
            and errors for each bin, including inverse values.
    """
    # 1. Define the binning structure to match the C++ script
    xf_bins = np.round(np.arange(0.0, 0.9, 0.05), 2)
    mass_bins = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7])

    results = [] # To store a dictionary of results from each bin

    # 2. Loop over both xF and Mass bin INDICES
    for ixf in range(len(xf_bins) - 1):
        for imass in range(len(mass_bins) - 1):

            # Get bin boundaries
            xf_low, xf_high = xf_bins[ixf], xf_bins[ixf + 1]
            mass_low, mass_high = mass_bins[imass], mass_bins[imass + 1]

            # 3. Construct the npz filename based on bin indices
            npz_file_for_bin = os.path.join(npz_dir, f"D2_Efficiency_xF{ixf}_mass{imass}.npz")

            print(f"--- Processing xF bin {ixf}, Mass bin {imass} using '{os.path.basename(npz_file_for_bin)}' ---")

            bin_results = {
                "xF_bin": f"[{xf_low}, {xf_high})",
                "mass_bin": f"[{mass_low}, {mass_high})",
                "N_events": 0,
                "avg_efficiency": np.nan,
                "stat_error": np.nan,
                "propagated_error": np.nan,
            }

            if not os.path.exists(npz_file_for_bin):
                print(f"Warning: File not found. Skipping.")
                results.append(bin_results)
                continue

            # 4. Filter the dataframe for events in the current bin
            df_cands = df[
                (df["xF"] >= xf_low) & (df["xF"] < xf_high) &
                (df["mass"] >= mass_low) & (df["mass"] < mass_high)
            ]

            df_cands = df_cands[~df_cands[["runID", "eventID"]].duplicated(keep="first")]
            df_cands = df_cands[df_cands[var] >= 8]

            if df_cands.empty:
                print("No events in this bin after cuts.")
                results.append(bin_results)
                continue

            # 5. Calculate efficiency for the events in this bin
            estimate, e_estimate = effi_model(
                df_cands[var], np.sqrt(df_cands[var]), npz_file_for_bin
            )

            eff_vals = estimate.iloc[:, 0].values

            # 6. Calculate averages and store them
            bin_results["N_events"] = len(df_cands)
            bin_results["avg_efficiency"] = np.mean(eff_vals)

            # Statistical error: sqrt( (<ε²> - <ε>²) / N )
            avg_eff_sq = np.mean(eff_vals**2)
            sq_avg_eff = bin_results["avg_efficiency"]**2
            bin_results["stat_error"] = np.sqrt(abs(avg_eff_sq - sq_avg_eff) / len(df_cands))

            # Propagated error
            bin_results["propagated_error"] = np.sqrt((e_estimate.iloc[:, 0]**2).sum()) / len(df_cands)

            results.append(bin_results)

    # 7. Convert the list of dictionaries to a final DataFrame
    df_results = pd.DataFrame(results)

    # 8. ADD NEW COLUMNS: Calculate inverse efficiency and its propagated error
    with np.errstate(divide='ignore', invalid='ignore'):
        # Calculate 1 / <ε>
        df_results['inv_avg_efficiency'] = 1.0 / df_results['avg_efficiency']

        # Propagate the error for 1 / <ε> which is δ(1/ε) = δε / ε²
        df_results['inv_propagated_error'] = (
            df_results['propagated_error'] / (df_results['avg_efficiency']**2)
        )

    # Replace infinite values resulting from division by zero with NaN
    df_results.replace([np.inf, -np.inf], np.nan, inplace=True)

    return df_results

# --- Main execution block ---
if __name__ == "__main__":
    # --- Configuration ---
    # Define the directory where your .npz files are located
    npz_dir = "/root/github/e906-development/src/kTrackerEfficiency/GenerateNPZFiles/D2_npz/" # Assumes a sub-directory named "D2_npz"

    # Define input and output files
    input_parquet_file = '/root/github/e906-development/ROOTFiles/Hugo/R008_roadset67_0_2111v42_tmp_noPhys_noOcc.parquet'
    output_csv_file = "average_efficiency_xF_mass_bins_RS67_LH2_only_targets_dataCut_dynamic_offset.csv"
    
    # --- Data Loading and Pre-selection ---
    print(f"Reading dimuon events from: {input_parquet_file}")
    df_events = pd.read_parquet(input_parquet_file)
    
    # Apply the full dataCut event selection with a dynamic BEAM_OFFSET
    df_events_filtered = apply_datacut(df_events)

    # Apply additional targetPos cut to only select LH2
    print("\nApplying target cut (targetPos == 1)...")
    initial_events = len(df_events_filtered)
    df_events_filtered = df_events_filtered[df_events_filtered["targetPos"] == 1]
    final_events = len(df_events_filtered)
    # Check for division by zero if initial_events is 0
    if initial_events > 0:
        print(f"Events after target cut: {final_events} ({final_events/initial_events:.2%} remaining).")
    else:
        print(f"Events after target cut: {final_events} (0 remaining).")


    # --- Main Calculation ---
    print("\nCalculating average efficiency for (xF, Mass) bins...")

    df_avg_efficiencies = calculate_average_efficiencies(df_events_filtered, npz_dir, var="D2")

    # --- Save Results ---
    # Reorder columns for clarity before saving
    column_order = [
        "xF_bin",
        "mass_bin",
        "N_events",
        "avg_efficiency",
        "stat_error",
        "propagated_error",
        "inv_avg_efficiency",
        "inv_propagated_error",
    ]
    df_avg_efficiencies = df_avg_efficiencies[column_order]

    df_avg_efficiencies.to_csv(output_csv_file, index=False)

    print("\n--- Calculation Complete ---")
    print("Results:")
    print(df_avg_efficiencies.to_string())
    print(f"\nResults have been successfully saved to: {output_csv_file}")
