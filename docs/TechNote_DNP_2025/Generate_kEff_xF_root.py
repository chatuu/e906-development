import ROOT
import pandas as pd
import numpy as np
import os

def create_efficiency_histograms(csv_file, root_file_name):
    """
    Reads an efficiency CSV file and creates a ROOT file containing TH1D histograms.

    For each unique xF bin in the CSV, this function creates a TH1D histogram binned
    in mass. The bin contents are set to the inverse average efficiency, and the bin
    errors are set to the corresponding inverse propagated error.

    A condition is applied: if the 'avg_efficiency' for a bin is 0 or greater than 1,
    the bin content and error for that bin are set to 0.

    Args:
        csv_file (str): Path to the input CSV file.
        root_file_name (str): Name of the output ROOT file to be created.
    """
    # --- 1. Define Binning and Load Data ---

    # Define the variable-width mass bins for the TH1D histograms
    mass_bins = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
    num_mass_bins = len(mass_bins) - 1

    print(f"Reading data from '{csv_file}'...")
    if not os.path.exists(csv_file):
        print(f"Error: Input file not found at '{csv_file}'")
        return

    df = pd.read_csv(csv_file)

    # --- 2. Clean and Prepare Data ---

    # Columns to be converted to numeric types.
    # The 'coerce' option will replace non-numeric values (like empty strings) with NaN.
    numeric_cols = ['avg_efficiency', 'inv_avg_efficiency', 'inv_propagated_error']
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Replace any resulting NaN values with 0.0, which is appropriate for bins with no data.
    df.fillna(0.0, inplace=True)

    # Get the unique xF bins from the data and sort them to ensure a consistent processing order.
    xf_bins_unique = sorted(df['xF_bin'].unique())
    print(f"Found {len(xf_bins_unique)} unique xF bins.")

    # --- 3. Create and Fill ROOT Histograms ---

    # Create the output ROOT file in RECREATE mode (overwrites if it exists).
    output_file = ROOT.TFile(root_file_name, "RECREATE")
    if not output_file.IsOpen():
        print(f"Error: Could not create ROOT file '{root_file_name}'")
        return

    print("Creating and filling histograms...")

    # Iterate over each unique xF bin to create a separate histogram.
    for ixf, xf_bin_str in enumerate(xf_bins_unique):
        # Define a unique name and title for each histogram.
        hist_name = f"kEff_xF_{ixf}"
        hist_title = f"Inverse Efficiency vs Mass for xF in {xf_bin_str}"

        # Create the TH1D histogram with the specified variable mass bins.
        hist = ROOT.TH1D(hist_name, hist_title, num_mass_bins, mass_bins)

        # Filter the DataFrame to get the data only for the current xF bin.
        df_subset = df[df['xF_bin'] == xf_bin_str]

        # Iterate over the rows in the subset. Each row corresponds to one mass bin.
        for _, row in df_subset.iterrows():
            # Find the correct bin index in the histogram for the current mass value.
            # We parse the lower edge of the mass bin string (e.g., "[4.2, 4.5)")
            mass_low_str = row['mass_bin'].split(',')[0].replace('[', '')
            mass_low = float(mass_low_str)
            
            # FindBin returns the global bin number for the given coordinate.
            # An epsilon is added to handle values exactly on the bin edge.
            bin_index = hist.FindBin(mass_low + 0.01)

            avg_eff = row['avg_efficiency']
            inv_avg_eff = row['inv_avg_efficiency']
            inv_prop_err = row['inv_propagated_error']

            # --- 4. Apply Condition and Set Bin Content/Error ---

            # If avg_efficiency is unphysical (0 or >1), set content and error to 0.
            if avg_eff == 0.0 or avg_eff > 1.0:
                print(f"Average Efficiency hits zero at the bin: {bin_index} at xF bin: {ixf}")
                hist.SetBinContent(bin_index, 0.0)
                hist.SetBinError(bin_index, 0.0)
            else:
                hist.SetBinContent(bin_index, inv_avg_eff)
                hist.SetBinError(bin_index, inv_prop_err)

        # Set axis titles for better readability.
        hist.GetXaxis().SetTitle("Mass (GeV)")
        hist.GetYaxis().SetTitle("1 / <#epsilon>")

        # Write the completed histogram to the ROOT file.
        hist.Write()

    # --- 5. Finalize ---
    output_file.Close()
    print(f"\nSuccessfully created ROOT file '{root_file_name}' with {len(xf_bins_unique)} histograms.")

# --- Main execution block ---
if __name__ == "__main__":
    # Define the name for the input CSV file.
    input_csv = "average_efficiency_xF_mass_bins_RS67_LH2_only_targets_dataCut_dynamic_offset.csv"
    
    # Create the CSV file from the multiline string data
    #with open(input_csv, "w") as f:
    #    f.write(csv_data)

    # Define the name for the output ROOT file
    output_root = "kEff_xF.root"

    # Execute the function to generate the ROOT file
    create_efficiency_histograms(input_csv, output_root)
