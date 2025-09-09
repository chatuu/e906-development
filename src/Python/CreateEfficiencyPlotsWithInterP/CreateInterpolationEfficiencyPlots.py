import numpy as np
import matplotlib.pyplot as plt
import glob
import os

def create_efficiency_plots():
    """
    Finds all 'interpolation_data_*.npz' files in the current directory,
    and for each file, it generates a plot of 'Efficiency vs. D2' with
    asymmetric error bars and a linear interpolation curve.
    """
    # --- 1. Find all the .npz data files ---
    # The sorted() function ensures the files are processed in a consistent order.
    file_list = sorted(glob.glob('interpolation_data_*.npz'))

    # If no files are found, create dummy files for demonstration purposes.
    if not file_list:
        print("No 'interpolation_data_*.npz' files found.")
        print("Creating dummy .npz files for demonstration...")
        for i in range(5):
            low_bin = i * 0.1
            high_bin = low_bin + 0.05
            filename = f'interpolation_data_{low_bin:.2f}to{high_bin:.2f}.npz'
            # Generate some plausible-looking random data
            D2 = np.sort(np.random.rand(10) * 4 + 1) # D2 values from ~1 to 5
            Efficiency = 0.12 * D2 + 0.2 + (np.random.rand(10) - 0.5) * 0.15 # Linear trend with noise
            y_error_low = np.random.rand(10) * 0.08 # Random lower error
            y_error_high = np.random.rand(10) * 0.08 # Random upper error
            # Save dummy data with the keys the script expects ('x', 'y', etc.)
            np.savez(filename, x=D2, y=Efficiency, y_error_low=y_error_low, y_error_high=y_error_high)
            print(f"Created '{filename}'")
        # Update the file list with the newly created dummy files
        file_list = sorted(glob.glob('interpolation_data_*.npz'))

    # --- 2. Create an output directory for the plots ---
    output_dir = "efficiency_plots"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"\nCreated directory: '{output_dir}'")

    print("\nProcessing files and generating plots...")
    # --- 3. Loop through each file, load data, and create a plot ---
    for filepath in file_list:
        try:
            # Load the data arrays from the .npz file
            with np.load(filepath) as data:
                D2 = data['x']
                efficiency = data['y']
                y_error_low = data['y_error_low']
                y_error_high = data['y_error_high']

            # --- 4. Prepare data for plotting ---
            # np.interp requires x-coordinates to be sorted
            sort_indices = np.argsort(D2)
            D2_sorted = D2[sort_indices]
            efficiency_sorted = efficiency[sort_indices]
            
            # The errorbar function can handle asymmetric errors when passed a (2, N) array-like object.
            asymmetric_error = [y_error_low[sort_indices], y_error_high[sort_indices]]

            # --- 5. Create the plot ---
            fig, ax = plt.subplots(figsize=(12, 7))

            # Plot the original data points with their error bars
            # CHANGE: Updated colors for data points (blue) and error bars (black)
            ax.errorbar(D2_sorted, efficiency_sorted, yerr=asymmetric_error,
                        fmt='o', color='blue', ecolor='black', elinewidth=3,
                        capsize=5, markersize=8, label='Data Points with Error')

            # --- 6. Perform and plot linear interpolation ---
            # Create a dense set of D2 values for a smooth interpolation curve
            d2_interp = np.linspace(D2_sorted.min(), D2_sorted.max(), 500)
            
            # Calculate the corresponding efficiency values using linear interpolation
            efficiency_interp = np.interp(d2_interp, D2_sorted, efficiency_sorted)
            
            # Plot the interpolation line
            # CHANGE: Updated interpolation curve color to red
            ax.plot(d2_interp, efficiency_interp, '-', color='red', linewidth=2, label='Linear Interpolation')

            # --- 7. Customize plot aesthetics and labels ---
            # Extract xF bin from the filename to use in the title
            filename_stem = os.path.basename(filepath)
            # Correctly parse the filename string for the title
            xF_bin_str = filename_stem.replace('interpolation_data_', '').replace('.npz', '').replace('to', ' to ')
            
            ax.set_title(f'Efficiency vs D2 for xF bin: {xF_bin_str}', fontsize=16, weight='bold')
            ax.set_xlabel('D2', fontsize=14)
            ax.set_ylabel('Efficiency', fontsize=14)
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax.legend(fontsize=12)
            plt.tight_layout()

            # --- 8. Save the figure to the output directory ---
            # CHANGE: Save the output file in PDF format
            output_filename = f"plot_xF_{xF_bin_str.replace(' to ', '_')}.pdf"
            output_path = os.path.join(output_dir, output_filename)
            plt.savefig(output_path, dpi=300)
            print(f"  - Saved plot to {output_path}")

            # Close the plot to free up memory before the next loop iteration
            plt.close(fig)

        except Exception as e:
            print(f"Could not process file {filepath}. Error: {e}")

    print(f"\n✅ All plots have been generated and saved in the '{output_dir}' directory.")


if __name__ == '__main__':
    create_efficiency_plots()


