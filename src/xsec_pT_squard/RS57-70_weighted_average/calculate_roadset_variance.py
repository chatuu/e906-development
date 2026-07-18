import uproot
import numpy as np
import csv

def calculate_roadset_variance(target="LH2", suffix="_geom_pT2", use_equal_weights=False):
    """
    Computes the bin-by-bin standard deviation of cross sections across multiple roadsets.
    Assumes ROOT files were generated dynamically with exact, roadset-specific POTs via config.py.
    """
    
    rs_files = {
        "RS57": "XSec_RS57_Objects.root",
        "RS59": "XSec_RS59_Objects.root",
        "RS62": "XSec_RS62_Objects.root",
        "RS67": "XSec_RS67_Objects.root",
        "RS70": "XSec_RS70_Objects.root"
    }
    
    # Target the specific histogram inside the cross-section directory
    hist_path = f"CrossSections_{target}/h1_xsec_{target}{suffix}"
    
    cross_sections = {}
    stat_errors = {}
    
    for rs, filepath in rs_files.items():
        try:
            with uproot.open(filepath) as f:
                hist = f[hist_path]
                vals = np.array(hist.values())
                errs = np.array(hist.errors())
                
                for i in range(len(vals)):
                    if i not in cross_sections:
                        cross_sections[i] = []
                        stat_errors[i] = []
                    cross_sections[i].append(vals[i])
                    stat_errors[i].append(errs[i])
        except Exception as e:
            print(f"Error reading {rs} from {filepath}: {e}")
            return

    # Prepare CSV Output
    csv_filename = f"Roadset_Sys_StdDev_{target}{suffix}.csv"
    csv_rows = []
    
    print(f"\n--- Roadset Standard Deviation ({target}, {suffix}) ---")
    print(f"{'Bin':<5} | {'Mean xSec':<12} | {'Std Dev':<12} | {'% Rel Error':<10}")
    print("-" * 50)
    
    for bin_idx in sorted(cross_sections.keys()):
        x = np.array(cross_sections[bin_idx])
        e = np.array(stat_errors[bin_idx])
        
        # Filter out zero bins to prevent math domain errors on empty bins
        mask = (x > 0) & (e > 0)
        x_valid = x[mask]
        e_valid = e[mask]
        
        N = len(x_valid)
        if N < 2:
            csv_rows.append([bin_idx, 0.0, 0.0, 0.0])
            continue
            
        if use_equal_weights:
            mean_val = np.mean(x_valid)
            std_dev = np.std(x_valid, ddof=1)
        else:
            # Calculate standard deviation using statistical inverse variance weighting
            w = 1.0 / (e_valid**2)
            V1 = np.sum(w)
            V2 = np.sum(w**2)
            mean_val = np.average(x_valid, weights=w)
            variance = (V1 / (V1**2 - V2)) * np.sum(w * (x_valid - mean_val)**2)
            std_dev = np.sqrt(variance)
            
        rel_err = (std_dev / mean_val) * 100 if mean_val > 0 else 0
        
        # Using scientific notation (e) for cross-sections to ensure smaller pT2 values format cleanly
        print(f"{bin_idx:<5} | {mean_val:<12.6e} | {std_dev:<12.6e} | {rel_err:.2f}%")
        csv_rows.append([bin_idx, mean_val, std_dev, rel_err])

    with open(csv_filename, "w", newline="") as f:
        writer = csv.writer(f)
        weight_type = "Equal Weights" if use_equal_weights else "Stat Inverse Weights"
        writer.writerow([f"Roadset Std Dev - {target} {suffix} - {weight_type}"])
        writer.writerow(["Bin Index", "Mean Cross Section", "Standard Deviation", "Relative Sys Error (%)"])
        writer.writerows(csv_rows)
        
    print(f"\nSaved variance table to {csv_filename}")

if __name__ == "__main__":
    # Execute the variance calculation for geometric pT2 bins for both targets
    calculate_roadset_variance("LH2", "_geom_pT2", use_equal_weights=False)
    calculate_roadset_variance("LD2", "_geom_pT2", use_equal_weights=False)