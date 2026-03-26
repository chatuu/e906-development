import numpy as np
import uproot
import awkward as ak
import ROOT
from scipy.interpolate import interp1d
import argparse
import os

# ==========================================
# ROOT Configuration
# ==========================================
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

# ==========================================
# CONFIGURATION
# ==========================================
# Fixed path for the interpolation data
input_npz_file  = "../GlobalEfficiencyCurve/interpolation_data_d1.npz"

# --- TOGGLES ---
# Set to False to skip covariance matrix generation and speed up the script
GENERATE_COVAR_MATRIX = True

# Limit the covariance matrix size to prevent Out-Of-Memory errors
max_covar_samples_per_bin = 5000 

# Binning Definitions (Match your analysis)
mass_bins_np = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
xf_bins_np = np.round(np.arange(0.0, 0.85, 0.05), 2)

# ==========================================
# EVENT SELECTION (CHUCK CUTS)
# ==========================================
def e906_chuck_cuts(tree: uproot.models.TTree.Model_TTree_v19, cut=4.2, beam_offset: float = 1.6):
    branch = tree.keys()
    events = tree.arrays(branch)

    dimuon_cut_2111_v42 = (
        (np.abs(events.dx) < 0.25) &
        (np.abs(events.dy - beam_offset) < 0.22) &
        (events.dz < -5.) &
        (events.dz > -280.) &
        (np.abs(events.dpx) < 1.8) &
        (np.abs(events.dpy) < 2.0) &
        (events.dpx * events.dpx + events.dpy * events.dpy < 5.) &
        (events.dpz < 116.) &
        (events.dpz > 38.) &
        (events.mass > cut) &
        (events.mass < 8.8) &
        (events.dx * events.dx + (events.dy - beam_offset) * (events.dy - beam_offset) < 0.06) &
        (events.xF < 0.95) &
        (events.xF > -0.1) &
        (events.xT > 0.05) &
        (events.xT <= 0.58) &
        (np.abs(events.costh) < 0.5) &
        (np.abs(events.trackSeparation) < 270.) &
        (events.chisq_dimuon < 18)
    )

    track1_cut_2111_v42 = (
        (events.chisq1_target < 15.) &
        (events.pz1_st1 > 9.) &
        (events.pz1_st1 < 75.) &
        (events.nHits1 > 13) &
        (events.x1_t * events.x1_t + (events.y1_t - beam_offset) * (events.y1_t - beam_offset) < 320.) &
        (events.x1_d * events.x1_d + (events.y1_d - beam_offset) * (events.y1_d -beam_offset) < 1100.) &
        (events.x1_d * events.x1_d + (events.y1_d - beam_offset) * (events.y1_d -beam_offset) > 16.) &
        (events.chisq1_target < 1.5 * events.chisq1_upstream) &
        (events.chisq1_target < 1.5 * events.chisq1_dump) &
        (events.z1_v < -5.) &
        (events.z1_v > -320.) &
        (events.chisq1/(events.nHits1 - 5) < 12) &
        ((events.y1_st1)/(events.y1_st3 ) < 1.) & 
        (np.abs(np.abs(events.px1_st1 - events.px1_st3) - 0.416) < 0.008) &
        (np.abs(events.py1_st1 - events.py1_st3) < 0.008) &
        (np.abs(events.pz1_st1 - events.pz1_st3) < 0.08) &
        ((events.y1_st1) * (events.y1_st3) > 0.) & 
        (np.abs(events.py1_st1) > 0.02)
    )

    track2_cut_2111_v42 = (
        (events.chisq2_target < 15.) &
        (events.pz2_st1 > 9.) &
        (events.pz2_st1 < 75.) &
        (events.nHits2 > 13) &
        (events.x2_t * events.x2_t + (events.y2_t - beam_offset) * (events.y2_t - beam_offset) < 320.) &
        (events.x2_d * events.x2_d + (events.y2_d - beam_offset) * (events.y2_d -beam_offset) < 1100.) &
        (events.x2_d * events.x2_d + (events.y2_d - beam_offset) * (events.y2_d -beam_offset) > 16.) &
        (events.chisq2_target < 1.5 * events.chisq2_upstream) &
        (events.chisq2_target < 1.5 * events.chisq2_dump) &
        (events.z2_v < -5.) &
        (events.z2_v > -320.) &
        (events.chisq2/(events.nHits2 - 5) < 12) &
        ((events.y2_st1 )/(events.y2_st3) < 1.) & 
        (np.abs(np.abs(events.px2_st1 - events.px2_st3) - 0.416) < 0.008) &
        (np.abs(events.py2_st1 - events.py2_st3) < 0.008) &
        (np.abs(events.pz2_st1 - events.pz2_st3) < 0.08) &
        ((events.y2_st1) * (events.y2_st3) > 0.) & 
        (np.abs(events.py2_st1) > 0.02)
    )

    tracks_cut_2111_v42 = (
        (np.abs(events.chisq1_target + events.chisq2_target - events.chisq_dimuon) < 2.) &
        ((events.y1_st3) * (events.y2_st3) < 0.) & 
        (events.nHits1 + events.nHits2 > 29) &
        (events.nHits1St1 + events.nHits2St1 > 8) &
        (np.abs(events.x1_st1 + events.x2_st1) < 42)
    )

    occ_cut_2111_v42 = (
        (events.D1 < 400) &
        (events.D2 < 400) &
        (events.D3 < 400) &
        (events.D1 + events.D2 + events.D3 < 1000)
    )

    events_cut = events[track1_cut_2111_v42 & track2_cut_2111_v42 & tracks_cut_2111_v42 & dimuon_cut_2111_v42 & occ_cut_2111_v42]
    
    return events_cut

# ==========================================
# ERROR PROPAGATION & CORRELATION LOGIC
# ==========================================
def get_bounds_vectorized(vals, ref_array):
    vals = np.asarray(vals)
    ref_array = np.asarray(ref_array)
    diff = vals[:, None] - ref_array[None, :]
    
    mask_lower = (diff >= 0)
    lower_candidates = np.where(mask_lower, ref_array[None, :], -np.inf)
    idx_lower = np.argmax(lower_candidates, axis=1)
    val_lower = ref_array[idx_lower]
    no_lower = np.all(~mask_lower, axis=1)
    idx_lower[no_lower] = 0 
    val_lower[no_lower] = ref_array[0]

    mask_upper = (diff <= 0)
    upper_candidates = np.where(mask_upper, ref_array[None, :], np.inf)
    idx_upper = np.argmin(upper_candidates, axis=1)
    val_upper = ref_array[idx_upper]
    no_upper = np.all(~mask_upper, axis=1)
    idx_upper[no_upper] = len(ref_array) - 1
    val_upper[no_upper] = ref_array[-1]
    
    return val_lower, idx_lower, val_upper, idx_upper

def calculate_recoeff_and_error(d1_vals, x_curve, y_curve, y_err_low, y_err_high):
    d1_vals = np.array(d1_vals)
    f_linear = interp1d(x_curve, y_curve, kind='linear', fill_value="extrapolate")
    track_effi = f_linear(d1_vals)
    
    d_minus, idx_d_minus, d_plus, idx_d_plus = get_bounds_vectorized(d1_vals, x_curve)
    effi_minus_val, idx_effi_minus, effi_plus_val, idx_effi_plus = get_bounds_vectorized(track_effi, y_curve)
    
    e_effi_minus = y_err_low[idx_effi_minus]
    e_effi_plus  = y_err_high[idx_effi_plus]
    
    d2 = d1_vals 
    delta_d = d_plus - d_minus
    mask_exact = (delta_d == 0)
    recoeff_error = np.zeros_like(track_effi)
    safe_delta_d = np.where(mask_exact, 1.0, delta_d) 
    
    de_minus = 1 + d2 * (1 / safe_delta_d)
    de_plus  = -d2 * (1 / safe_delta_d)
    dd2      = -(effi_plus_val - effi_minus_val) / safe_delta_d
    dd_plus  = d2 * (effi_plus_val - effi_minus_val) / (safe_delta_d**2)
    dd_minus = -d2 * (effi_plus_val - effi_minus_val) / (safe_delta_d**2)
    
    variance = (
        (de_minus**2) * (e_effi_minus**2) +
        (de_plus**2)  * (e_effi_plus**2) +
        (dd2**2)      * d2 +
        (dd_plus**2)  * d_plus +
        (dd_minus**2) * d_minus
    )
    
    recoeff_error[~mask_exact] = np.sqrt(variance[~mask_exact])
    recoeff_error[mask_exact] = e_effi_minus[mask_exact]
    
    return track_effi, recoeff_error, idx_d_minus

def generate_covariance_matrix_per_bin(events_cut, loc_data_array, recoeff_error_array, target_name, tree_name):
    """
    Loops through double-differential bins, extracts events inside that bin,
    and generates/saves a specific covariance matrix as a TH2D.
    """
    # Map the tree name to 'total' or 'mix'
    tree_label = "mix" if tree_name == "result_mix" else "total"
    
    out_dir = f"CovarianceMatrices_{target_name}_{tree_label}"
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"    -> Generating per-bin Covariance Matrices (Saving to '{out_dir}/')...")
    
    mass_vals = np.asarray(events_cut.mass)
    xf_vals = np.asarray(events_cut.xF)
    
    inner_correl = 1.0
    neighbor_correl = 1.0

    bins_processed = 0
    
    for i_m in range(len(mass_bins_np) - 1):
        m_low, m_high = mass_bins_np[i_m], mass_bins_np[i_m+1]
        
        for i_x in range(len(xf_bins_np) - 1):
            x_low, x_high = xf_bins_np[i_x], xf_bins_np[i_x+1]
            
            # Mask to find events in this specific kinematic bin
            mask = (mass_vals >= m_low) & (mass_vals < m_high) & (xf_vals >= x_low) & (xf_vals < x_high)
            
            locs_in_bin = loc_data_array[mask]
            errs_in_bin = recoeff_error_array[mask]
            
            N = len(locs_in_bin)
            if N == 0:
                continue 
                
            bins_processed += 1
            n_samples = min(N, max_covar_samples_per_bin)
            
            locs = locs_in_bin[:n_samples]
            errs = errs_in_bin[:n_samples]
            
            # Create correlation matrix logic
            diff_matrix = np.abs(locs[:, None] - locs[None, :])
            correl_matrix = np.zeros((n_samples, n_samples))
            correl_matrix[diff_matrix == 0] = inner_correl      
            correl_matrix[diff_matrix == 1] = neighbor_correl   
        
            # Generate final Covariance Matrix in NumPy
            covar_matrix = correl_matrix * np.outer(errs, errs)
            
            # Create ROOT TH2D
            hist_name = f"covar_{target_name}_{tree_label}_xf{i_x}_m{i_m}"
            hist_title = f"{target_name} ({tree_label}) Covariance | M: [{m_low:.1f}, {m_high:.1f}) xF: [{x_low:.2f}, {x_high:.2f});Event i;Event j"
            
            h_cov = ROOT.TH2D(hist_name, hist_title, n_samples, 0, n_samples, n_samples, 0, n_samples)
            h_cov.SetStats(0)
            h_cov.GetXaxis().CenterTitle()
            h_cov.GetYaxis().CenterTitle()
            h_cov.GetXaxis().SetTitleOffset(1.2)
            h_cov.GetYaxis().SetTitleOffset(1.2)

            # Highly optimized filling of sparse ROOT TH2D 
            nonzero_i, nonzero_j = np.nonzero(covar_matrix)
            for i, j in zip(nonzero_i, nonzero_j):
                # TH2D indices are 1-based in ROOT
                h_cov.SetBinContent(int(i) + 1, int(j) + 1, float(covar_matrix[i, j]))

            # Setup Canvas and save as PDF
            c = ROOT.TCanvas(f"c_{hist_name}", "", 800, 800)
            c.SetRightMargin(0.15) # Leave space for COLZ palette
            
            h_cov.Draw("COLZ")
            
            pdf_filename = f"{out_dir}/CovarianceMatrix_{target_name}_{tree_label}_XfBin_{i_x}_MassBin_{i_m}.pdf"
            c.SaveAs(pdf_filename)
            c.Close()
            
    print(f"    -> Finished. Saved {bins_processed} covariance matrices to {out_dir}/.")

# ==========================================
# MAIN EXECUTION
# ==========================================
def process_file(input_root_file, output_root_file, target_name):
    print(f"\n==========================================")
    print(f"Processing Target: {target_name}")
    print(f"Input File: {input_root_file}")
    print(f"Output File: {output_root_file}")
    print(f"==========================================\n")

    try:
        print(f"Loading .npz data from {input_npz_file}...")
        npz_data = np.load(input_npz_file)
        x_curve = npz_data['x']
        y_curve = npz_data['y']
        y_err_low = npz_data['y_error_low']
        y_err_high = npz_data['y_error_high']
        
        print(f"Opening input ROOT file: {input_root_file}")
        with uproot.open(input_root_file) as file:
            tree_names = ['result', 'result_mix']
            output_trees = {}
            
            for t_name in tree_names:
                if t_name not in file:
                    print(f"Warning: Tree '{t_name}' not found. Skipping.")
                    continue
                
                print(f"\nProcessing tree: {t_name}...")
                tree = file[t_name]
                
                # 1. Apply cuts first
                events_cut = e906_chuck_cuts(tree)
                print(f"  - Events after cuts: {len(events_cut)}")
                
                if len(events_cut) == 0:
                    print("  - No events passed cuts. Skipping efficiency calculations.")
                    continue
                
                # 2. Extract D1 values
                if t_name == 'result_mix':
                    if 'ptrk_D1' in events_cut.fields and 'ntrk_D1' in events_cut.fields:
                        d1_values = 0.5 * (events_cut.ptrk_D1 + events_cut.ntrk_D1)
                    else:
                        d1_values = events_cut.D1
                else:
                    d1_values = events_cut.D1

                # 3. Calculate efficiencies
                print(f"  - Calculating efficiencies...")
                recoeff, recoeff_error, loc_data = calculate_recoeff_and_error(
                    np.asarray(d1_values), 
                    x_curve, y_curve, y_err_low, y_err_high
                )
                
                events_cut["recoeff"] = recoeff
                events_cut["recoeff_error"] = recoeff_error
                
                output_trees[t_name] = {field: events_cut[field] for field in events_cut.fields}

                # 4. Generate the covariance matrices per bin
                if GENERATE_COVAR_MATRIX:
                    generate_covariance_matrix_per_bin(events_cut, loc_data, recoeff_error, target_name, t_name)
                else:
                    print("    -> Skipping Covariance Matrix generation (Flag is set to False).")

        # Write Output
        if output_trees:
            print(f"\nWriting output to {output_root_file}...")
            with uproot.recreate(output_root_file) as f_out:
                for t_name, data_dict in output_trees.items():
                    print(f"  - Writing tree {t_name}...")
                    f_out[t_name] = data_dict
            print("Done! File saved successfully.")
        else:
            print("\nNo trees were processed or survived the cuts.")
            
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate ROOT Files with Reco Efficiencies")
    parser.add_argument("-i", "--input", required=True, help="Path to the input ROOT file")
    parser.add_argument("-o", "--output", required=True, help="Path to save the output ROOT file")
    parser.add_argument("-t", "--target", required=True, choices=["LH2", "LD2", "Flask"], help="Target name (LH2, LD2, Flask)")
    
    args = parser.parse_args()
    
    process_file(args.input, args.output, args.target)