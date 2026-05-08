import os
import uproot
import numpy as np
import ROOT
import array

class EventNamespace:
    """Helper class to allow dot-notation access to dictionary keys."""
    def __init__(self, data):
        self.__dict__.update(data)

def apply_cuts(tree, cut=0.0):
    """Applies standard physics cuts to the given TTree arrays in numpy."""
    events = tree.arrays(library="np")
    e = EventNamespace(events)

    # Dynamic beam offset based on runID
    bo = np.where(e.runID >= 11000, 1.6, 0.4)

    dimuon_cut = (
        (np.abs(e.dx) < 0.25) & (np.abs(e.dy - bo) < 0.22) &
        (e.dz < -5.) & (e.dz > -280.) & (np.abs(e.dpx) < 1.8) & (np.abs(e.dpy) < 2.0) &
        (e.dpx * e.dpx + e.dpy * e.dpy < 5.) & (e.dpz < 116.) & (e.dpz > 38.) &
        (e.mass > cut) & (e.mass < 8.8) &
        (e.dx * e.dx + (e.dy - bo) * (e.dy - bo) < 0.06) &
        (e.xF < 0.95) & (e.xF > -0.1) & (e.xT > 0.05) & (e.xT <= 0.58) &
        (np.abs(e.costh) < 0.5) & (np.abs(e.trackSeparation) < 270.) &
        (e.chisq_dimuon < 18)
    )

    track1_cut = (
        (e.chisq1_target < 15.) & (e.pz1_st1 > 9.) & (e.pz1_st1 < 75.) & (e.nHits1 > 13) &
        (e.x1_t * e.x1_t + (e.y1_t - bo) * (e.y1_t - bo) < 320.) &
        (e.x1_d * e.x1_d + (e.y1_d - bo) * (e.y1_d - bo) < 1100.) &
        (e.x1_d * e.x1_d + (e.y1_d - bo) * (e.y1_d - bo) > 16.) &
        (e.chisq1_target < 1.5 * e.chisq1_upstream) & (e.chisq1_target < 1.5 * e.chisq1_dump) &
        (e.z1_v < -5.) & (e.z1_v > -320.) & (e.chisq1 / (e.nHits1 - 5) < 12) &
        ((e.y1_st1) / (e.y1_st3) < 1.) & (np.abs(np.abs(e.px1_st1 - e.px1_st3) - 0.416) < 0.008) &
        (np.abs(e.py1_st1 - e.py1_st3) < 0.008) & (np.abs(e.pz1_st1 - e.pz1_st3) < 0.08) &
        ((e.y1_st1) * (e.y1_st3) > 0.) & (np.abs(e.py1_st1) > 0.02)
    )

    track2_cut = (
        (e.chisq2_target < 15.) & (e.pz2_st1 > 9.) & (e.pz2_st1 < 75.) & (e.nHits2 > 13) &
        (e.x2_t * e.x2_t + (e.y2_t - bo) * (e.y2_t - bo) < 320.) &
        (e.x2_d * e.x2_d + (e.y2_d - bo) * (e.y2_d - bo) < 1100.) &
        (e.x2_d * e.x2_d + (e.y2_d - bo) * (e.y2_d - bo) > 16.) &
        (e.chisq2_target < 1.5 * e.chisq2_upstream) & (e.chisq2_target < 1.5 * e.chisq2_dump) &
        (e.z2_v < -5.) & (e.z2_v > -320.) & (e.chisq2 / (e.nHits2 - 5) < 12) &
        ((e.y2_st1) / (e.y2_st3) < 1.) & (np.abs(np.abs(e.px2_st1 - e.px2_st3) - 0.416) < 0.008) &
        (np.abs(e.py2_st1 - e.py2_st3) < 0.008) & (np.abs(e.pz2_st1 - e.pz2_st3) < 0.08) &
        ((e.y2_st1) * (e.y2_st3) > 0.) & (np.abs(e.py2_st1) > 0.02)
    )

    tracks_cut = (
        (np.abs(e.chisq1_target + e.chisq2_target - e.chisq_dimuon) < 2.) &
        ((e.y1_st3) * (e.y2_st3) < 0.) & (e.nHits1 + e.nHits2 > 29) &
        (e.nHits1St1 + e.nHits2St1 > 8) & (np.abs(e.x1_st1 + e.x2_st1) < 42)
    )

    occ_cut = (
        (e.D1 < 400) & (e.D2 < 400) & (e.D3 < 400) & (e.D1 + e.D2 + e.D3 < 1000)
    )

    total_cut_mask = (track1_cut & track2_cut & tracks_cut & dimuon_cut & occ_cut)

    filtered_events = {}
    for key, val in events.items():
        filtered_events[key] = val[total_cut_mask]
        
    return filtered_events

def get_concatenated_events(file_paths, tree_name):
    """Iterates through file paths, applies cuts, and merges the numpy arrays in memory."""
    all_events = None
    for fp in file_paths:
        if not os.path.exists(fp):
            print(f"Warning: File {fp} does not exist, skipping.")
            continue
        try:
            with uproot.open(fp) as f:
                if tree_name not in f:
                    print(f"Warning: Tree '{tree_name}' not found in {fp}, skipping.")
                    continue
                filtered = apply_cuts(f[tree_name])
                if all_events is None:
                    all_events = {k: [v] for k, v in filtered.items()}
                else:
                    for k, v in filtered.items():
                        if k in all_events:
                            all_events[k].append(v)
        except Exception as e:
            print(f"Error reading {fp}: {e}")
            
    if all_events is None:
        raise RuntimeError(f"No valid data found for tree '{tree_name}' in provided files.")
    
    return {k: np.concatenate(v) for k, v in all_events.items()}

def fill_root_histogram(hist, data_array):
    """Efficiently fills a ROOT TH1 from a NumPy array."""
    data_list = np.asarray(data_array).astype('float64')
    weights = np.ones(len(data_list), dtype='float64')
    hist.FillN(len(data_list), array.array('d', data_list), array.array('d', weights))

def process_target_canvas(target_name, file_paths, var_name, bins, x_min, x_max, output_root_file):
    """Reads multiple trees, applies cuts, draws to TCanvas, and saves to ROOT file directory."""
    print(f"Processing {target_name} target across {len(file_paths)} files...")
    
    # Calculate bin width for the Y-axis label
    bin_width = (x_max - x_min) / bins
    
    # Extract, cut, and concatenate events across all files
    events_tot = get_concatenated_events(file_paths, "result")
    events_mix = get_concatenated_events(file_paths, "result_mix")

    # Extract target variable array
    data_tot = events_tot[var_name]
    data_mix = events_mix[var_name]

    # Initialize ROOT Histograms
    h_tot = ROOT.TH1D(f"Y_total_{target_name}", f"{target_name} Yields", bins, x_min, x_max)
    h_mix = ROOT.TH1D(f"Y_mix_{target_name}", f"{target_name} Y_{{mix}}", bins, x_min, x_max)

    # Calculate statistical errors
    h_tot.Sumw2()
    h_mix.Sumw2()

    # Fill histograms
    fill_root_histogram(h_tot, data_tot)
    fill_root_histogram(h_mix, data_mix)

    # Calculate difference
    h_diff = h_tot.Clone(f"Y_diff_{target_name}")
    h_diff.Add(h_mix, -1.0)

    # Configure Axis Titles & Centering (Done on the first drawn histogram)
    h_tot.GetXaxis().SetTitle("Invariant Mass (GeV)")
    h_tot.GetXaxis().CenterTitle(True)
    h_tot.GetYaxis().SetTitle(f"Counts per {bin_width:.2f} GeV")
    h_tot.GetYaxis().CenterTitle(True)

    # Set up TCanvas with required ticks (no grid)
    canvas = ROOT.TCanvas(f"Canvas_{target_name}", f"{target_name} Analysis", 800, 600)
    canvas.SetTickx(1)
    canvas.SetTicky(1)

    # Styling
    h_tot.SetLineColor(ROOT.kBlack)
    h_tot.SetMarkerColor(ROOT.kBlack)
    h_tot.SetMarkerStyle(20)

    h_mix.SetLineColor(ROOT.kBlue)
    h_mix.SetMarkerColor(ROOT.kBlue)
    h_mix.SetMarkerStyle(21)

    h_diff.SetLineColor(ROOT.kRed)
    h_diff.SetMarkerColor(ROOT.kRed)
    h_diff.SetMarkerStyle(22)

    # Draw everything
    h_tot.Draw("E1")
    h_mix.Draw("E1 SAME")
    h_diff.Draw("E1 SAME")

    # Add Updated Legend (Coordinates adjusted slightly to fit longer text)
    legend = ROOT.TLegend(0.45, 0.70, 0.88, 0.88)
    legend.AddEntry(h_tot, "Data (Signal + Backgrounds)", "lep")
    legend.AddEntry(h_mix, "Background Estimate", "lep")
    legend.AddEntry(h_diff, "Extracted Signal (DY, J/#psi, #psi')", "lep")
    legend.SetBorderSize(0)
    legend.Draw()

    # Save to TDirectory in the output ROOT file
    if output_root_file:
        target_dir = output_root_file.mkdir(target_name)
        target_dir.cd()
        h_tot.Write()
        h_mix.Write()
        h_diff.Write()
        canvas.Write()

    # Prevent garbage collection
    ROOT.SetOwnership(canvas, False)
    ROOT.SetOwnership(h_tot, False)
    ROOT.SetOwnership(h_mix, False)
    ROOT.SetOwnership(h_diff, False)
    ROOT.SetOwnership(legend, False)

    return canvas

def main():
    # Globally disable statistics boxes for all plots
    ROOT.gStyle.SetOptStat(0)

    base_dir = "/root/github/e906-development/ROOTFiles/Harsha/ROOTFilesWithoutRoadIDAdded"
    mixed_dir = "/root/github/e906-development/ROOTFiles/MixedEvents"

    # Input File paths for Roadsets 57, 59, 62, 67, 70
    lh2_files = [
        f"{base_dir}/merged_RS57_LH2_1_1138.root",
        f"{base_dir}/merged_RS59_LH2_1_465.root",
        f"{base_dir}/merged_RS62_LH2_1_1234.root",
        f"{mixed_dir}/merged_RS67_3089LH2.root",
        f"{base_dir}/merged_RS70_LH2_1_264.root"
    ]

    ld2_files = [
        f"{base_dir}/merged_RS57_LD2_3_1138.root",
        f"{base_dir}/merged_RS59_LD2_3_466.root",
        f"{base_dir}/merged_RS62_LD2_3_1237.root",
        f"{mixed_dir}/merged_RS67_3089LD2.root",
        f"{base_dir}/merged_RS70_LD2_3_266.root"
    ]
    
    # Output ROOT File
    output_filename = "Target_Yields_Analysis.root"
    out_file = ROOT.TFile(output_filename, "RECREATE")

    # Plot parameters
    variable_to_plot = "mass"
    bins = 70     # Adjusted to 70 for an exact 0.10 GeV bin width across a 7.0 GeV range
    x_min = 2.0   # Updated lower bound
    x_max = 9.0   # Updated upper bound

    # Generate canvases & write root file contents using the file lists
    canvas_lh2 = process_target_canvas("LH2", lh2_files, variable_to_plot, bins, x_min, x_max, out_file)
    canvas_ld2 = process_target_canvas("LD2", ld2_files, variable_to_plot, bins, x_min, x_max, out_file)

    # Save to PDF files
    canvas_lh2.SaveAs("LH2_Yield_Comparison.pdf")
    canvas_ld2.SaveAs("LD2_Yield_Comparison.pdf")

    # Close output root file
    out_file.Close()

    print(f"Processing complete. Canvases saved as PDFs. All TObjects saved to {output_filename}")

if __name__ == "__main__":
    main()