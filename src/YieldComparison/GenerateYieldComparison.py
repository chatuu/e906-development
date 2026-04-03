import uproot
import numpy as np
import ROOT
import array

def e906_chuck_cuts(tree: uproot.models.TTree.Model_TTree_v19, cut=4.2, beam_offset: float = 1.6):
    """Event selection criteria."""
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

def fill_root_histogram(hist, data_array):
    """Efficiently fills a ROOT TH1 from a NumPy array."""
    data_list = np.asarray(data_array).astype('float64')
    weights = np.ones(len(data_list), dtype='float64')
    hist.FillN(len(data_list), array.array('d', data_list), array.array('d', weights))

def process_target_canvas(target_name, file_path, var_name, bins, x_min, x_max, output_root_file):
    """Reads trees, applies cuts, draws to TCanvas, and saves to ROOT file directory."""
    print(f"Processing {target_name} target...")
    
    # Calculate bin width for the Y-axis label
    bin_width = (x_max - x_min) / bins
    
    # Open trees
    tree_tot = uproot.open(f"{file_path}:result")
    tree_mix = uproot.open(f"{file_path}:result_mix")

    # Apply selection criteria
    events_tot = e906_chuck_cuts(tree_tot)
    events_mix = e906_chuck_cuts(tree_mix)

    # Extract target variable
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

    # Input File paths
    file_lh2 = "/root/github/e906-development/ROOTFiles/MixedEvents/merged_RS67_3089LH2.root"
    file_ld2 = "/root/github/e906-development/ROOTFiles/MixedEvents/merged_RS67_3089LD2.root"
    
    # Output ROOT File
    output_filename = "Target_Yields_Analysis.root"
    out_file = ROOT.TFile(output_filename, "RECREATE")

    # Plot parameters
    variable_to_plot = "mass"
    bins = 70     # Adjusted to 70 for an exact 0.10 GeV bin width across a 7.0 GeV range
    x_min = 2.0   # Updated lower bound
    x_max = 9.0   # Updated upper bound

    # Generate canvases & write root file contents
    canvas_lh2 = process_target_canvas("LH2", file_lh2, variable_to_plot, bins, x_min, x_max, out_file)
    canvas_ld2 = process_target_canvas("LD2", file_ld2, variable_to_plot, bins, x_min, x_max, out_file)

    # Save to PDF files
    canvas_lh2.SaveAs("LH2_Yield_Comparison.pdf")
    canvas_ld2.SaveAs("LD2_Yield_Comparison.pdf")

    # Close output root file
    out_file.Close()

    print(f"Processing complete. Canvases saved as PDFs. All TObjects saved to {output_filename}")

if __name__ == "__main__":
    main()