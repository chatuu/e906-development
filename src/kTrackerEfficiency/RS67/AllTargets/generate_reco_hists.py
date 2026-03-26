import uproot
import numpy as np
import ROOT
import os
import argparse

# ==========================================
# ROOT Configuration
# ==========================================
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

# ==========================================
# Global Configuration
# ==========================================
input_npz_file  = "../GlobalEfficiencyCurve/interpolation_data_d1.npz"

# ==========================================
# Binning Definitions
# ==========================================
mass_bins_np = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
xf_bins_np = np.round(np.arange(0.0, 0.85, 0.05), 2)

# ==========================================
# Cut Function
# ==========================================
def e906_chuck_cuts(tree, cut=4.2, beam_offset=1.6):
    events = tree.arrays(library="np")
    class EventNamespace:
        def __init__(self, data):
            self.__dict__.update(data)
    e = EventNamespace(events)

    dimuon_cut = (
        (np.abs(e.dx) < 0.25) & (np.abs(e.dy - beam_offset) < 0.22) &
        (e.dz < -5.) & (e.dz > -280.) &
        (np.abs(e.dpx) < 1.8) & (np.abs(e.dpy) < 2.0) &
        (e.dpx**2 + e.dpy**2 < 5.) &
        (e.dpz < 116.) & (e.dpz > 38.) &
        (e.mass > cut) & (e.mass < 8.8) &
        (e.dx**2 + (e.dy - beam_offset)**2 < 0.06) &
        (e.xF < 0.95) & (e.xF > -0.1) &
        (e.xT > 0.05) & (e.xT <= 0.58) &
        (np.abs(e.costh) < 0.5) &
        (np.abs(e.trackSeparation) < 270.) &
        (e.chisq_dimuon < 18)
    )

    track1_cut = (
        (e.chisq1_target < 15.) & (e.pz1_st1 > 9.) & (e.pz1_st1 < 75.) &
        (e.nHits1 > 13) & (e.x1_t**2 + (e.y1_t - beam_offset)**2 < 320.) &
        (e.x1_d**2 + (e.y1_d - beam_offset)**2 < 1100.) &
        (e.x1_d**2 + (e.y1_d - beam_offset)**2 > 16.) &
        (e.chisq1_target < 1.5 * e.chisq1_upstream) &
        (e.chisq1_target < 1.5 * e.chisq1_dump) &
        (e.z1_v < -5.) & (e.z1_v > -320.) &
        (e.chisq1/(e.nHits1 - 5) < 12) & ((e.y1_st1)/(e.y1_st3 ) < 1.) & 
        (np.abs(np.abs(e.px1_st1 - e.px1_st3) - 0.416) < 0.008) & 
        (np.abs(e.py1_st1 - e.py1_st3) < 0.008) &
        (np.abs(e.pz1_st1 - e.pz1_st3) < 0.08) &
        ((e.y1_st1) * (e.y1_st3) > 0.) & (np.abs(e.py1_st1) > 0.02)
    )

    track2_cut = (
        (e.chisq2_target < 15.) & (e.pz2_st1 > 9.) & (e.pz2_st1 < 75.) &
        (e.nHits2 > 13) & (e.x2_t**2 + (e.y2_t - beam_offset)**2 < 320.) &
        (e.x2_d**2 + (e.y2_d - beam_offset)**2 < 1100.) &
        (e.x2_d**2 + (e.y2_d - beam_offset)**2 > 16.) &
        (e.chisq2_target < 1.5 * e.chisq2_upstream) &
        (e.chisq2_target < 1.5 * e.chisq2_dump) &
        (e.z2_v < -5.) & (e.z2_v > -320.) &
        (e.chisq2/(e.nHits2 - 5) < 12) & ((e.y2_st1 )/(e.y2_st3) < 1.) & 
        (np.abs(np.abs(e.px2_st1 - e.px2_st3) - 0.416) < 0.008) &
        (np.abs(e.py2_st1 - e.py2_st3) < 0.008) &
        (np.abs(e.pz2_st1 - e.pz2_st3) < 0.08) &
        ((e.y2_st1) * (e.y2_st3) > 0.) & (np.abs(e.py2_st1) > 0.02)
    )

    tracks_cut = (
        (np.abs(e.chisq1_target + e.chisq2_target - e.chisq_dimuon) < 2.) &
        ((e.y1_st3) * (e.y2_st3) < 0.) & 
        (e.nHits1 + e.nHits2 > 29) & (e.nHits1St1 + e.nHits2St1 > 8) &
        (np.abs(e.x1_st1 + e.x2_st1) < 42)
    )

    occ_cut = (
        (e.D1 < 400) & (e.D2 < 400) & (e.D3 < 400) &
        (e.D1 + e.D2 + e.D3 < 1000)
    )

    total_cut_mask = (track1_cut & track2_cut & tracks_cut & dimuon_cut & occ_cut)

    filtered_events = {}
    for key, val in events.items():
        filtered_events[key] = val[total_cut_mask]
        
    return filtered_events

# ==========================================
# Correlation & Plotting Logic
# ==========================================
def get_correlated_mean_and_error(eff_arr, err_arr, loc_arr):
    N = len(eff_arr)
    if N == 0: 
        return 0.0, 0.0
    
    mean_eff = np.mean(eff_arr)
    
    diff_matrix = np.abs(loc_arr[:, None] - loc_arr[None, :])
    correl_matrix = np.zeros((N, N))
    correl_matrix[diff_matrix == 0] = 1.0  
    correl_matrix[diff_matrix == 1] = 1.0  
    
    covar_matrix = correl_matrix * np.outer(err_arr, err_arr)
    total_variance = np.sum(covar_matrix)
    
    err_on_mean = np.sqrt(total_variance) / N
    return mean_eff, err_on_mean

def add_latex_to_bin(hist, x_center, y_center, value, error):
    val_f = float(value)
    err_f = float(error)
    
    # Check for NaN or Inf from empty bins and set them to 0.0
    if np.isnan(val_f) or np.isinf(val_f): val_f = 0.0
    if np.isnan(err_f) or np.isinf(err_f): err_f = 0.0
    
    latex_text = f"#splitline{{{val_f:.3f}}}{{#pm {err_f:.3f}}}"
    l = ROOT.TLatex(x_center, y_center, latex_text)
    l.SetTextSize(0.015)
    l.SetTextAlign(22)
    l.SetTextColor(ROOT.kBlack)
    hist.GetListOfFunctions().Add(l)

def process_and_plot_tree(data_dict, loc_data_array, hist_name, hist_title, out_pdf):
    print(f"  -> Processing {hist_name}...")
    h_reco = ROOT.TH2D(hist_name, f"{hist_title};Mass [GeV];x_{{F}}", 
                       len(mass_bins_np)-1, mass_bins_np, 
                       len(xf_bins_np)-1, xf_bins_np)
    h_reco.Sumw2()
    h_reco.SetStats(0)
    h_reco.GetXaxis().CenterTitle()
    h_reco.GetYaxis().CenterTitle()
    h_reco.GetXaxis().SetTitleOffset(1.2)
    h_reco.GetYaxis().SetTitleOffset(1.2)

    for i_m in range(len(mass_bins_np) - 1):
        m_low, m_high = mass_bins_np[i_m], mass_bins_np[i_m+1]
        m_center = (m_low + m_high) / 2.0
        root_x = i_m + 1
        
        for i_x in range(len(xf_bins_np) - 1):
            x_low, x_high = xf_bins_np[i_x], xf_bins_np[i_x+1]
            x_center = (x_low + x_high) / 2.0
            root_y = i_x + 1
            
            mask = (data_dict["mass"] >= m_low) & (data_dict["mass"] < m_high) & \
                   (data_dict["xF"] >= x_low) & (data_dict["xF"] < x_high)
            
            eff_arr = data_dict["recoeff"][mask]
            err_arr = data_dict["recoeff_error"][mask]
            loc_arr = loc_data_array[mask]
            
            mean_eff, corr_err = get_correlated_mean_and_error(eff_arr, err_arr, loc_arr)
            
            h_reco.SetBinContent(root_x, root_y, mean_eff)
            h_reco.SetBinError(root_x, root_y, corr_err)
            
            add_latex_to_bin(h_reco, m_center, x_center, mean_eff, corr_err)

    print(f"  -> Drawing and saving to {out_pdf}...")
    c = ROOT.TCanvas(f"c_{hist_name}", "", 1200, 900)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTickx(1)
    c.SetTicky(1)
    
    h_reco.Draw("COLZ")
    c.SaveAs(out_pdf)
    c.Close()

def generate_recoeff_histograms(input_root_file, target_str):
    print(f"\n==========================================")
    print(f"Processing Target: {target_str}")
    print(f"Input File: {input_root_file}")
    print(f"==========================================")

    print(f"Loading data from {input_root_file}...")
    try:
        f_in = uproot.open(input_root_file)
        
        print("  - Applying cuts to 'result' tree...")
        data_tot = e906_chuck_cuts(f_in["result"])
        
        print("  - Applying cuts to 'result_mix' tree...")
        data_mix = e906_chuck_cuts(f_in["result_mix"])
        
    except Exception as e:
        print(f"Error reading root file: {e}")
        return

    print("\nLoading NPZ file to recreate interpolation bins...")
    npz_data = np.load(input_npz_file)
    x_curve = npz_data['x']
    
    def get_d1_array(events_dict):
        if 'ptrk_D1' in events_dict and 'ntrk_D1' in events_dict:
            return 0.5 * (events_dict['ptrk_D1'] + events_dict['ntrk_D1'])
        return events_dict['D1']

    loc_data_tot = np.digitize(get_d1_array(data_tot), x_curve) - 1
    loc_data_mix = np.digitize(get_d1_array(data_mix), x_curve) - 1

    print(f"\nCalculating correlated errors per bin for {target_str} Total Yield...")
    process_and_plot_tree(
        data_dict=data_tot, 
        loc_data_array=loc_data_tot, 
        hist_name=f"E_total_reco_{target_str}", 
        hist_title=f"Avg Reco Eff with Correlated Errors (Total) ({target_str})", 
        out_pdf=f"E_total_reco_{target_str}_correlated.pdf"
    )

    print(f"\nCalculating correlated errors per bin for {target_str} Mix Yield...")
    process_and_plot_tree(
        data_dict=data_mix, 
        loc_data_array=loc_data_mix, 
        hist_name=f"E_mix_reco_{target_str}", 
        hist_title=f"Avg Reco Eff with Correlated Errors (Mix) ({target_str})", 
        out_pdf=f"E_mix_reco_{target_str}_correlated.pdf"
    )
    
    print("\nDone! Both Total and Mix PDFs generated.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate 2D Reco Efficiency Histograms with Correlated Errors")
    parser.add_argument("-i", "--input", required=True, help="Path to the input ROOT file")
    parser.add_argument("-t", "--target", required=True, choices=["LH2", "LD2", "Flask"], help="Target string (LH2, LD2, Flask)")
    
    args = parser.parse_args()
    
    generate_recoeff_histograms(args.input, args.target)