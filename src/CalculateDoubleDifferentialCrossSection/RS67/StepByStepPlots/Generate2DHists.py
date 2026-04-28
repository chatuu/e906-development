import uproot
import numpy as np
import ROOT
import sys
import os
import math

# Set ROOT batch mode to avoid opening windows
ROOT.gROOT.SetBatch(True)
# Remove the statistics box
ROOT.gStyle.SetOptStat(0)
# Set color palette
ROOT.gStyle.SetPalette(ROOT.kBird)

# ==========================================
# 1. Define Physics Constants
# ==========================================
PROTONS_ON_TARGET_LH2 = 1.57319e+17
PROTONS_ON_TARGET_FLASK = 3.57904e+16
LH2_TARGET_DENSITY_MOL_CM2 = 3.597
AVOGADRO_CONSTANT = 6.022e23
PROTONS_PER_NUCLEON_LH2 = 1.008
XF_BIN_WIDTH = 0.05

# Beam Attenuation
val_exp = -(50.8 * 0.0708) / 52.0
BEAM_ATTENUATION = (52.0 / (0.0708 * 50.8)) * (1.0 - math.exp(val_exp))

# Global Normalization Constant
GLOBAL_CONSTANT = (PROTONS_PER_NUCLEON_LH2 * 1e33) / (
    LH2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LH2 * BEAM_ATTENUATION * XF_BIN_WIDTH
)

# Flask Normalization
FLASK_NORM = PROTONS_ON_TARGET_LH2 / PROTONS_ON_TARGET_FLASK

print("=== Physics Constants ===")
print(f"PoT LH2: {PROTONS_ON_TARGET_LH2:.4e}")
print(f"Flask Norm: {FLASK_NORM:.4f}")
print(f"Beam Attenuation: {BEAM_ATTENUATION:.4f}")
print(f"Global Constant: {GLOBAL_CONSTANT:.4e}")
print("=========================")


# ==========================================
# 2. Define the User's Cut Function
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
# 3. Setup Binning and Histograms
# ==========================================

mass_bins_np = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
xf_bins_np = np.round(np.arange(0.0, 0.85, 0.05), 2)

def create_histograms(file_label):
    hists = {}
    
    def make_th2(name, title):
        h = ROOT.TH2D(f"{name}_{file_label}", f"{title} ({file_label});Mass [GeV];x_{{F}}", 
                      len(mass_bins_np)-1, mass_bins_np, 
                      len(xf_bins_np)-1, xf_bins_np)
        h.Sumw2()
        h.SetStats(0)
        h.GetXaxis().CenterTitle()
        h.GetYaxis().CenterTitle()
        h.GetXaxis().SetTitleOffset(1.2)
        h.GetYaxis().SetTitleOffset(1.2)
        return h

    hists["Y_total"] = make_th2("Y_total", "Total Yield (result)")
    hists["Y_mix"] = make_th2("Y_mix", "Mix Yield (result_mix)")
    hists["E_total_reco"] = make_th2("E_total_reco", "Avg Reco Eff (Total)")
    hists["E_mix_reco"] = make_th2("E_mix_reco", "Avg Reco Eff (Mix)")
    hists["E_total_hodo"] = make_th2("E_total_hodo", "Avg Hodo Eff (Total)")
    hists["E_mix_hodo"] = make_th2("E_mix_hodo", "Avg Hodo Eff (Mix)")
    hists["E_total_final"] = make_th2("E_total_final", "Avg Final Eff (Total)")
    hists["E_mix_final"] = make_th2("E_mix_final", "Avg Final Eff (Mix)")
    hists["E_final_signal"] = make_th2("E_final_signal", "Avg Signal Efficiency")
    
    hists["Y_corrected"] = make_th2("Y_corrected", "Corrected Yield (Total Error)")
    hists["Y_corrected_stat"] = make_th2("Y_corrected_stat", "Corrected Yield (Stat Error)")
    hists["Y_corrected_sys"] = make_th2("Y_corrected_sys", "Corrected Yield (Sys Error)")

    return hists

# ==========================================
# 4. Calculation Helpers
# ==========================================

def get_weighted_mean_and_error(eff_arr, err_arr):
    if len(eff_arr) == 0: return 0.0, 0.0
    mean = np.mean(eff_arr)
    err_on_mean = np.sqrt(np.sum(err_arr**2)) / len(eff_arr)
    return mean, err_on_mean

def add_latex_to_bin(hist, x_center, y_center, value, error):
    val_f = float(value)
    err_f = float(error)
    if np.isnan(val_f) or np.isinf(val_f): val_f = 0.0
    if np.isnan(err_f) or np.isinf(err_f): err_f = 0.0
    
    latex_text = f"#splitline{{{val_f:.3f}}}{{#pm {err_f:.3f}}}"
    l = ROOT.TLatex(x_center, y_center, latex_text)
    l.SetTextSize(0.015)
    l.SetTextAlign(22)
    l.SetTextColor(ROOT.kBlack)
    hist.GetListOfFunctions().Add(l)

# ==========================================
# 5. Table Generation (Combined)
# ==========================================

def generate_combined_table_latex(hists, file_label):
    filename = f"Table_Combined_{file_label}.tex"
    print(f"  Generating combined table: {filename}")
    if "Y_total" not in hists: return

    h_tot = hists["Y_total"]
    h_mix = hists["Y_mix"]
    h_eff = hists["E_final_signal"]
    h_corr = hists["Y_corrected"] 

    content = r"""\begin{longtable}{|c|c|c|c|c|c|}
\caption{Combined Data Table for %s} \label{tab:combined_%s} \\
\hline
\textbf{xF bin} & \textbf{Mass (GeV)} & \textbf{$Y_{total} \pm \delta Y_{total}$} & \textbf{$Y_{mix} \pm \delta Y_{mix}$} & \textbf{$\langle E_{signal} \rangle$} & \textbf{$Y_{corr} \pm \delta Y_{corr}$} \\
\hline
\endfirsthead

\multicolumn{6}{c}%%
{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
\hline
\textbf{xF bin} & \textbf{Mass (GeV)} & \textbf{$Y_{total} \pm \delta Y_{total}$} & \textbf{$Y_{mix} \pm \delta Y_{mix}$} & \textbf{$\langle E_{signal} \rangle$} & \textbf{$Y_{corr} \pm \delta Y_{corr}$} \\
\hline
\endhead

\hline \multicolumn{6}{|r|}{{Continued on next page}} \\ \hline
\endfoot

\hline
\endlastfoot
""" % (file_label, file_label)

    for iy in range(1, h_tot.GetNbinsY() + 1):
        xf_min = h_tot.GetYaxis().GetBinLowEdge(iy)
        xf_max = xf_min + h_tot.GetYaxis().GetBinWidth(iy)
        for ix in range(1, h_tot.GetNbinsX() + 1):
            mass_min = h_tot.GetXaxis().GetBinLowEdge(ix)
            mass_max = mass_min + h_tot.GetXaxis().GetBinWidth(ix)
            
            val_tot = h_tot.GetBinContent(ix, iy)
            err_tot = h_tot.GetBinError(ix, iy)
            val_mix = h_mix.GetBinContent(ix, iy)
            err_mix = h_mix.GetBinError(ix, iy)
            val_eff = h_eff.GetBinContent(ix, iy)
            err_eff = h_eff.GetBinError(ix, iy)
            val_corr = h_corr.GetBinContent(ix, iy)
            err_corr = h_corr.GetBinError(ix, iy)
            
            if val_tot == 0 and err_tot == 0: continue
            
            s_xf = f"[{xf_min:.2f}, {xf_max:.2f})"
            s_mass = f"[{mass_min:.2f}, {mass_max:.2f})"
            
            s_tot = f"{int(round(val_tot))} $\\pm$ {err_tot:.3f}"
            s_mix = f"{int(round(val_mix))} $\\pm$ {err_mix:.3f}"
            s_eff = f"{val_eff:.4f} $\\pm$ {err_eff:.4f}"
            s_corr = f"{val_corr:.3f} $\\pm$ {err_corr:.3f}"
            
            line = f"{s_xf} & {s_mass} & {s_tot} & {s_mix} & {s_eff} & {s_corr} \\\\ \n"
            content += line
            content += r"\hline" + "\n"

    content += r"\end{longtable}" + "\n"
    with open(filename, "w") as f: f.write(content)

# ==========================================
# 6. Basic Processing
# ==========================================

def process_file_data(filename, file_label, output_file):
    print(f"\nProcessing {filename} as {file_label}...")
    try: f = uproot.open(filename)
    except Exception as e:
        print(f"Error opening {filename}: {e}")
        return {}

    output_file.cd()
    hists = create_histograms(file_label)
    
    try:
        data_total = e906_chuck_cuts(f["result"])
        data_mix = e906_chuck_cuts(f["result_mix"])
    except KeyError as e:
        print(f"  Missing TTree in file: {e}")
        return {}

    for i_m in range(len(mass_bins_np) - 1):
        m_low = mass_bins_np[i_m]
        m_high = mass_bins_np[i_m+1]
        m_center = (m_low + m_high) / 2.0
        
        for i_x in range(len(xf_bins_np) - 1):
            x_low = xf_bins_np[i_x]
            x_high = xf_bins_np[i_x+1]
            x_center = (x_low + x_high) / 2.0
            
            mask_tot = (data_total["mass"] >= m_low) & (data_total["mass"] < m_high) & \
                       (data_total["xF"] >= x_low) & (data_total["xF"] < x_high)
            mask_mix = (data_mix["mass"] >= m_low) & (data_mix["mass"] < m_high) & \
                       (data_mix["xF"] >= x_low) & (data_mix["xF"] < x_high)
            
            N_tot = np.sum(mask_tot)
            N_mix = np.sum(mask_mix)
            err_N_tot = np.sqrt(N_tot)
            err_N_mix = np.sqrt(N_mix)

            eff_reco_tot = data_total["recoeff"][mask_tot]
            err_reco_tot = data_total["recoeff_error"][mask_tot]
            eff_hodo_tot = data_total["hodoeff"][mask_tot]
            err_hodo_tot = data_total["hodoeff_error"][mask_tot]
            eff_final_tot = eff_reco_tot * eff_hodo_tot
            err_final_tot = np.sqrt( (eff_hodo_tot * err_reco_tot)**2 + (eff_reco_tot * err_hodo_tot)**2 )

            eff_reco_mix = data_mix["recoeff"][mask_mix]
            err_reco_mix = data_mix["recoeff_error"][mask_mix]
            eff_hodo_mix = data_mix["hodoeff"][mask_mix]
            err_hodo_mix = data_mix["hodoeff_error"][mask_mix]
            eff_final_mix = eff_reco_mix * eff_hodo_mix
            err_final_mix = np.sqrt( (eff_hodo_mix * err_reco_mix)**2 + (eff_reco_mix * err_hodo_mix)**2 )

            mean_reco_tot, e_mean_reco_tot = get_weighted_mean_and_error(eff_reco_tot, err_reco_tot)
            mean_hodo_tot, e_mean_hodo_tot = get_weighted_mean_and_error(eff_hodo_tot, err_hodo_tot)
            mean_final_tot, e_mean_final_tot = get_weighted_mean_and_error(eff_final_tot, err_final_tot)
            
            mean_reco_mix, e_mean_reco_mix = get_weighted_mean_and_error(eff_reco_mix, err_reco_mix)
            mean_hodo_mix, e_mean_hodo_mix = get_weighted_mean_and_error(eff_hodo_mix, err_hodo_mix)
            mean_final_mix, e_mean_final_mix = get_weighted_mean_and_error(eff_final_mix, err_final_mix)

            diff_yield = N_tot - N_mix
            val_sig_eff = 0.0
            err_sig_eff = 0.0
            
            if diff_yield > 0:
                numerator = (N_tot * mean_final_tot) - (N_mix * mean_final_mix)
                val_sig_eff = numerator / diff_yield
                term1 = (N_tot * e_mean_final_tot)**2
                term2 = (N_mix * e_mean_final_mix)**2
                err_sig_eff = (1.0 / diff_yield) * np.sqrt(term1 + term2)
            
            val_corr_yield = 0.0
            err_corr_stat = 0.0
            err_corr_sys = 0.0
            err_corr_total = 0.0

            if val_sig_eff > 0 and diff_yield > 0:
                val_corr_yield = diff_yield / val_sig_eff
                err_corr_stat = np.sqrt(N_tot + N_mix) / val_sig_eff
                err_corr_sys = val_corr_yield * (err_sig_eff / val_sig_eff)
                err_corr_total = np.sqrt(err_corr_stat**2 + err_corr_sys**2)

            def fill_bin(key, val, err):
                hists[key].SetBinContent(i_m + 1, i_x + 1, val)
                hists[key].SetBinError(i_m + 1, i_x + 1, err)
                if "stat" not in key and "sys" not in key:
                    add_latex_to_bin(hists[key], m_center, x_center, val, err)

            fill_bin("Y_total", N_tot, err_N_tot)
            fill_bin("Y_mix", N_mix, err_N_mix)
            fill_bin("E_total_reco", mean_reco_tot, e_mean_reco_tot)
            fill_bin("E_mix_reco", mean_reco_mix, e_mean_reco_mix)
            fill_bin("E_total_hodo", mean_hodo_tot, e_mean_hodo_tot)
            fill_bin("E_mix_hodo", mean_hodo_mix, e_mean_hodo_mix)
            fill_bin("E_total_final", mean_final_tot, e_mean_final_tot)
            fill_bin("E_mix_final", mean_final_mix, e_mean_final_mix)
            fill_bin("E_final_signal", val_sig_eff, err_sig_eff)
            
            fill_bin("Y_corrected", val_corr_yield, err_corr_total)
            fill_bin("Y_corrected_stat", val_corr_yield, err_corr_stat)
            fill_bin("Y_corrected_sys", val_corr_yield, err_corr_sys)

    print(f"  Saving PDFs for {file_label}...")
    for name, hist in hists.items():
        if "stat" in name or "sys" in name: 
            continue 
        
        c_name = f"c_{name}_{file_label}"
        c = ROOT.TCanvas(c_name, c_name, 1200, 900)
        c.SetRightMargin(0.15)
        c.SetLeftMargin(0.12)
        c.SetBottomMargin(0.12)
        c.SetTickx(1)
        c.SetTicky(1)
        
        local_min = sys.float_info.max
        local_max = -sys.float_info.max
        has_data = False
        
        for ix in range(1, hist.GetNbinsX() + 1):
            for iy in range(1, hist.GetNbinsY() + 1):
                val = hist.GetBinContent(ix, iy)
                if val != 0.0:
                    has_data = True
                    if val < local_min: local_min = val
                    if val > local_max: local_max = val
        if has_data:
            hist.SetMinimum(local_min)
            hist.SetMaximum(local_max)
        else:
            hist.SetMinimum(0.0)
            hist.SetMaximum(1.0)
        
        hist.Draw("COLZ")
        c.SaveAs(f"{name}_{file_label}.pdf")
        c.Close()

    generate_combined_table_latex(hists, file_label)
    return hists

# ==========================================
# 7. Subtracted Yield Calculation
# ==========================================

def generate_subtracted_plot(hists_lh2, hists_flask, output_file):
    print("\nGenerating Subtracted Plot: Y_corrected_LH2 - Y_corrected_Flask...")
    
    if "Y_corrected_stat" not in hists_lh2 or "Y_corrected_stat" not in hists_flask:
        print("Error: Y_corrected_stat histogram missing.")
        return None

    output_file.cd()
    
    h_lh2_stat = hists_lh2["Y_corrected_stat"]
    h_lh2_sys = hists_lh2["Y_corrected_sys"]
    h_flask_stat = hists_flask["Y_corrected_stat"]
    h_flask_sys = hists_flask["Y_corrected_sys"]
    
    name = "Y_corrected_Subtracted"
    title = "Corrected Yield (LH2 - Flask)"
    h_sub = ROOT.TH2D(name, f"{title};Mass [GeV];x_{{F}}", 
                      len(mass_bins_np)-1, mass_bins_np, 
                      len(xf_bins_np)-1, xf_bins_np)
    h_sub.Sumw2()
    h_sub.SetStats(0)
    h_sub.GetXaxis().CenterTitle()
    h_sub.GetYaxis().CenterTitle()
    h_sub.GetXaxis().SetTitleOffset(1.2)
    h_sub.GetYaxis().SetTitleOffset(1.2)

    h_sub_stat = h_sub.Clone("Y_corrected_Subtracted_stat")
    h_sub_sys = h_sub.Clone("Y_corrected_Subtracted_sys")

    for i_m in range(len(mass_bins_np) - 1):
        m_center = (mass_bins_np[i_m] + mass_bins_np[i_m+1]) / 2.0
        bin_x = i_m + 1
        
        for i_x in range(len(xf_bins_np) - 1):
            x_center = (xf_bins_np[i_x] + xf_bins_np[i_x+1]) / 2.0
            bin_y = i_x + 1
            
            y_lh2 = h_lh2_stat.GetBinContent(bin_x, bin_y)
            e_lh2_stat = h_lh2_stat.GetBinError(bin_x, bin_y)
            e_lh2_sys = h_lh2_sys.GetBinError(bin_x, bin_y)
            
            y_flask = h_flask_stat.GetBinContent(bin_x, bin_y)
            e_flask_stat = h_flask_stat.GetBinError(bin_x, bin_y)
            e_flask_sys = h_flask_sys.GetBinError(bin_x, bin_y)
            
            val_sub = y_lh2 - (FLASK_NORM * y_flask)
            
            err_sub_stat = np.sqrt(e_lh2_stat**2 + (FLASK_NORM * e_flask_stat)**2)
            err_sub_sys = np.sqrt(e_lh2_sys**2 + (FLASK_NORM * e_flask_sys)**2)
            err_sub_total = np.sqrt(err_sub_stat**2 + err_sub_sys**2)
            
            h_sub.SetBinContent(bin_x, bin_y, val_sub)
            h_sub.SetBinError(bin_x, bin_y, err_sub_total)
            add_latex_to_bin(h_sub, m_center, x_center, val_sub, err_sub_total)

            h_sub_stat.SetBinContent(bin_x, bin_y, val_sub)
            h_sub_stat.SetBinError(bin_x, bin_y, err_sub_stat)

            h_sub_sys.SetBinContent(bin_x, bin_y, val_sub)
            h_sub_sys.SetBinError(bin_x, bin_y, err_sub_sys)

    c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTickx(1)
    c.SetTicky(1)
    
    local_min = sys.float_info.max
    local_max = -sys.float_info.max
    has_data = False
    for ix in range(1, h_sub.GetNbinsX() + 1):
        for iy in range(1, h_sub.GetNbinsY() + 1):
            val = h_sub.GetBinContent(ix, iy)
            if val != 0.0:
                has_data = True
                if val < local_min: local_min = val
                if val > local_max: local_max = val
    if has_data:
        h_sub.SetMinimum(local_min)
        h_sub.SetMaximum(local_max)
    else:
        h_sub.SetMinimum(0.0)
        h_sub.SetMaximum(1.0)
        
    h_sub.Draw("COLZ")
    c.SaveAs(f"{name}.pdf")
    c.Close()
    print("  Subtracted plot saved.")
    
    return {"stat": h_sub_stat, "sys": h_sub_sys}

# ==========================================
# 8. Cross Section Calculation & Plotting
# ==========================================

def calculate_and_plot_cross_section(h_sub_dict, output_file):
    print("\nCalculating and Plotting Double Differential Cross-Sections...")
    acc_path = "acceptance_mass_xF.root"
    theory_ct18_path = "CT18_xFnew_p_1sigma.root"
    theory_nnpdf_path = "NNPDF40_xFnew_p.root"
    
    if not os.path.exists(acc_path):
        raise FileNotFoundError(f"CRITICAL ERROR: Acceptance file '{acc_path}' not found!")

    # Open Files
    try:
        acc_file = ROOT.TFile.Open(acc_path)
        if not acc_file or acc_file.IsZombie():
            raise OSError(f"Acceptance file corrupted.")
            
        f_ct18 = None
        if os.path.exists(theory_ct18_path):
            f_ct18 = ROOT.TFile.Open(theory_ct18_path)
            
        f_nnpdf = None
        if os.path.exists(theory_nnpdf_path):
            f_nnpdf = ROOT.TFile.Open(theory_nnpdf_path)
            
    except Exception as e:
        print(f"Error opening files: {e}")
        sys.exit(1)

    output_file.cd()
    
    h_sub_stat = h_sub_dict["stat"]
    h_sub_sys = h_sub_dict["sys"]

    n_xf_bins = len(xf_bins_np) - 1
    n_mass_bins = len(mass_bins_np) - 1
    
    for i_x in range(n_xf_bins):
        xf_bin_index = i_x 
        root_bin_y = i_x + 1
        
        xf_min = xf_bins_np[i_x]
        xf_max = xf_bins_np[i_x+1]
        
        acc_hist_name = f"h_ratio_LH2_xF_bin{xf_bin_index}"
        h_acc = acc_file.Get(acc_hist_name)
        
        if not h_acc:
            print(f"Warning: Acceptance histogram '{acc_hist_name}' not found.")
            continue
        
        print(f"\n--- xF Bin {xf_bin_index}: {xf_min:.2f} < xF < {xf_max:.2f} ---")

        # --- Data Graph ---
        g_xsec = ROOT.TGraphErrors()
        g_xsec.SetName(f"g_xsec_xF_{xf_bin_index}")
        g_xsec.SetTitle(f"Cross Section for {xf_min:.2f} < x_{{F}} < {xf_max:.2f};Mass [GeV];M^{{3}} d^{{2}}\\sigma / dM dx_{{F}} [nb GeV^{{2}}]")
        
        # --- Systematic Error Band Graph ---
        g_sys = ROOT.TGraphErrors()
        g_sys.SetName(f"g_sys_xF_{xf_bin_index}")
        
        point_idx = 0
        
        # MODIFIED: Variables to track dynamic min/max for the Y-axis
        y_min_data = sys.float_info.max
        y_max_data = -sys.float_info.max
        
        for i_m in range(n_mass_bins):
            root_bin_x = i_m + 1
            
            mass_min = mass_bins_np[i_m]
            mass_max = mass_bins_np[i_m+1]
            mass_center = (mass_min + mass_max) / 2.0
            mass_width = mass_max - mass_min
            
            Y_sub = h_sub_stat.GetBinContent(root_bin_x, root_bin_y)
            Y_sub_stat_err = h_sub_stat.GetBinError(root_bin_x, root_bin_y)
            Y_sub_sys_err = h_sub_sys.GetBinError(root_bin_x, root_bin_y)
            
            if Y_sub <= 0: continue
                
            acc_bin = h_acc.FindBin(mass_center)
            acceptance = h_acc.GetBinContent(acc_bin)
            acceptance_err = h_acc.GetBinError(acc_bin)
            
            print(f"  Mass Bin {mass_center:.2f} GeV: Acceptance = {acceptance:.5f} +/- {acceptance_err:.5f}")
            
            if acceptance <= 0: continue
            
            numerator = GLOBAL_CONSTANT * Y_sub
            denominator = mass_width * acceptance
            xsec = numerator / denominator
            
            # --- Statistical Error Calculation (Bars) ---
            stat_unc = (Y_sub_stat_err / Y_sub) * xsec
            
            # --- Systematic Error Calculation (Band) ---
            sys_acc = (acceptance_err / acceptance) * xsec
            sys_tot_yield = (Y_sub_sys_err / Y_sub) * xsec
            sys_psip_cont = 0.0 # Placeholder
            total_systematic_unc = np.sqrt(sys_acc**2 + sys_tot_yield**2 + sys_psip_cont**2)
            
            # --- Apply M^3 Scaling ---
            scaling_factor = mass_center**3
            scaled_xsec = xsec * scaling_factor
            scaled_stat_err = stat_unc * scaling_factor
            scaled_sys_err = total_systematic_unc * scaling_factor
            
            # Fill Data Graph (Statistical Error Bars)
            g_xsec.SetPoint(point_idx, mass_center, scaled_xsec)
            g_xsec.SetPointError(point_idx, 0, scaled_stat_err)
            
            # Fill Systematic Graph (Error Band)
            g_sys.SetPoint(point_idx, mass_center, scaled_xsec)
            g_sys.SetPointError(point_idx, mass_width/2.0, scaled_sys_err)
            
            # MODIFIED: Track min/max for dynamic Y-axis range
            max_err_for_range = max(scaled_stat_err, scaled_sys_err)
            y_high = scaled_xsec + max_err_for_range
            y_low  = scaled_xsec - max_err_for_range
            
            if y_low <= 0:
                y_low = scaled_xsec * 0.5 # Log scale safeguard to prevent <= 0 bounds
                
            if y_high > y_max_data:
                y_max_data = y_high
            if y_low < y_min_data:
                y_min_data = y_low
            
            point_idx += 1
            
        # --- Plotting ---
        if g_xsec.GetN() > 0:
            c_xsec = ROOT.TCanvas(f"c_xsec_{xf_bin_index}", f"Cross Section xF {xf_bin_index}", 800, 600)
            c_xsec.SetLogy()
            c_xsec.SetTickx(1)
            c_xsec.SetTicky(1)
            
            mg = ROOT.TMultiGraph()
            mg.SetTitle(f"Cross Section for {xf_min:.2f} < x_{{F}} < {xf_max:.2f};Mass [GeV];M^{{3}} d^{{2}}\\sigma / dM dx_{{F}} [nb GeV^{{2}}]")
            
            # MODIFIED: Dynamically Set Y-axis Range based on tracked data
            if y_min_data < y_max_data:
                mg.SetMinimum(y_min_data * 0.2)  # Space at bottom
                mg.SetMaximum(y_max_data * 5.0)  # Space at top for legend
            else:
                # Fallback just in case
                mg.SetMinimum(1e-3)
                mg.SetMaximum(3.0)
            
            # --- Prepare Legend ---
            leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
            leg.SetBorderSize(0)
            
            gr_name = f"gr_xFbin{xf_bin_index}"
            
            # --- 1. Add Theoretical Curves FIRST (Background) ---
            if f_ct18:
                g_ct18 = f_ct18.Get(gr_name)
                if g_ct18:
                    g_ct18.SetLineColor(ROOT.kGreen + 2)
                    g_ct18.SetFillColorAlpha(ROOT.kGreen - 5, 0.5)
                    g_ct18.SetFillStyle(3002)
                    mg.Add(g_ct18, "L3") 
                    leg.AddEntry(g_ct18, "CT18 NLO", "lf") 
            
            if f_nnpdf:
                g_nnpdf = f_nnpdf.Get(gr_name)
                if g_nnpdf:
                    g_nnpdf.SetLineColor(ROOT.kBlue + 2)
                    g_nnpdf.SetFillColorAlpha(ROOT.kAzure + 1, 0.5)
                    g_nnpdf.SetFillStyle(3002)
                    mg.Add(g_nnpdf, "L3") 
                    leg.AddEntry(g_nnpdf, "NNPDF4.0 NLO", "lf") 
            
            # --- 2. Add Systematic Error Band (Middle) ---
            g_sys.SetMarkerSize(0)
            g_sys.SetLineColor(ROOT.kRed)
            g_sys.SetFillColorAlpha(ROOT.kPink - 9, 0.5)
            g_sys.SetFillStyle(1001)
            mg.Add(g_sys, "2") # "2" draws error rectangles
            leg.AddEntry(g_sys, "Systematic Unc.", "f")

            # --- 3. Add Data Points LAST (Foreground) ---
            g_xsec.SetMarkerStyle(20)
            g_xsec.SetMarkerColor(ROOT.kRed)
            g_xsec.SetLineColor(ROOT.kRed)
            mg.Add(g_xsec, "P") 
            leg.AddEntry(g_xsec, "Data (LH2)", "lep")
            
            mg.Draw("A") 
            mg.GetXaxis().CenterTitle()
            mg.GetYaxis().CenterTitle()
            mg.GetXaxis().SetLimits(4.2, 8.7) # Explicit x-axis range
            
            leg.Draw()
            
            plot_name = f"CrossSection_xF_{xf_min:.2f}_{xf_max:.2f}.pdf"
            c_xsec.SaveAs(plot_name)
            c_xsec.Close()
            print(f"  Saved plot: {plot_name}")
        else:
            print(f"  No valid points for xF bin {xf_bin_index}.")

    acc_file.Close()
    if f_ct18: f_ct18.Close()
    if f_nnpdf: f_nnpdf.Close()


# ==========================================
# Main Execution
# ==========================================
def main():
    output_filename = "Processed_Kinematic_Hists.root"
    print(f"Creating ROOT output file: {output_filename}")
    out_file = ROOT.TFile(output_filename, "RECREATE")

    hists_lh2 = process_file_data("../../../HodoEfficiency/RS67/LH2/merged_RS67_3089_LH2_recoeff_hodoeff.root", "LH2", out_file)
    hists_flask = process_file_data("../../../HodoEfficiency/RS67/EmptyFlask/merged_RS67_3089_Flask_recoeff_hodoeff.root", "Flask", out_file)
    hists_ld2 = process_file_data("../../../HodoEfficiency/RS67/LD2/merged_RS67_3089_LD2_recoeff_hodoeff.root", "LD2", out_file)

    h_sub_dict = None
    if hists_lh2 and hists_flask:
        h_sub_dict = generate_subtracted_plot(hists_lh2, hists_flask, out_file)
    else:
        print("\nSkipping subtraction plot because LH2 or Flask histograms were not generated successfully.")

    if h_sub_dict:
        calculate_and_plot_cross_section(h_sub_dict, out_file)

    out_file.Close()
    print("\nAll histograms, tables, and cross-section plots generated successfully.")

if __name__ == "__main__":
    main()