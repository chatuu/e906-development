import uproot
import numpy as np
import ROOT
import sys
import os
import math
import csv

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
PROTONS_ON_TARGET_LD2 = 7.51113e+16
PROTONS_ON_TARGET_FLASK = 3.57904e+16

LH2_TARGET_DENSITY_MOL_CM2 = 3.5966
LD2_TARGET_DENSITY_MOL_CM2 = 8.0431

LH2_TARGET_LENGTH_CM = 50.8
LD2_TARGET_LENGTH_CM = 50.8

LH2_TARGET_DENSITY_MOL_CM3 = 0.0708
LD2_TARGET_DENSITY_MOL_CM3 = 0.163

AVOGADRO_CONSTANT = 6.022e23
NUCLEONS_PER_NUCLEUS_LH2 = 1.008
NUCLEONS_PER_NUCLEUS_LD2 = 2.014

NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2 = 52.0
NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2 = 54.7

XF_BIN_WIDTH = 0.05
THD_THH_RATIO = 0.1084/3.5966  # Ratio of target thicknesses T_HD / T_HH

# Beam Attenuation
val_exp_LH2 = -(LH2_TARGET_LENGTH_CM * LH2_TARGET_DENSITY_MOL_CM3) / NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2
BEAM_ATTENUATION_LH2 = (NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2 / (LH2_TARGET_DENSITY_MOL_CM3 * LH2_TARGET_LENGTH_CM)) * (1.0 - math.exp(val_exp_LH2))

val_exp_LD2 = -(LD2_TARGET_LENGTH_CM * LD2_TARGET_DENSITY_MOL_CM3) / NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2
BEAM_ATTENUATION_LD2 = (NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2 / (LD2_TARGET_DENSITY_MOL_CM3 * LD2_TARGET_LENGTH_CM)) * (1.0 - math.exp(val_exp_LD2))

# Global Normalization Constants
GLOBAL_CONSTANT_LH2 = (NUCLEONS_PER_NUCLEUS_LH2 * 1e33) / (
    LH2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LH2 * BEAM_ATTENUATION_LH2 * XF_BIN_WIDTH
)

GLOBAL_CONSTANT_LD2 = (NUCLEONS_PER_NUCLEUS_LD2 * 1e33) / (
    LD2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LD2 * BEAM_ATTENUATION_LD2 * XF_BIN_WIDTH
)

# Flask Normalization Factors
FLASK_NORM_LH2 = PROTONS_ON_TARGET_LH2 / PROTONS_ON_TARGET_FLASK
FLASK_NORM_LD2 = PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_FLASK
LH2_TO_LD2_NORM = THD_THH_RATIO * (PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_LH2)

print("=== Physics Constants ===")
print(f"PoT LH2:   {PROTONS_ON_TARGET_LH2:.4e}")
print(f"PoT LD2:   {PROTONS_ON_TARGET_LD2:.4e}")
print(f"PoT Flask: {PROTONS_ON_TARGET_FLASK:.4e}")
print(f"Flask Norm LH2: {FLASK_NORM_LH2:.4f}")
print(f"Flask Norm LD2: {FLASK_NORM_LD2:.4f}")
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
    
    hists["Mass_Centroid"] = make_th2("Mass_Centroid", "Data-Driven Mass Centroid")

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

def extract_target_stats(data_tot, data_mix, m_low, m_high, x_low, x_high):
    """Extract yields, efficiency, and array values for a given bin."""
    mask_tot = (data_tot["mass"] >= m_low) & (data_tot["mass"] < m_high) & \
               (data_tot["xF"] >= x_low) & (data_tot["xF"] < x_high)
    mask_mix = (data_mix["mass"] >= m_low) & (data_mix["mass"] < m_high) & \
               (data_mix["xF"] >= x_low) & (data_mix["xF"] < x_high)

    N_tot = np.sum(mask_tot)
    N_mix = np.sum(mask_mix)
    sum_mass_tot = np.sum(data_tot["mass"][mask_tot]) if N_tot > 0 else 0.0
    sum_mass_mix = np.sum(data_mix["mass"][mask_mix]) if N_mix > 0 else 0.0

    eff_reco_tot = data_tot["recoeff"][mask_tot]
    err_reco_tot = data_tot["recoeff_error"][mask_tot]
    eff_hodo_tot = data_tot["hodoeff"][mask_tot]
    err_hodo_tot = data_tot["hodoeff_error"][mask_tot]
    eff_final_tot = eff_reco_tot * eff_hodo_tot
    err_final_tot = np.sqrt( (eff_hodo_tot * err_reco_tot)**2 + (eff_reco_tot * err_hodo_tot)**2 ) if len(eff_reco_tot) > 0 else []

    eff_reco_mix = data_mix["recoeff"][mask_mix]
    err_reco_mix = data_mix["recoeff_error"][mask_mix]
    eff_hodo_mix = data_mix["hodoeff"][mask_mix]
    err_hodo_mix = data_mix["hodoeff_error"][mask_mix]
    eff_final_mix = eff_reco_mix * eff_hodo_mix
    err_final_mix = np.sqrt( (eff_hodo_mix * err_reco_mix)**2 + (eff_reco_mix * err_hodo_mix)**2 ) if len(eff_reco_mix) > 0 else []

    means = {
        "r_tot": get_weighted_mean_and_error(eff_reco_tot, err_reco_tot),
        "h_tot": get_weighted_mean_and_error(eff_hodo_tot, err_hodo_tot),
        "f_tot": get_weighted_mean_and_error(eff_final_tot, err_final_tot),
        "r_mix": get_weighted_mean_and_error(eff_reco_mix, err_reco_mix),
        "h_mix": get_weighted_mean_and_error(eff_hodo_mix, err_hodo_mix),
        "f_mix": get_weighted_mean_and_error(eff_final_mix, err_final_mix),
    }

    val_sig_eff, err_sig_eff = 0.0, 0.0
    diff_yield = N_tot - N_mix
    if diff_yield != 0:
        numerator = (N_tot * means["f_tot"][0]) - (N_mix * means["f_mix"][0])
        val_sig_eff = numerator / diff_yield
        term1 = (N_tot * means["f_tot"][1])**2
        term2 = (N_mix * means["f_mix"][1])**2
        err_sig_eff = (1.0 / diff_yield) * np.sqrt(term1 + term2)

    return N_tot, N_mix, sum_mass_tot, sum_mass_mix, val_sig_eff, err_sig_eff, diff_yield, means

def fill_histograms(hists, root_x, root_y, N_tot, N_mix, m_center, x_center, sig_eff, err_sig_eff, diff_yield, means, final_centroid):
    err_N_tot = np.sqrt(N_tot)
    err_N_mix = np.sqrt(N_mix)

    val_corr_yield, err_corr_stat, err_corr_sys, err_corr_total = 0.0, 0.0, 0.0, 0.0
    if sig_eff > 0 and diff_yield > 0:
        val_corr_yield = diff_yield / sig_eff
        err_corr_stat = np.sqrt(N_tot + N_mix) / sig_eff
        err_corr_sys = val_corr_yield * (err_sig_eff / sig_eff)
        err_corr_total = np.sqrt(err_corr_stat**2 + err_corr_sys**2)

    def set_bin(key, val, err):
        hists[key].SetBinContent(root_x, root_y, val)
        hists[key].SetBinError(root_x, root_y, err)
        if "stat" not in key and "sys" not in key and "Centroid" not in key:
            add_latex_to_bin(hists[key], m_center, x_center, val, err)

    set_bin("Y_total", N_tot, err_N_tot)
    set_bin("Y_mix", N_mix, err_N_mix)
    set_bin("E_total_reco", means["r_tot"][0], means["r_tot"][1])
    set_bin("E_mix_reco", means["r_mix"][0], means["r_mix"][1])
    set_bin("E_total_hodo", means["h_tot"][0], means["h_tot"][1])
    set_bin("E_mix_hodo", means["h_mix"][0], means["h_mix"][1])
    set_bin("E_total_final", means["f_tot"][0], means["f_tot"][1])
    set_bin("E_mix_final", means["f_mix"][0], means["f_mix"][1])
    set_bin("E_final_signal", sig_eff, err_sig_eff)
    set_bin("Y_corrected", val_corr_yield, err_corr_total)
    set_bin("Y_corrected_stat", val_corr_yield, err_corr_stat)
    set_bin("Y_corrected_sys", val_corr_yield, err_corr_sys)
    
    hists["Mass_Centroid"].SetBinContent(root_x, root_y, final_centroid)


# ==========================================
# 5. Combined Target Processing
# ==========================================
def process_all_datasets(file_lh2, file_ld2, file_fl, output_file):
    print("\nExtracting all data sets into memory to compute correct combined mass centroids...")
    
    # Ensure plots directory exists
    os.makedirs("./MassBinCentroids", exist_ok=True)
    
    f_lh2 = uproot.open(file_lh2)
    f_ld2 = uproot.open(file_ld2)
    f_fl = uproot.open(file_fl)

    data_lh2_tot = e906_chuck_cuts(f_lh2["result"])
    data_lh2_mix = e906_chuck_cuts(f_lh2["result_mix"])
    data_ld2_tot = e906_chuck_cuts(f_ld2["result"])
    data_ld2_mix = e906_chuck_cuts(f_ld2["result_mix"])
    data_fl_tot = e906_chuck_cuts(f_fl["result"])
    data_fl_mix = e906_chuck_cuts(f_fl["result_mix"])

    output_file.cd()
    hists_lh2 = create_histograms("LH2")
    hists_ld2 = create_histograms("LD2")
    hists_fl = create_histograms("Flask")

    # Colors for the 1D Mass Plots
    colors = {
        "LD2_tot": ROOT.kRed, "LD2_mix": ROOT.kBlue,
        "LH2_tot": ROOT.kRed,  "LH2_mix": ROOT.kBlue,
        "Fl_tot": ROOT.kGreen+2, "Fl_mix": ROOT.kOrange+1
    }

    for i_x in range(len(xf_bins_np) - 1):
        x_low, x_high = xf_bins_np[i_x], xf_bins_np[i_x+1]
        x_center = (x_low + x_high) / 2.0
        root_y = i_x + 1

        csv_filename_lh2 = f"Table_Kinematics_LH2_xF_{x_low:.2f}_{x_high:.2f}.csv"
        csv_filename_ld2 = f"Table_Kinematics_LD2_xF_{x_low:.2f}_{x_high:.2f}.csv"
        
        csv_rows_lh2 = []
        csv_rows_ld2 = []

        h1_ld2_tot = ROOT.TH1D(f"h1_ld2_tot_{i_x}", "", len(mass_bins_np)-1, mass_bins_np)
        h1_ld2_mix = ROOT.TH1D(f"h1_ld2_mix_{i_x}", "", len(mass_bins_np)-1, mass_bins_np)
        h1_lh2_tot = ROOT.TH1D(f"h1_lh2_tot_{i_x}", "", len(mass_bins_np)-1, mass_bins_np)
        h1_lh2_mix = ROOT.TH1D(f"h1_lh2_mix_{i_x}", "", len(mass_bins_np)-1, mass_bins_np)
        h1_fl_tot = ROOT.TH1D(f"h1_fl_tot_{i_x}", "", len(mass_bins_np)-1, mass_bins_np)
        h1_fl_mix = ROOT.TH1D(f"h1_fl_mix_{i_x}", "", len(mass_bins_np)-1, mass_bins_np)

        h1s = {"LD2_tot": h1_ld2_tot, "LD2_mix": h1_ld2_mix, "LH2_tot": h1_lh2_tot, 
               "LH2_mix": h1_lh2_mix, "Fl_tot": h1_fl_tot, "Fl_mix": h1_fl_mix}
        for k, h in h1s.items():
            h.SetLineColor(colors[k])
            h.SetLineWidth(2)
        
        latex_draws_lh2 = []
        latex_draws_ld2 = []
        
        max_label_y_lh2 = 0.5
        max_label_y_ld2 = 0.5

        for i_m in range(len(mass_bins_np) - 1):
            m_low, m_high = mass_bins_np[i_m], mass_bins_np[i_m+1]
            m_center = (m_low + m_high) / 2.0
            root_x = i_m + 1

            N_l2_t, N_l2_m, M_l2_t, M_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2 = extract_target_stats(data_ld2_tot, data_ld2_mix, m_low, m_high, x_low, x_high)
            N_lh_t, N_lh_m, M_lh_t, M_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh = extract_target_stats(data_lh2_tot, data_lh2_mix, m_low, m_high, x_low, x_high)
            N_fl_t, N_fl_m, M_fl_t, M_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl = extract_target_stats(data_fl_tot, data_fl_mix, m_low, m_high, x_low, x_high)

            h1_ld2_tot.SetBinContent(root_x, N_l2_t)
            h1_ld2_mix.SetBinContent(root_x, N_l2_m)
            h1_lh2_tot.SetBinContent(root_x, N_lh_t)
            h1_lh2_mix.SetBinContent(root_x, N_lh_m)
            h1_fl_tot.SetBinContent(root_x, N_fl_t)
            h1_fl_mix.SetBinContent(root_x, N_fl_m)

            if eps_lh2 > 0:
                cY_lh_t, cM_lh_t = N_lh_t / eps_lh2, M_lh_t / eps_lh2
                cY_lh_m, cM_lh_m = N_lh_m / eps_lh2, M_lh_m / eps_lh2
            else:
                cY_lh_t, cM_lh_t, cY_lh_m, cM_lh_m = 0.0, 0.0, 0.0, 0.0

            if eps_ld2 > 0:
                cY_l2_t, cM_l2_t = N_l2_t / eps_ld2, M_l2_t / eps_ld2
                cY_l2_m, cM_l2_m = N_l2_m / eps_ld2, M_l2_m / eps_ld2
            else:
                cY_l2_t, cM_l2_t, cY_l2_m, cM_l2_m = 0.0, 0.0, 0.0, 0.0

            if eps_fl > 0:
                cY_fl_t_lh, cM_fl_t_lh = (FLASK_NORM_LH2 * N_fl_t) / eps_fl, (FLASK_NORM_LH2 * M_fl_t) / eps_fl
                cY_fl_m_lh, cM_fl_m_lh = (FLASK_NORM_LH2 * N_fl_m) / eps_fl, (FLASK_NORM_LH2 * M_fl_m) / eps_fl
                cY_fl_t_l2, cM_fl_t_l2 = (FLASK_NORM_LD2 * N_fl_t) / eps_fl, (FLASK_NORM_LD2 * M_fl_t) / eps_fl
                cY_fl_m_l2, cM_fl_m_l2 = (FLASK_NORM_LD2 * N_fl_m) / eps_fl, (FLASK_NORM_LD2 * M_fl_m) / eps_fl
            else:
                cY_fl_t_lh, cM_fl_t_lh, cY_fl_m_lh, cM_fl_m_lh = 0.0, 0.0, 0.0, 0.0
                cY_fl_t_l2, cM_fl_t_l2, cY_fl_m_l2, cM_fl_m_l2 = 0.0, 0.0, 0.0, 0.0

            num_LH2 = (cM_lh_t - cM_lh_m) - (cM_fl_t_lh - cM_fl_m_lh)
            den_LH2 = (cY_lh_t - cY_lh_m) - (cY_fl_t_lh - cY_fl_m_lh)
            cent_LH2 = num_LH2 / den_LH2 if den_LH2 != 0 else m_center

            num_LD2 = (cM_l2_t - cM_l2_m) - (cM_fl_t_l2 - cM_fl_m_l2) - (LH2_TO_LD2_NORM * num_LH2)
            den_LD2 = (cY_l2_t - cY_l2_m) - (cY_fl_t_l2 - cY_fl_m_l2) - (LH2_TO_LD2_NORM * den_LH2)
            cent_LD2 = num_LD2 / den_LD2 if den_LD2 != 0 else m_center

            fill_histograms(hists_ld2, root_x, root_y, N_l2_t, N_l2_m, m_center, x_center, eps_ld2, err_eps_ld2, diff_l2, mean_l2, cent_LD2)
            fill_histograms(hists_lh2, root_x, root_y, N_lh_t, N_lh_m, m_center, x_center, eps_lh2, err_eps_lh2, diff_lh, mean_lh, cent_LH2)
            num_fl = (M_fl_t - M_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
            den_fl = (N_fl_t - N_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
            cent_fl = num_fl / den_fl if den_fl != 0 else m_center
            fill_histograms(hists_fl, root_x, root_y, N_fl_t, N_fl_m, m_center, x_center, eps_fl, err_eps_fl, diff_fl, mean_fl, cent_fl)

            jitter_x = [-0.05, 0.05]
            
            comp_data_lh2 = [
                ("LH2_tot", cM_lh_t, N_lh_t), ("LH2_mix", cM_lh_m, N_lh_m),
                ("Fl_tot", cM_fl_t_lh, N_fl_t), ("Fl_mix", cM_fl_m_lh, N_fl_m)
            ]
            valid_lh2 = [(lbl, cM, N) for lbl, cM, N in comp_data_lh2 if N > 0]
            valid_lh2.sort(key=lambda x: x[2])
            
            last_y = 1e-9
            for j, (lbl, cM, N) in enumerate(valid_lh2):
                desired_y = N * 1.5
                actual_y = max(desired_y, last_y * 2.5)
                x_pos = m_center + jitter_x[j % 2]
                txt = ROOT.TLatex(x_pos, actual_y, f"#splitline{{{cM:.2f}}}{{({N})}}")
                txt.SetTextSize(0.025)
                txt.SetTextColor(colors[lbl])
                txt.SetTextAlign(22)
                latex_draws_lh2.append(txt)
                last_y = actual_y
                max_label_y_lh2 = max(max_label_y_lh2, actual_y)


            comp_data_ld2 = [
                ("LD2_tot", cM_l2_t, N_l2_t), ("LD2_mix", cM_l2_m, N_l2_m),
                ("Fl_tot", cM_fl_t_l2, N_fl_t), ("Fl_mix", cM_fl_m_l2, N_fl_m)
            ]
            valid_ld2 = [(lbl, cM, N) for lbl, cM, N in comp_data_ld2 if N > 0]
            valid_ld2.sort(key=lambda x: x[2])
            
            last_y = 1e-9
            for j, (lbl, cM, N) in enumerate(valid_ld2):
                desired_y = N * 1.5
                actual_y = max(desired_y, last_y * 2.5)
                x_pos = m_center + jitter_x[j % 2]
                txt = ROOT.TLatex(x_pos, actual_y, f"#splitline{{{cM:.2f}}}{{({N})}}")
                txt.SetTextSize(0.025)
                txt.SetTextColor(colors[lbl])
                txt.SetTextAlign(22)
                latex_draws_ld2.append(txt)
                last_y = actual_y
                max_label_y_ld2 = max(max_label_y_ld2, actual_y)

            # Record CSV data
            csv_rows_lh2.append({
                "Mass Bin": f"[{m_low:.2f}, {m_high:.2f})",
                "Mass Center": m_center,
                "LH2 Mass Bin Average": cent_LH2,
                "N_LH2_total": N_lh_t,
                "N_LH2_mixed": N_lh_m,
                "N_flask_total": N_fl_t,
                "N_flask_mixed": N_fl_m,
                "eps_LH2": eps_lh2,
                "eps_LD2": eps_ld2,
                "eps_flask": eps_fl,
                "Corrected Total Mass LH2 (num)": num_LH2,
                "Corrected Yield LH2 (denom)": den_LH2
            })

            csv_rows_ld2.append({
                "Mass Bin": f"[{m_low:.2f}, {m_high:.2f})",
                "Mass Center": m_center,
                "LD2 Mass Bin Average": cent_LD2,
                "N_LD2_total": N_l2_t,
                "N_LD2_mixed": N_l2_m,
                "N_LH2_total": N_lh_t,
                "N_LH2_mixed": N_lh_m,
                "N_flask_total": N_fl_t,
                "N_flask_mixed": N_fl_m,
                "eps_LH2": eps_lh2,
                "eps_LD2": eps_ld2,
                "eps_flask": eps_fl,
                "Corrected Total Mass LD2 (num)": num_LD2,
                "Corrected Yield LD2 (denom)": den_LD2
            })

        # === Draw 1D Mass Plot for LH2 ===
        c_lh2 = ROOT.TCanvas(f"c_1d_mass_lh2_{i_x}", "", 1000, 700)
        c_lh2.SetRightMargin(0.05)
        c_lh2.SetLogy() 
        c_lh2.SetTickx(1) 
        c_lh2.SetTicky(1) 
        
        max_y_hist_lh2 = max([h1_lh2_tot.GetMaximum(), h1_lh2_mix.GetMaximum(), h1_fl_tot.GetMaximum(), h1_fl_mix.GetMaximum()])
        overall_max_lh2 = max(max_y_hist_lh2, max_label_y_lh2)
        if overall_max_lh2 <= 0: overall_max_lh2 = 1.0
        h1_lh2_tot.SetMaximum(overall_max_lh2 * 10.0) 
        h1_lh2_tot.SetMinimum(0.5) 
        
        h1_lh2_tot.SetTitle(f"Mass Distributions (LH2) {x_low:.2f} <= x_{{F}} < {x_high:.2f};Mass [GeV];Counts")
        h1_lh2_tot.GetXaxis().CenterTitle(True)
        h1_lh2_tot.GetYaxis().CenterTitle(True)
        h1_lh2_tot.Draw("HIST")
        h1_lh2_mix.Draw("HIST SAME")
        h1_fl_tot.Draw("HIST SAME")
        h1_fl_mix.Draw("HIST SAME")
        
        for l in latex_draws_lh2: l.Draw()
        
        leg_lh2 = ROOT.TLegend(0.75, 0.72, 0.93, 0.88)
        leg_lh2.SetBorderSize(0)
        leg_lh2.AddEntry(h1_lh2_tot, "LH2 Total", "l")
        leg_lh2.AddEntry(h1_lh2_mix, "LH2 Mix", "l")
        leg_lh2.AddEntry(h1_fl_tot,  "Flask Total", "l")
        leg_lh2.AddEntry(h1_fl_mix,  "Flask Mix", "l")
        leg_lh2.Draw()
        
        c_lh2.SaveAs(f"./MassBinCentroids/MassDist_LH2_xF_{x_low:.2f}_{x_high:.2f}.pdf")
        c_lh2.Close()

        # === Draw 1D Mass Plot for LD2 ===
        c_ld2 = ROOT.TCanvas(f"c_1d_mass_ld2_{i_x}", "", 1000, 700)
        c_ld2.SetRightMargin(0.05)
        c_ld2.SetLogy() 
        c_ld2.SetTickx(1) 
        c_ld2.SetTicky(1) 
        
        max_y_hist_ld2 = max([h1_ld2_tot.GetMaximum(), h1_ld2_mix.GetMaximum(), h1_fl_tot.GetMaximum(), h1_fl_mix.GetMaximum()])
        overall_max_ld2 = max(max_y_hist_ld2, max_label_y_ld2)
        if overall_max_ld2 <= 0: overall_max_ld2 = 1.0
        h1_ld2_tot.SetMaximum(overall_max_ld2 * 10.0) 
        h1_ld2_tot.SetMinimum(0.5)
        
        h1_ld2_tot.SetTitle(f"Mass Distributions (LD2) {x_low:.2f} <= x_{{F}} < {x_high:.2f};Mass [GeV];Counts")
        h1_ld2_tot.GetXaxis().CenterTitle(True)
        h1_ld2_tot.GetYaxis().CenterTitle(True)
        h1_ld2_tot.Draw("HIST")
        h1_ld2_mix.Draw("HIST SAME")
        h1_fl_tot.Draw("HIST SAME")
        h1_fl_mix.Draw("HIST SAME")
        
        for l in latex_draws_ld2: l.Draw()
        
        leg_ld2 = ROOT.TLegend(0.75, 0.72, 0.93, 0.88)
        leg_ld2.SetBorderSize(0)
        leg_ld2.AddEntry(h1_ld2_tot, "LD2 Total", "l")
        leg_ld2.AddEntry(h1_ld2_mix, "LD2 Mix", "l")
        leg_ld2.AddEntry(h1_fl_tot,  "Flask Total", "l")
        leg_ld2.AddEntry(h1_fl_mix,  "Flask Mix", "l")
        leg_ld2.Draw()
        
        c_ld2.SaveAs(f"./MassBinCentroids/MassDist_LD2_xF_{x_low:.2f}_{x_high:.2f}.pdf")
        c_ld2.Close()

        with open(csv_filename_lh2, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=csv_rows_lh2[0].keys())
            writer.writeheader()
            writer.writerows(csv_rows_lh2)
            
        with open(csv_filename_ld2, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=csv_rows_ld2[0].keys())
            writer.writeheader()
            writer.writerows(csv_rows_ld2)
            
        print(f"  Saved LH2 & LD2 1D Mass Dists and individual CSV Tables for xF: [{x_low:.2f}, {x_high:.2f})")

    # Save 2D PDFs
    def save_2d_pdfs(hists_dict, label):
        for name, hist in hists_dict.items():
            if "stat" in name or "sys" in name or "Centroid" in name: continue 
            c = ROOT.TCanvas(f"c_{name}_{label}", "", 1200, 900)
            c.SetRightMargin(0.15); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12)
            c.SetTickx(1); c.SetTicky(1)
            
            # Ensure 0.0% bins remain white in COLZ when there are negative bins
            if hist.GetMinimum() < 0:
                hist.SetMinimum(0)
                
            hist.Draw("COLZ")
            c.SaveAs(f"{name}_{label}.pdf")
            c.Close()

    print("Saving 2D Histogram PDFs...")
    save_2d_pdfs(hists_lh2, "LH2")
    save_2d_pdfs(hists_ld2, "LD2")
    save_2d_pdfs(hists_fl, "Flask")

    return hists_lh2, hists_ld2, hists_fl


# ==========================================
# 6. Subtracted Yield Calculations
# ==========================================
def generate_subtracted_plot(hists_target, hists_flask, flask_norm, target_label, output_file):
    print(f"\nGenerating Subtracted Plot: Y_corrected_{target_label} - Y_corrected_Flask...")
    
    if "Y_corrected_stat" not in hists_target or "Y_corrected_stat" not in hists_flask:
        print("Error: Y_corrected_stat histogram missing.")
        return None

    output_file.cd()
    
    h_target_stat = hists_target["Y_corrected_stat"]
    h_target_sys = hists_target["Y_corrected_sys"]
    h_flask_stat = hists_flask["Y_corrected_stat"]
    h_flask_sys = hists_flask["Y_corrected_sys"]
    
    name = f"Y_corrected_Subtracted_{target_label}"
    title = f"Corrected Yield ({target_label} - Flask)"
    h_sub = ROOT.TH2D(name, f"{title};Mass [GeV];x_{{F}}", 
                      len(mass_bins_np)-1, mass_bins_np, 
                      len(xf_bins_np)-1, xf_bins_np)
    h_sub.Sumw2()
    h_sub.SetStats(0)
    h_sub.GetXaxis().CenterTitle()
    h_sub.GetYaxis().CenterTitle()
    h_sub.GetXaxis().SetTitleOffset(1.2)
    h_sub.GetYaxis().SetTitleOffset(1.2)

    h_sub_stat = h_sub.Clone(f"{name}_stat")
    h_sub_sys = h_sub.Clone(f"{name}_sys")
    
    h_sub_centroid = hists_target["Mass_Centroid"].Clone(f"{name}_Mass_Centroid")

    for i_m in range(len(mass_bins_np) - 1):
        m_center = (mass_bins_np[i_m] + mass_bins_np[i_m+1]) / 2.0
        bin_x = i_m + 1
        
        for i_x in range(len(xf_bins_np) - 1):
            x_center = (xf_bins_np[i_x] + xf_bins_np[i_x+1]) / 2.0
            bin_y = i_x + 1
            
            y_target = h_target_stat.GetBinContent(bin_x, bin_y)
            e_target_stat = h_target_stat.GetBinError(bin_x, bin_y)
            e_target_sys = h_target_sys.GetBinError(bin_x, bin_y)
            
            y_flask = h_flask_stat.GetBinContent(bin_x, bin_y)
            e_flask_stat = h_flask_stat.GetBinError(bin_x, bin_y)
            e_flask_sys = h_flask_sys.GetBinError(bin_x, bin_y)
            
            val_sub = y_target - (flask_norm * y_flask)
            
            err_sub_stat = np.sqrt(e_target_stat**2 + (flask_norm * e_flask_stat)**2)
            err_sub_sys = np.sqrt(e_target_sys**2 + (flask_norm * e_flask_sys)**2)
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
        h_sub.SetMinimum(0.0 if local_min < 0 else local_min)
        h_sub.SetMaximum(local_max)
    else:
        h_sub.SetMinimum(0.0)
        h_sub.SetMaximum(1.0)
        
    h_sub.Draw("COLZ")
    c.SaveAs(f"{name}.pdf")
    c.Close()
    print(f"  Subtracted plot saved for {target_label}.")
    
    return {"stat": h_sub_stat, "sys": h_sub_sys, "centroid": h_sub_centroid}

def generate_pd_subtracted_plot(hists_ld2, hists_lh2, hists_flask, output_file):
    print("\nGenerating Subtracted Plot for LD2 pd formula...")
    output_file.cd()

    h_ld2_stat = hists_ld2["Y_corrected_stat"]
    h_ld2_sys = hists_ld2["Y_corrected_sys"]
    
    h_lh2_stat = hists_lh2["Y_corrected_stat"]
    h_lh2_sys = hists_lh2["Y_corrected_sys"]
    
    h_flask_stat = hists_flask["Y_corrected_stat"]
    h_flask_sys = hists_flask["Y_corrected_sys"]

    c_lh2 = THD_THH_RATIO * (PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_LH2)
    c_flask_sub = FLASK_NORM_LD2 - (c_lh2 * FLASK_NORM_LH2)

    name = "Y_corrected_Subtracted_LD2"
    title = "Corrected Yield (LD2 - LH2 - Flask)"
    h_pd = ROOT.TH2D(name, f"{title};Mass [GeV];x_{{F}}", 
                      len(mass_bins_np)-1, mass_bins_np, 
                      len(xf_bins_np)-1, xf_bins_np)
    h_pd.Sumw2()
    h_pd.SetStats(0)
    h_pd.GetXaxis().CenterTitle()
    h_pd.GetYaxis().CenterTitle()
    h_pd.GetXaxis().SetTitleOffset(1.2)
    h_pd.GetYaxis().SetTitleOffset(1.2)

    h_pd_stat = h_pd.Clone(f"{name}_stat")
    h_pd_sys = h_pd.Clone(f"{name}_sys")
    
    h_pd_centroid = hists_ld2["Mass_Centroid"].Clone(f"{name}_Mass_Centroid")

    for i_m in range(len(mass_bins_np) - 1):
        m_center = (mass_bins_np[i_m] + mass_bins_np[i_m+1]) / 2.0
        bin_x = i_m + 1
        
        for i_x in range(len(xf_bins_np) - 1):
            x_center = (xf_bins_np[i_x] + xf_bins_np[i_x+1]) / 2.0
            bin_y = i_x + 1
            
            y_ld2 = h_pd_stat.GetBinContent(bin_x, bin_y) if 'h_pd_stat' in locals() else h_ld2_stat.GetBinContent(bin_x, bin_y)
            y_lh2 = h_lh2_stat.GetBinContent(bin_x, bin_y)
            y_flask = h_flask_stat.GetBinContent(bin_x, bin_y)
            
            y_ld2 = h_ld2_stat.GetBinContent(bin_x, bin_y)

            e_ld2_stat = h_ld2_stat.GetBinError(bin_x, bin_y)
            e_lh2_stat = h_lh2_stat.GetBinError(bin_x, bin_y)
            e_flask_stat = h_flask_stat.GetBinError(bin_x, bin_y)
            
            e_ld2_sys = h_ld2_sys.GetBinError(bin_x, bin_y)
            e_lh2_sys = h_lh2_sys.GetBinError(bin_x, bin_y)
            e_flask_sys = h_flask_sys.GetBinError(bin_x, bin_y)
            
            val_pd = y_ld2 - (c_lh2 * y_lh2) - (c_flask_sub * y_flask)
            
            err_pd_stat = np.sqrt(e_ld2_stat**2 + (c_lh2 * e_lh2_stat)**2 + (c_flask_sub * e_flask_stat)**2)
            err_pd_sys = np.sqrt(e_ld2_sys**2 + (c_lh2 * e_lh2_sys)**2 + (c_flask_sub * e_flask_sys)**2)
            err_pd_total = np.sqrt(err_pd_stat**2 + err_pd_sys**2)
            
            h_pd.SetBinContent(bin_x, bin_y, val_pd)
            h_pd.SetBinError(bin_x, bin_y, err_pd_total)
            add_latex_to_bin(h_pd, m_center, x_center, val_pd, err_pd_total)

            h_pd_stat.SetBinContent(bin_x, bin_y, val_pd)
            h_pd_stat.SetBinError(bin_x, bin_y, err_pd_stat)

            h_pd_sys.SetBinContent(bin_x, bin_y, val_pd)
            h_pd_sys.SetBinError(bin_x, bin_y, err_pd_sys)

    c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTickx(1)
    c.SetTicky(1)
    
    local_min = sys.float_info.max
    local_max = -sys.float_info.max
    has_data = False
    for ix in range(1, h_pd.GetNbinsX() + 1):
        for iy in range(1, h_pd.GetNbinsY() + 1):
            val = h_pd.GetBinContent(ix, iy)
            if val != 0.0:
                has_data = True
                if val < local_min: local_min = val
                if val > local_max: local_max = val
    if has_data:
        h_pd.SetMinimum(0.0 if local_min < 0 else local_min)
        h_pd.SetMaximum(local_max)
    else:
        h_pd.SetMinimum(0.0)
        h_pd.SetMaximum(1.0)
        
    h_pd.Draw("COLZ")
    c.SaveAs(f"{name}.pdf")
    c.Close()
    print("  Subtracted plot saved for LD2 pd equation.")

    return {"stat": h_pd_stat, "sys": h_pd_sys, "centroid": h_pd_centroid}


# ==========================================
# 7. Cross Section Calculation & Plotting
# ==========================================
def calculate_and_plot_cross_section(h_sub_dict, target_label, global_constant, output_file):
    print(f"\nCalculating and Plotting Double Differential Cross-Sections for {target_label}...")
    acc_path = "acceptance_mass_xF.root"
    psip_path = "All_PsiP_Contaminations.root" 
    
    if target_label == "LD2":
        theory_ct18_path = "CT18_xFnew_d_1sigma.root"
        theory_nnpdf_path = "NNPDF40_xFnew_d.root"
    else:
        theory_ct18_path = "CT18_xFnew_p_1sigma.root"
        theory_nnpdf_path = "NNPDF40_xFnew_p.root"
        
    if not os.path.exists(acc_path):
        raise FileNotFoundError(f"CRITICAL ERROR: Acceptance file '{acc_path}' not found!")

    try:
        acc_file = ROOT.TFile.Open(acc_path)
        if not acc_file or acc_file.IsZombie():
            raise OSError(f"Acceptance file corrupted.")
            
        f_ct18 = ROOT.TFile.Open(theory_ct18_path) if os.path.exists(theory_ct18_path) else None
        f_nnpdf = ROOT.TFile.Open(theory_nnpdf_path) if os.path.exists(theory_nnpdf_path) else None
        f_psip = ROOT.TFile.Open(psip_path) if os.path.exists(psip_path) else None
            
    except Exception as e:
        print(f"Error opening files: {e}")
        sys.exit(1)

    output_file.cd()
    
    h_sub_stat = h_sub_dict["stat"]
    h_sub_sys = h_sub_dict["sys"]
    h_centroid = h_sub_dict["centroid"]

    n_xf_bins = len(xf_bins_np) - 1
    n_mass_bins = len(mass_bins_np) - 1
    
    latex_psip_table_content = r"""\begin{longtable}{|c|c|c|c|c|}
\caption{$\psi'$ Contamination Table for %s} \label{tab:psip_contamination_%s} \\
\hline
\textbf{xF bin} & \textbf{Mass bin (GeV)} & \textbf{$\psi'$ contamination} & \textbf{$\sigma \pm \delta \sigma^{stat.} \pm \delta \sigma^{syst.}$} & \textbf{$\delta\sigma_{\psi'}^{\rm syst.}$ (nb GeV$^2$)} \\
\hline
\endfirsthead

\multicolumn{5}{c}%%
{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
\hline
\textbf{xF bin} & \textbf{Mass bin (GeV)} & \textbf{$\psi'$ contamination} & \textbf{$\sigma \pm \delta \sigma^{stat.} \pm \delta \sigma^{syst.}$} & \textbf{$\delta\sigma_{\psi'}^{\rm syst.}$ (nb GeV$^2$)} \\
\hline
\endhead

\hline \multicolumn{5}{|r|}{{Continued on next page}} \\ \hline
\endfoot

\hline
\endlastfoot
""" % (target_label, target_label)

    # Helper function to avoid duplicating ROOT canvas boilerplate
    def draw_and_save_canvas(plot_type, g_xsec, g_sys, xf_min, xf_max, xf_bin_index, y_min_data, y_max_data):
        c_xsec = ROOT.TCanvas(f"c_xsec_{target_label}_{xf_bin_index}_{plot_type}", "", 800, 600)
        c_xsec.SetLogy(); c_xsec.SetTickx(1); c_xsec.SetTicky(1)
        
        mg = ROOT.TMultiGraph()
        mg.SetTitle(f";Mass [GeV];M^{{3}} d^{{2}}\\sigma / dM dx_{{F}} [nb GeV^{{2}}/Nucleus]")
        
        if y_min_data < y_max_data:
            mg.SetMinimum(y_min_data * 0.2)
            mg.SetMaximum(y_max_data * 5.0)
        else:
            mg.SetMinimum(1e-3)
            mg.SetMaximum(3.0)
        
        leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
        leg.SetBorderSize(0)
        gr_name = f"gr_xFbin{xf_bin_index}"
        
        if f_ct18:
            g_ct18 = f_ct18.Get(gr_name)
            if g_ct18:
                g_ct18_clone = g_ct18.Clone()
                g_ct18_clone.SetLineColor(ROOT.kGreen + 2); g_ct18_clone.SetFillColorAlpha(ROOT.kGreen - 5, 0.5); g_ct18_clone.SetFillStyle(3002)
                mg.Add(g_ct18_clone, "L3")
                leg.AddEntry(g_ct18_clone, "CT18 NLO", "lf") 
        
        if f_nnpdf:
            g_nnpdf = f_nnpdf.Get(gr_name)
            if g_nnpdf:
                g_nnpdf_clone = g_nnpdf.Clone()
                g_nnpdf_clone.SetLineColor(ROOT.kBlue + 2); g_nnpdf_clone.SetFillColorAlpha(ROOT.kAzure + 1, 0.5); g_nnpdf_clone.SetFillStyle(3002)
                mg.Add(g_nnpdf_clone, "L3")
                leg.AddEntry(g_nnpdf_clone, "NNPDF4.0 NLO", "lf") 
        
        g_sys.SetMarkerSize(0); g_sys.SetLineColor(ROOT.kRed); g_sys.SetFillColorAlpha(ROOT.kPink - 9, 0.5); g_sys.SetFillStyle(1001)
        mg.Add(g_sys, "2"); leg.AddEntry(g_sys, "Systematic Unc.", "f")

        g_xsec.SetMarkerStyle(20); g_xsec.SetMarkerColor(ROOT.kRed); g_xsec.SetLineColor(ROOT.kRed)
        mg.Add(g_xsec, "P"); leg.AddEntry(g_xsec, f"Data ({target_label})", "lep")
        
        mg.Draw("A"); mg.GetXaxis().CenterTitle(); mg.GetYaxis().CenterTitle(); mg.GetXaxis().SetLimits(4.2, 8.7)
        leg.Draw()

        target_prefix = "pp" if target_label == "LH2" else "pd" if target_label == "LD2" else target_label
        internal_title = ROOT.TLatex()
        internal_title.SetNDC(True); internal_title.SetTextFont(42); internal_title.SetTextSize(0.04); internal_title.SetTextAlign(13)
        internal_title.DrawLatex(0.14, 0.86, f"Drell-Yan process in {target_prefix} at {xf_min:.2f} #leq x_{{F}} < {xf_max:.2f}")

        prelim = ROOT.TLatex()
        prelim.SetNDC(True); prelim.SetTextColor(ROOT.kBlue); prelim.SetTextAlign(33); prelim.SetTextSize(0.05)
        prelim.DrawLatex(0.82, 0.6, "Preliminary")

        plot_y_min = y_min_data * 0.2 if y_min_data < y_max_data else 1e-3
        plot_y_max = y_max_data * 5.0 if y_min_data < y_max_data else 3.0
        
        log_min = math.log10(plot_y_min)
        log_max = math.log10(plot_y_max)
        dynamic_y = 10**(log_min + 0.05 * (log_max - log_min))

        lumi_note = ROOT.TLatex()
        lumi_note.SetNDC(False)
        lumi_note.SetTextFont(42)
        lumi_note.SetTextColor(ROOT.kBlack)
        lumi_note.SetTextAlign(11) 
        lumi_note.SetTextSize(0.025)
        lumi_note.DrawLatex(4.3, dynamic_y, "10% global uncertainty due to the integrated luminosity is not included in the error bands")
        
        plot_name = f"CrossSection_{target_label}_xF_{xf_min:.2f}_{xf_max:.2f}_{plot_type}.pdf"
        c_xsec.SaveAs(plot_name)
        
        g_xsec.Write(f"g_xsec_{target_label}_{xf_bin_index}_{plot_type}")
        g_sys.Write(f"g_sys_{target_label}_{xf_bin_index}_{plot_type}")
        c_xsec.Close()
        print(f"  Saved {plot_type} plot: {plot_name}")

    for i_x in range(n_xf_bins):
        xf_bin_index = i_x 
        root_bin_y = i_x + 1
        
        xf_min, xf_max = xf_bins_np[i_x], xf_bins_np[i_x+1]
        
        acc_hist_name = f"h_ratio_{target_label}_xF_bin{xf_bin_index}"
        h_acc = acc_file.Get(acc_hist_name)
        if not h_acc:
            fallback_name = f"h_ratio_LH2_xF_bin{xf_bin_index}"
            h_acc = acc_file.Get(fallback_name)

        if not h_acc:
            print(f"Warning: Acceptance histogram for xF bin {xf_bin_index} not found. Skipping bin.")
            continue
        
        h_ratio_psip = f_psip.Get(f"hRatio_PsiP_DY_xF_{xf_bin_index}") if f_psip else None
        
        print(f"\n--- {target_label} xF Bin {xf_bin_index}: {xf_min:.2f} < xF < {xf_max:.2f} ---")

        # 1. Setup Data-Driven Centroid Graphs
        g_xsec_cent = ROOT.TGraphErrors()
        g_xsec_cent.SetName(f"g_xsec_cent_{target_label}_xF_{xf_bin_index}")
        g_sys_cent = ROOT.TGraphErrors()
        g_sys_cent.SetName(f"g_sys_cent_{target_label}_xF_{xf_bin_index}")

        # 2. Setup Geometric Center Graphs
        g_xsec_geo = ROOT.TGraphErrors()
        g_xsec_geo.SetName(f"g_xsec_geo_{target_label}_xF_{xf_bin_index}")
        g_sys_geo = ROOT.TGraphErrors()
        g_sys_geo.SetName(f"g_sys_geo_{target_label}_xF_{xf_bin_index}")
        
        point_idx = 0
        y_min_data = sys.float_info.max
        y_max_data = -sys.float_info.max
        
        for i_m in range(n_mass_bins):
            root_bin_x = i_m + 1
            
            mass_min, mass_max = mass_bins_np[i_m], mass_bins_np[i_m+1]
            geometric_center = (mass_min + mass_max) / 2.0
            mass_width = mass_max - mass_min
            
            actual_mass_center = h_centroid.GetBinContent(root_bin_x, root_bin_y)
            if actual_mass_center < mass_min or actual_mass_center > mass_max:
                actual_mass_center = geometric_center 
            
            Y_sub = h_sub_stat.GetBinContent(root_bin_x, root_bin_y)
            Y_sub_stat_err = h_sub_stat.GetBinError(root_bin_x, root_bin_y)
            Y_sub_sys_err = h_sub_sys.GetBinError(root_bin_x, root_bin_y)
            
            if Y_sub <= 0: continue
            
            psip_ratio, psip_ratio_err = 0.0, 0.0
            if h_ratio_psip:
                ratio_bin = h_ratio_psip.FindBin(actual_mass_center)
                psip_ratio = h_ratio_psip.GetBinContent(ratio_bin)
                psip_ratio_err = h_ratio_psip.GetBinError(ratio_bin)

            if psip_ratio > 1.0: psip_ratio, psip_ratio_err = 0.0, 0.0

            acc_bin = h_acc.FindBin(actual_mass_center)
            acceptance = h_acc.GetBinContent(acc_bin)
            acceptance_err = h_acc.GetBinError(acc_bin)
            
            print(f"  Mass Centroid {actual_mass_center:.3f} GeV: Acceptance = {acceptance:.5f} +/- {acceptance_err:.5f}")
            if acceptance <= 0: continue
            
            numerator = global_constant * Y_sub
            denominator = mass_width * acceptance
            xsec = numerator / denominator
            
            stat_unc = (Y_sub_stat_err / Y_sub) * xsec
            sys_acc = (acceptance_err / acceptance) * xsec
            sys_tot_yield = (Y_sub_sys_err / Y_sub) * xsec
            sys_psip_cont = psip_ratio * xsec if Y_sub > 0 else 0.0
            total_systematic_unc = np.sqrt(sys_acc**2 + sys_tot_yield**2 + sys_psip_cont**2)
            
            scaling_factor = actual_mass_center**3
            scaled_xsec = xsec * scaling_factor
            scaled_stat_err = stat_unc * scaling_factor
            scaled_sys_err = total_systematic_unc * scaling_factor
            scaled_sys_psip = sys_psip_cont * scaling_factor
            
            if psip_ratio > 0.0:
                s_xf = f"[{xf_min:.2f}, {xf_max:.2f})"
                s_mass = f"[{mass_min:.2f}, {mass_max:.2f})"
                s_ratio = f"{psip_ratio:.4f} $\\pm$ {psip_ratio_err:.4f}"
                s_sigma_psip_col = f"{scaled_xsec:.4f} $\\pm$ {scaled_stat_err:.4f} $\\pm$ {scaled_sys_err:.4f}"
                latex_psip_table_content += f"{s_xf} & {s_mass} & {s_ratio} & {s_sigma_psip_col} & {scaled_sys_psip:.4f} \\\\ \n\\hline\n"
            
            # Populate Data-Driven Centroid Graphs
            g_sys_cent.SetPoint(point_idx, geometric_center, scaled_xsec) # The box spans the full bin width
            g_sys_cent.SetPointError(point_idx, mass_width/2.0, scaled_sys_err)
            g_xsec_cent.SetPoint(point_idx, actual_mass_center, scaled_xsec)
            g_xsec_cent.SetPointError(point_idx, 0.0, scaled_stat_err)

            # Populate Geometric Center Graphs
            g_sys_geo.SetPoint(point_idx, geometric_center, scaled_xsec)
            g_sys_geo.SetPointError(point_idx, mass_width/2.0, scaled_sys_err)
            g_xsec_geo.SetPoint(point_idx, geometric_center, scaled_xsec)
            g_xsec_geo.SetPointError(point_idx, 0.0, scaled_stat_err)
            
            max_err_for_range = max(scaled_stat_err, scaled_sys_err)
            y_high = scaled_xsec + max_err_for_range
            y_low  = scaled_xsec - max_err_for_range
            
            if y_low <= 0: y_low = scaled_xsec * 0.5 
            if y_high > y_max_data: y_max_data = y_high
            if y_low < y_min_data: y_min_data = y_low
            
            point_idx += 1
            
        if g_xsec_cent.GetN() > 0:
            # Draw and save original Centroid plot
            draw_and_save_canvas("Centroid", g_xsec_cent, g_sys_cent, xf_min, xf_max, xf_bin_index, y_min_data, y_max_data)
            
            # Draw and save new Geometric Center plot
            draw_and_save_canvas("GeoCenter", g_xsec_geo, g_sys_geo, xf_min, xf_max, xf_bin_index, y_min_data, y_max_data)
        else:
            print(f"  No valid points for xF bin {xf_bin_index}.")

    with open(f"Table_PsiP_Contamination_{target_label}.tex", "w") as f:
        f.write(latex_psip_table_content + r"\end{longtable}" + "\n")

    acc_file.Close()
    if f_ct18: f_ct18.Close()
    if f_nnpdf: f_nnpdf.Close()
    if f_psip: f_psip.Close()

# ==========================================
# 8. Overlay All Cross-Sections on One Canvas
# ==========================================
def scale_tgrapherrors(g, scale):
    """Helper to scale a TGraphErrors object's Y values and Y errors using buffer access."""
    if not g: return
    x_buf = g.GetX()
    y_buf = g.GetY()
    for i in range(g.GetN()):
        x = x_buf[i]
        y = y_buf[i]
        ex = g.GetErrorX(i)
        ey = g.GetErrorY(i)
        g.SetPoint(i, x, y * scale)
        g.SetPointError(i, ex, ey * scale)

def generate_overlay_plot(out_file, target_label, plot_type):
    print(f"\nGenerating full overlay plot for {target_label} ({plot_type})...")
    
    ROOT.gStyle.SetTitleAlign(23)
    ROOT.gStyle.SetTitleX(0.5)
    ROOT.gStyle.SetTitleY(0.99)
    ROOT.gStyle.SetTitleH(0.04)
    ROOT.gStyle.SetTitleBorderSize(0)

    # Canvas height is 1800 (2x taller than original 900)
    canvas = ROOT.TCanvas(f"canvas_overlay_{target_label}_{plot_type}", "Cross-Section Comparison", 1200, 1800)
    canvas.SetLogy()
    
    # Left margin at 0.15 to ensure the Y-axis label doesn't get cut off
    canvas.SetLeftMargin(0.15) 
    canvas.SetBottomMargin(0.12)

    # 1. Add Pad Ticks if it's the LD2 target overlay
    if target_label == "LD2":
        canvas.SetTickx(1)
        canvas.SetTicky(1)

    legend = ROOT.TLegend(0.75, 0.45, 0.9, 0.9)
    legend.SetHeader("x_{F} Bins")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(43)   # Pixel-based font
    legend.SetTextSize(18)   # Exactly 18 pixels tall

    colors = [
        ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kCyan,
        ROOT.kOrange + 7, ROOT.kSpring + 5, ROOT.kTeal + 5, ROOT.kAzure + 1,
        ROOT.kGray + 2, ROOT.kViolet - 5, ROOT.kYellow + 2, ROOT.kPink + 1,
        ROOT.kGreen - 9, ROOT.kRed - 9
    ]

    h_frame = canvas.DrawFrame(2.8, 1e-3, 9.0, 1e35)
    h_frame.SetTitle(f"DY Absolute Cross-Section Vs Mass for x_{{F}} bins ({target_label})")
    h_frame.GetXaxis().SetTitle("Invariant Mass (GeV)")
    h_frame.GetXaxis().CenterTitle()
    h_frame.GetXaxis().SetTitleOffset(1.2)
    h_frame.GetYaxis().SetTitle("M^{3} #frac{d^{2}#sigma}{dMdx_{F}} (nb GeV^{2}) #times 5 #times 10^{2 #times x_{F} bin index}")
    h_frame.GetYaxis().CenterTitle()
    h_frame.GetYaxis().SetTitleOffset(1.8)  
    
    latex_labels = []
    
    n_xf_bins = len(xf_bins_np) - 1
    # Iterate and extract graphs directly from the current output file
    for i in range(n_xf_bins):
        g_xsec = out_file.Get(f"g_xsec_{target_label}_{i}_{plot_type}")
        g_sys = out_file.Get(f"g_sys_{target_label}_{i}_{plot_type}")
        
        if not g_xsec or not g_sys:
            continue
            
        # Clone to avoid overwriting the original file contents
        g_xsec_clone = g_xsec.Clone(f"g_xsec_clone_{i}")
        g_sys_clone = g_sys.Clone(f"g_sys_clone_{i}")
        
        scale_factor = 5 * (10**(2*i))
        scale_tgrapherrors(g_xsec_clone, scale_factor)
        scale_tgrapherrors(g_sys_clone, scale_factor)

        color = colors[i % len(colors)]

        # Set visual styles
        g_sys_clone.SetLineColor(color)
        g_sys_clone.SetFillColorAlpha(color, 0.35)
        g_sys_clone.SetFillStyle(1001)
        g_sys_clone.SetMarkerSize(0)

        g_xsec_clone.SetLineColor(color)
        g_xsec_clone.SetMarkerColor(color)
        g_xsec_clone.SetMarkerStyle(ROOT.kFullCircle)
        g_xsec_clone.SetMarkerSize(1.0)

        # Draw graphs
        g_sys_clone.Draw("2 SAME")
        g_xsec_clone.Draw("P SAME")

        # Add TLatex label for the xF bin range
        low_edge, high_edge = xF_bin_ranges[i-2] if i -2 >= 0 else (xf_bins_np[i], xf_bins_np[i+1])
        y_pos = 1.0 * scale_factor
        
        # Draw in data coordinates (to align with the curves) but sized in pixels
        latex = ROOT.TLatex(3.1, y_pos, f"{low_edge:.2f} #leq x_{{F}} < {high_edge:.2f}")
        latex.SetTextFont(43)    # Pixel-based font
        latex.SetTextSize(20)    # Lock to 20 pixels
        latex.SetTextColor(color)
        latex.Draw()
        latex_labels.append(latex) 

        legend.AddEntry(g_xsec_clone, f"x_{{F}} bin {i}", "pl")

    # 2. Conditionally Draw Legend (Hide it for LD2)
    if target_label != "LD2":
        legend.Draw()

    if target_label == "LD2":
        # 3. Draw Blue Preliminary label on Top Right (leaving space for SeaQuest logo)
        prelim = ROOT.TLatex()
        prelim.SetNDC(True)       # Canvas percentages
        prelim.SetTextFont(43)    # Pixel-based font
        prelim.SetTextSize(32)    # Larger text for Preliminary
        prelim.SetTextColor(ROOT.kBlue)
        prelim.SetTextAlign(33)   # Right-aligned
        
        # Uncomment the line below to draw it when you are ready
        # prelim_drawn = prelim.DrawLatex(0.85, 0.70, "Preliminary")
        # latex_labels.append(prelim_drawn)

        # 4. Add the 10% global uncertainty text using NDC so it never floats away
        lumi_note = ROOT.TLatex()
        lumi_note.SetNDC(True)    # Canvas percentages
        lumi_note.SetTextFont(43) # Pixel-based font
        lumi_note.SetTextSize(24) # Lock to 16 pixels
        lumi_note.SetTextColor(ROOT.kBlack)
        lumi_note.SetTextAlign(11) # Left-aligned
        
        # Place firmly at 16.5% from left edge, 14% from bottom edge
        lumi_drawn = lumi_note.DrawLatex(0.165, 0.15, "10% global uncertainty due to the integrated luminosity is not included in the error bands")
        latex_labels.append(lumi_drawn)

    canvas.Update()
    
    out_pdf = f"cross_section_overlay_{target_label}_{plot_type}.pdf"
    canvas.SaveAs(out_pdf)
    
    # Save the canvas to the active root file so you can adjust it later
    out_file.cd()
    canvas.Write(f"canvas_overlay_{target_label}_{plot_type}")
    print(f"  Overlay plot saved to {out_pdf} and written to output ROOT file.")
# ==========================================
# 9. Latex Appendix Generation
# ==========================================
def generate_latex_appendix():
    print("\nGenerating LaTeX Appendix...")
    latex_filename = "Appendix_MassCentroids.tex"
    
    # NOTE: Set this to the directory where the plots are generated relative to your main.tex
    IMAGE_PATH_PREFIX = "./MassBinCentroids/" 
    
    with open(latex_filename, "w") as tex_file:
        
        # Using a raw string (r"") prevents Python from interpreting \b in \begin or \f in \frac as escape characters.
        intro_text = r"""
\section{Appendix: Determination of Mass Bin Centroids}
\label{sec:appendix_mass_centroids}

In order to accurately plot the double-differential cross-sections, the horizontal placement of the data points must reflect the true physical average of the invariant mass within each finite mass bin, rather than the simple geometric center. Because the dimuon mass distribution falls exponentially, the true centroid is typically skewed toward the lower edge of the bin.

To account for combinatoric backgrounds, efficiency corrections, and flask contributions, the centroid calculation must be weighted using the fully corrected equations for the yields. 

For the LH$_2$ target, the corrected yield is given by:
\begin{equation}
    \label{eq:appendix_yield_lh2}
    Y_{\text{corrected}}^{\text{LH}_2} = \frac{Y^{\text{LH}_2}_{\text{total}} - Y^{\text{LH}_2}_{\text{mixed}}}{\langle\epsilon^{\text{LH}_2}_{\text{signal}}\rangle} 
        - \frac{I_{\text{LH}_2}}{I_{\text{flask}}} 
        \left( \frac{Y^{\text{flask}}_{\text{total}} - Y^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right)
\end{equation}

For the LD$_2$ target, the corrected yield accounts for both the empty flask and the hydrogen contribution:
\begin{align}
    \label{eq:appendix_yield_ld2}
    Y_{\text{corrected}}^{\text{LD}_2} &= 
        \frac{Y^{\text{LD}_2}_{\text{total}} - Y^{\text{LD}_2}_{\text{mixed}}}{\langle\epsilon^{\text{LD}_2}_{\text{signal}}\rangle} 
        - \frac{I_{\text{LD}_2}}{I_{\text{flask}}} 
        \left( \frac{Y^{\text{flask}}_{\text{total}} - Y^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right) \nonumber \\
    &\quad - \frac{T_{\text{HD}}}{T_{\text{HH}}} \frac{I_{\text{LD}_2}}{I_{\text{LH}_2}}\left[ 
        \frac{Y^{\text{LH}_2}_{\text{total}} - Y^{\text{LH}_2}_{\text{mixed}}}{\langle\epsilon^{\text{LH}_2}_{\text{signal}}\rangle} 
        - \frac{I_{\text{LH}_2}}{I_{\text{flask}}} 
        \left( \frac{Y^{\text{flask}}_{\text{total}} - Y^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right)
    \right]
\end{align}

The mass bin centroid for the LH$_2$ target is then calculated as the efficiency-corrected, background-subtracted mass sum divided by the corrected yield:
\begin{equation}
    \langle M \rangle_{\text{LH}_2} = \frac{1}{Y^{\text{LH}_2}_{\text{corrected}}} \left[ \frac{\sum M^{\text{LH}_2}_{\text{total}} - \sum M^{\text{LH}_2}_{\text{mixed}}}{\langle\epsilon^{\text{LH}_2}_{\text{signal}}\rangle} - \frac{I_{\text{LH}_2}}{I_{\text{flask}}}\left( \frac{\sum M^{\text{flask}}_{\text{total}} - \sum M^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right) \right]
\end{equation}

Similarly, we calculate the mass bin centroid for the LD$_2$ target using the fully corrected yield and mass sums:
\begin{align}
    \langle M \rangle_{\text{LD}_2} &= \frac{1}{Y^{\text{LD}_2}_{\text{corrected}}} \Bigg\{ \frac{\sum M^{\text{LD}_2}_{\text{total}} - \sum M^{\text{LD}_2}_{\text{mixed}}}{\langle\epsilon^{\text{LD}_2}_{\text{signal}}\rangle} - \frac{I_{\text{LD}_2}}{I_{\text{flask}}}\left( \frac{\sum M^{\text{flask}}_{\text{total}} - \sum M^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right) \nonumber \\
    &\quad - \frac{T_{\text{HD}}}{T_{\text{HH}}} \frac{I_{\text{LD}_2}}{I_{\text{LH}_2}}\left[\frac{\sum M^{\text{LH}_2}_{\text{total}} - \sum M^{\text{LH}_2}_{\text{mixed}}}{\langle\epsilon^{\text{LH}_2}_{\text{signal}}\rangle} - \frac{I_{\text{LH}_2}}{I_{\text{flask}}}\left( \frac{\sum M^{\text{flask}}_{\text{total}} - \sum M^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right)\right] \Bigg\}
\end{align}

Here, we plot dimuon yield distributions for LH$_2$ cross-section calculations for the following components:
\begin{itemize}
    \item LH$_2$ total (red)
    \item LH$_2$ mixed (blue)
    \item Flask total (green)
    \item Flask mixed (orange)
\end{itemize}

Then, the table of each dimuon component, efficiency corrections, and resulting mass bin centroids are presented for each $x_F$ bin. Similarly, for LD$_2$ cross-section calculations, we plot the following components:
\begin{itemize}
    \item LD$_2$ total (red)
    \item LD$_2$ mixed (blue)
    \item Flask total (green)
    \item Flask mixed (orange)
\end{itemize}

These represent the total, mixed background, and flask components used in these calculations. Alongside the distributions are the precise kinematic values, raw counts, and resulting mass bin centroids for each $x_F$ bin. Finally, the table of each dimuon component, efficiency corrections, and resulting mass bin centroids are presented for each $x_F$ bin for the LD$_2$ target as well. This detailed breakdown allows for a clear understanding of how the final cross-section points are derived from the raw data and corrections.
"""
        tex_file.write(intro_text)

        def fnum(val_str, fmt=".4f"):
            try: return f"{float(val_str):{fmt}}"
            except ValueError: return "-"
                
        def fint(val_str):
            try: return f"{int(float(val_str))}"
            except ValueError: return "-"

        for i_x in range(len(xf_bins_np) - 1):
            x_low, x_high = xf_bins_np[i_x], xf_bins_np[i_x+1]
            
            csv_lh2 = f"Table_Kinematics_LH2_xF_{x_low:.2f}_{x_high:.2f}.csv"
            csv_ld2 = f"Table_Kinematics_LD2_xF_{x_low:.2f}_{x_high:.2f}.csv"
            plot_lh2 = f"{IMAGE_PATH_PREFIX}MassDist_LH2_xF_{x_low:.2f}_{x_high:.2f}.pdf"
            plot_ld2 = f"{IMAGE_PATH_PREFIX}MassDist_LD2_xF_{x_low:.2f}_{x_high:.2f}.pdf"
            
            sec_title = f"\\texorpdfstring{{Mass bin centroid calculations in $x_F \\in [{x_low:.2f}, {x_high:.2f})$}}{{Mass bin centroid calculations in x\_F in [{x_low:.2f}, {x_high:.2f})]}}"
            tex_file.write(f"\n\\clearpage\n\\subsection{{{sec_title}}}\n")
            
            lh2_figure_block = f"""
\\begin{{figure}}[htbp]
    \\centering
    \\includegraphics[width=0.7\\textwidth]{{{plot_lh2}}}
    \\caption{{Invariant Mass distributions for LH$_2$ in the range $x_F \\in [{x_low:.2f}, {x_high:.2f})$. Components shown include target total (red), target mix (blue), flask total (green), and flask mix (orange). The values listed above each bin represent the efficiency-weighted corrected mass $\\sum M / \\epsilon$ and the raw counts $(N)$.}}
    \\label{{fig:mass_dist_lh2_xf_{i_x}}}
\\end{{figure}}
"""
            tex_file.write(lh2_figure_block)

            if os.path.exists(csv_lh2):
                tex_file.write("\\begin{table}[htbp]\n    \\centering\n")
                tex_file.write(f"    \\caption{{LH$_2$ raw counts, efficiency corrections, yields, and final data-driven centroids for $x_F \\in [{x_low:.2f}, {x_high:.2f})$.}}\n")
                tex_file.write(f"    \\label{{tab:kinematics_lh2_xf_{i_x}}}\n")
                tex_file.write("    \\resizebox{\\textwidth}{!}{\n    \\begin{tabular}{cc|c|cccc|ccc|cc}\n        \\hline\\hline\n")
                tex_file.write("        Mass & Geo. &  $\\langle M \\rangle$ & $N_{\\text{tot}}^{\\text{LH}_2}$ & $N_{\\text{mix}}^{\\text{LH}_2}$ & $N_{\\text{tot}}^{\\text{fl}}$ & $N_{\\text{mix}}^{\\text{fl}}$ & $\\langle\\epsilon^{\\text{LH}_2}_{\\text{sig.}}\\rangle$ & $\\langle\\epsilon^{\\text{LD}_2}_{\\text{sig.}}\\rangle$ & $\\langle\\epsilon^{\\text{fl}}_{\\text{sig.}}\\rangle$ & $\\sum M_{\\text{corr}}^{\\text{LH}_2}$ & $Y_{\\text{corr}}^{\\text{LH}_2}$ \\\\\n")
                tex_file.write("        Bin (GeV) & $\\langle M \\rangle$ & (LH$_2$) & & & & & & & & (Num.) & (Denom.) \\\\\n        \\hline\n")
                
                with open(csv_lh2, 'r') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        m_bin = row["Mass Bin"]
                        m_geo = fnum(row["Mass Center"], ".2f")
                        m_lh2 = fnum(row["LH2 Mass Bin Average"], ".3f")
                        n_lh2_t = fint(row["N_LH2_total"])
                        n_lh2_m = fint(row["N_LH2_mixed"])
                        n_fl_t  = fint(row["N_flask_total"])
                        n_fl_m  = fint(row["N_flask_mixed"])
                        e_lh2 = fnum(row["eps_LH2"], ".4f")
                        e_ld2 = fnum(row["eps_LD2"], ".4f")
                        e_fl  = fnum(row["eps_flask"], ".4f")
                        num_lh2 = fnum(row["Corrected Total Mass LH2 (num)"], ".4f")
                        den_lh2 = fnum(row["Corrected Yield LH2 (denom)"], ".4f")
                        
                        tex_file.write(f"        {{{m_bin}}} & {m_geo} & {m_lh2} & {n_lh2_t} & {n_lh2_m} & {n_fl_t} & {n_fl_m} & {e_lh2} & {e_ld2} & {e_fl} & {num_lh2} & {den_lh2} \\\\\n")
                
                tex_file.write("        \\hline\\hline\n    \\end{tabular}\n    }\n\\end{table}\n\n")
            else:
                tex_file.write(f"\n\\textit{{Table data missing: {csv_lh2} not found.}}\n\n")

            ld2_figure_block = f"""
\\begin{{figure}}[htbp]
    \\centering
    \\includegraphics[width=0.7\\textwidth]{{{plot_ld2}}}
    \\caption{{Invariant Mass distributions for LD$_2$ in the range $x_F \\in [{x_low:.2f}, {x_high:.2f})$. Components shown include target total (red), target mix (blue), flask total (green), and flask mix (orange). The values listed above each bin represent the efficiency-weighted corrected mass $\\sum M / \\epsilon$ and the raw counts $(N)$.}}
    \\label{{fig:mass_dist_ld2_xf_{i_x}}}
\\end{{figure}}
"""
            tex_file.write(ld2_figure_block)

            if os.path.exists(csv_ld2):
                tex_file.write("\\begin{table}[htbp]\n    \\centering\n")
                tex_file.write(f"    \\caption{{LD$_2$ raw counts, efficiency corrections, yields, and final data-driven centroids for $x_F \in [{x_low:.2f}, {x_high:.2f})$.}}\n")
                tex_file.write(f"    \\label{{tab:kinematics_ld2_xf_{i_x}}}\n")
                tex_file.write("    \\resizebox{\\textwidth}{!}{\n    \\begin{tabular}{cc|c|cccccc|ccc|cc}\n        \\hline\\hline\n")
                tex_file.write("        Mass & Geo. &  $\\langle M \\rangle$ & $N_{\\rm tot}^{\\rm LD2}$ & $N_{\\rm mix}^{\\rm LD2}$ & $N_{\\rm tot}^{\\rm LH2}$ & $N_{\\rm mix}^{\\rm LH2}$ & $N_{\\rm tot}^{\\rm fl}$ & $N_{\\rm mix}^{\\rm fl}$ & $\\langle\\epsilon^{\\rm LH2}_{\\rm sig.}\\rangle$ & $\\langle\\epsilon^{\\rm LD2}_{\\rm sig.}\\rangle$ & $\\langle\\epsilon^{\\rm fl}_{\\rm sig.}\\rangle$ & $\\sum M_{\\rm corr}^{\\rm LD2}$ & $Y_{\\rm corr}^{\\rm LD2}$ \\\\\n")
                tex_file.write("        Mass Bin (GeV) & $\\langle M \\rangle$ & (LD2) & & & & & & & & & & (Num.) & (Denom.) \\\\\n        \\hline\n")
                
                with open(csv_ld2, 'r') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        m_bin = row["Mass Bin"]
                        m_geo = fnum(row["Mass Center"], ".2f")
                        m_ld2 = fnum(row["LD2 Mass Bin Average"], ".3f")
                        n_ld2_t = fint(row["N_LD2_total"])
                        n_ld2_m = fint(row["N_LD2_mixed"])
                        n_lh2_t = fint(row["N_LH2_total"])
                        n_lh2_m = fint(row["N_LH2_mixed"])
                        n_fl_t  = fint(row["N_flask_total"])
                        n_fl_m  = fint(row["N_flask_mixed"])
                        e_lh2 = fnum(row["eps_LH2"], ".4f")
                        e_ld2 = fnum(row["eps_LD2"], ".4f")
                        e_fl  = fnum(row["eps_flask"], ".4f")
                        num_ld2 = fnum(row["Corrected Total Mass LD2 (num)"], ".4f")
                        den_ld2 = fnum(row["Corrected Yield LD2 (denom)"], ".4f")
                        
                        tex_file.write(f"        {{{m_bin}}} & {m_geo} & {m_ld2} & {n_ld2_t} & {n_ld2_m} & {n_lh2_t} & {n_lh2_m} & {n_fl_t} & {n_fl_m} & {e_lh2} & {e_ld2} & {e_fl} & {num_ld2} & {den_ld2} \\\\\n")
                
                tex_file.write("        \\hline\\hline\n    \\end{tabular}\n    }\n\\end{table}\n")
            else:
                tex_file.write(f"\n\\textit{{Table data missing: {csv_ld2} not found.}}\n")

    print(f"LaTeX appendix successfully generated: {latex_filename}")


# ==========================================
# Main Execution
# ==========================================
xF_bin_ranges = [
    (0.10, 0.15), (0.15, 0.20), (0.20, 0.25), (0.25, 0.30),
    (0.30, 0.35), (0.35, 0.40), (0.40, 0.45), (0.45, 0.50),
    (0.50, 0.55), (0.55, 0.60), (0.60, 0.65), (0.65, 0.70),
    (0.70, 0.75), (0.75, 0.80), (0.80, 0.85), (0.85, 0.90)
]

def main():
    output_filename = "Processed_Kinematic_Hists.root"
    print(f"Creating ROOT output file: {output_filename}")
    out_file = ROOT.TFile(output_filename, "RECREATE")

    hists_lh2, hists_ld2, hists_flask = process_all_datasets(
        "../../../HodoEfficiency/RS67/LH2/merged_RS67_3089_LH2_recoeff_hodoeff.root",
        "../../../HodoEfficiency/RS67/LD2/merged_RS67_3089_LD2_recoeff_hodoeff.root",
        "../../../HodoEfficiency/RS67/EmptyFlask/merged_RS67_3089_Flask_recoeff_hodoeff.root",
        out_file
    )

    if hists_lh2 and hists_flask:
        h_sub_dict_lh2 = generate_subtracted_plot(hists_lh2, hists_flask, FLASK_NORM_LH2, "LH2", out_file)
        if h_sub_dict_lh2:
            calculate_and_plot_cross_section(h_sub_dict_lh2, "LH2", GLOBAL_CONSTANT_LH2, out_file)

    if hists_ld2 and hists_lh2 and hists_flask:
        h_sub_dict_pd = generate_pd_subtracted_plot(hists_ld2, hists_lh2, hists_flask, out_file)
        if h_sub_dict_pd:
            calculate_and_plot_cross_section(h_sub_dict_pd, "LD2", GLOBAL_CONSTANT_LD2, out_file)

    # Generate the overlay comparison plots for LD2 at the geometric bin center 
    generate_overlay_plot(out_file, "LD2", "GeoCenter")

    out_file.Write()
    out_file.Close()
    print("\nAll histograms, tables, cross-section plots, and overlays generated successfully.")

    # Generate the LaTeX Appendix automatically at the end
    generate_latex_appendix()

if __name__ == "__main__":
    main()