#!/usr/bin/env python3
import os
import csv
import math
import numpy as np
import ROOT

def get_val_errs(g_xsec, g_sys, target_x, tol=0.01):
    if not g_xsec or not g_sys:
        return None, None, None
    for i in range(g_xsec.GetN()):
        if abs(g_xsec.GetX()[i] - target_x) < tol:
            return g_xsec.GetY()[i], g_xsec.GetErrorY(i), g_sys.GetErrorY(i)
    return None, None, None

def format_cell(val, stat, sys):
    if val is None:
        return "-"
    return f"{val:.6f} \\pm {stat:.6f} \\pm {sys:.6f}"

def main():
    # Setup ROOT to run in batch mode and format stats box
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning  # Silences Info messages like TCanvas::Print
    ROOT.gStyle.SetOptStat(1110)
    ROOT.gStyle.SetStatX(0.88); ROOT.gStyle.SetStatY(0.88)
    
    FILE_ROOT = "All_XSec_Objects.root"
    CSV_OUTPUT = "Unfolding_Systematics_Bootstrap.csv"
    PLOT_TYPE = "GeoCenter"
    HIST_DIR = "Systematics_Histograms"

    os.makedirs(HIST_DIR, exist_ok=True)

    XF_BINS = np.round(np.arange(-0.05, 0.90, 0.05), 2)
    MASS_BINS = np.array([3.9, 4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.8, 10.0], dtype=float)

    if not os.path.exists(FILE_ROOT):
        return

    f = ROOT.TFile.Open(FILE_ROOT, "READ")
    csv_data = []

    headers = [
        r"Target", r"xF bin", r"Mass bin",
        r"Raw Cross Section \pm stat. error \pm syst. error",
        r"Unfolded (clean) Cross Section (iter=3) \pm stat. error \pm syst. error",
        r"Unfolded (messy) Cross Section (iter=3) \pm stat. error \pm syst. error",
        r"Unfolded (clean) Cross Section (iter=4) \pm stat. error \pm syst. error",
        r"Unfolded (messy) Cross Section (iter=4) \pm stat. error \pm syst. error",
        r"Unfolding syst (clean) [%]", r"Unfolding syst (messy) [%]"
    ]

    f_out = ROOT.TFile("Unfolding_Sys_Hists.root", "RECREATE")
    canvas = ROOT.TCanvas("c_sys", "Systematics", 800, 600)
    canvas.SetLeftMargin(0.12); canvas.SetBottomMargin(0.12)

    for target in ["LH2", "LD2"]:
        for i_x in range(len(XF_BINS) - 1):
            xf_min, xf_max = XF_BINS[i_x], XF_BINS[i_x + 1]
            xf_str = f"{xf_min:.2f} <= xF < {xf_max:.2f}"
            dir_path = f"CrossSections_{target}"
            
            g_raw_xsec = f.Get(f"{dir_path}/g_xsec_{target}_{i_x}_Raw_{PLOT_TYPE}")
            g_raw_sys  = f.Get(f"{dir_path}/g_sys_{target}_{i_x}_Raw_{PLOT_TYPE}")
            g_c3_xsec = f.Get(f"{dir_path}/g_xsec_{target}_{i_x}_UnfClean_{PLOT_TYPE}")
            g_c3_sys  = f.Get(f"{dir_path}/g_sys_{target}_{i_x}_UnfClean_{PLOT_TYPE}")
            g_m3_xsec = f.Get(f"{dir_path}/g_xsec_{target}_{i_x}_UnfMessy_{PLOT_TYPE}")
            g_m3_sys  = f.Get(f"{dir_path}/g_sys_{target}_{i_x}_UnfMessy_{PLOT_TYPE}")

            h2_c3 = f.Get(f"Unfolding_Toys/h2_toy_yields_{target}_xF_{i_x}_Clean_iter3")
            h2_c4 = f.Get(f"Unfolding_Toys/h2_toy_yields_{target}_xF_{i_x}_Clean_iter4")
            h2_m3 = f.Get(f"Unfolding_Toys/h2_toy_yields_{target}_xF_{i_x}_Messy_iter3")
            h2_m4 = f.Get(f"Unfolding_Toys/h2_toy_yields_{target}_xF_{i_x}_Messy_iter4")

            if not g_raw_xsec or g_raw_xsec.GetN() == 0:
                continue

            for pt in range(g_raw_xsec.GetN()):
                mass_center = g_raw_xsec.GetX()[pt]
                
                mass_str = "Unknown"
                m_bin_idx = -1
                for i_m in range(len(MASS_BINS) - 1):
                    m_min, m_max = MASS_BINS[i_m], MASS_BINS[i_m+1]
                    geo_center = (m_min + m_max) / 2.0
                    if abs(geo_center - mass_center) < 0.01:
                        mass_str = f"{m_min:.2f} <= Mass < {m_max:.2f}"
                        m_bin_idx = i_m + 1
                        break
                
                if m_bin_idx == -1: continue

                v_raw, e_raw_stat, e_raw_sys = g_raw_xsec.GetY()[pt], g_raw_xsec.GetErrorY(pt), g_raw_sys.GetErrorY(pt)
                v_c3, e_c3_stat, e_c3_sys = get_val_errs(g_c3_xsec, g_c3_sys, mass_center)
                v_m3, e_m3_stat, e_m3_sys = get_val_errs(g_m3_xsec, g_m3_sys, mass_center)

                sys_clean, sys_messy = "-", "-"
                
                if h2_c3 and h2_c4:
                    h_sys_c = ROOT.TH1D(f"h_sys_clean_{target}_{i_x}_{m_bin_idx}", f"Clean Sys: {target} xF [{xf_min:.2f}, {xf_max:.2f}) Mass [{m_min:.2f}, {m_max:.2f})", 50, 0, 25)
                    n_toys = h2_c3.GetNbinsY()
                    for toy in range(1, n_toys + 1):
                        y3 = h2_c3.GetBinContent(m_bin_idx, toy)
                        y4 = h2_c4.GetBinContent(m_bin_idx, toy)
                        if y3 > 0: h_sys_c.Fill(abs(y4 - y3)/y3 * 100.0)
                    
                    if h_sys_c.GetEntries() > 0:
                        sys_clean = f"{h_sys_c.GetMean():.2f}"
                        f_out.cd(); h_sys_c.Write()
                        
                        h_sys_c.SetFillColor(ROOT.kAzure+1)
                        h_sys_c.GetXaxis().SetTitle("Relative Difference [%]")
                        h_sys_c.GetYaxis().SetTitle("Pseudo-experiments (Toys)")
                        h_sys_c.GetXaxis().SetTitleSize(0.045); h_sys_c.GetYaxis().SetTitleSize(0.045)
                        h_sys_c.Draw("HIST")
                        canvas.SaveAs(f"{HIST_DIR}/SysHist_Clean_{target}_xF_{i_x}_Mass_{m_bin_idx}.pdf")

                if h2_m3 and h2_m4:
                    h_sys_m = ROOT.TH1D(f"h_sys_messy_{target}_{i_x}_{m_bin_idx}", f"Messy Sys: {target} xF [{xf_min:.2f}, {xf_max:.2f}) Mass [{m_min:.2f}, {m_max:.2f})", 50, 0, 25)
                    n_toys = h2_m3.GetNbinsY()
                    for toy in range(1, n_toys + 1):
                        y3 = h2_m3.GetBinContent(m_bin_idx, toy)
                        y4 = h2_m4.GetBinContent(m_bin_idx, toy)
                        if y3 > 0: h_sys_m.Fill(abs(y4 - y3)/y3 * 100.0)
                        
                    if h_sys_m.GetEntries() > 0:
                        sys_messy = f"{h_sys_m.GetMean():.2f}"
                        f_out.cd(); h_sys_m.Write()
                        
                        h_sys_m.SetFillColor(ROOT.kRed-4)
                        h_sys_m.GetXaxis().SetTitle("Relative Difference [%]")
                        h_sys_m.GetYaxis().SetTitle("Pseudo-experiments (Toys)")
                        h_sys_m.GetXaxis().SetTitleSize(0.045); h_sys_m.GetYaxis().SetTitleSize(0.045)
                        h_sys_m.Draw("HIST")
                        canvas.SaveAs(f"{HIST_DIR}/SysHist_Messy_{target}_xF_{i_x}_Mass_{m_bin_idx}.pdf")

                row = {
                    r"Target": target, r"xF bin": xf_str, r"Mass bin": mass_str,
                    r"Raw Cross Section \pm stat. error \pm syst. error": format_cell(v_raw, e_raw_stat, e_raw_sys),
                    r"Unfolded (clean) Cross Section (iter=3) \pm stat. error \pm syst. error": format_cell(v_c3, e_c3_stat, e_c3_sys),
                    r"Unfolded (messy) Cross Section (iter=3) \pm stat. error \pm syst. error": format_cell(v_m3, e_m3_stat, e_m3_sys),
                    r"Unfolded (clean) Cross Section (iter=4) \pm stat. error \pm syst. error": "-",
                    r"Unfolded (messy) Cross Section (iter=4) \pm stat. error \pm syst. error": "-",
                    r"Unfolding syst (clean) [%]": sys_clean, r"Unfolding syst (messy) [%]": sys_messy
                }
                csv_data.append(row)

    f_out.Close()
    f.Close()

    with open(CSV_OUTPUT, "w", newline='') as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=headers)
        writer.writeheader()
        writer.writerows(csv_data)

if __name__ == "__main__":
    main()