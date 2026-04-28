import ROOT
import os
import sys
import numpy as np

def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # ==========================================
    # 1. Setup Bins and Open Files
    # ==========================================
    # Binning matching your master script
    mass_bins_np = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
    xf_bins_np = [round(i * 0.05, 2) for i in range(17)]

    kin_file_path = "Processed_Kinematic_Hists.root"
    if not os.path.exists(kin_file_path):
        print(f"Error: {kin_file_path} not found. Run your master script first.")
        sys.exit(1)
        
    f_kin = ROOT.TFile.Open(kin_file_path)

    # Theoretical Files for LH2 (Proton target)
    theory_ct18_path = "CT18_xFnew_p_1sigma.root"
    theory_nnpdf_path = "NNPDF40_xFnew_p.root"
    
    f_ct18 = ROOT.TFile.Open(theory_ct18_path) if os.path.exists(theory_ct18_path) else None
    f_nnpdf = ROOT.TFile.Open(theory_nnpdf_path) if os.path.exists(theory_nnpdf_path) else None

    # ==========================================
    # 2. Loop over xF Bins
    # ==========================================
    for i_x in range(16):
        xf_min, xf_max = xf_bins_np[i_x], xf_bins_np[i_x+1]
        
        # --- A. Load "Current" TGraphs (Red) ---
        g_new_stat_orig = f_kin.Get(f"g_xsec_LH2_{i_x}")
        g_new_sys_orig  = f_kin.Get(f"g_sys_LH2_{i_x}")
        
        if not g_new_stat_orig or not g_new_sys_orig:
            print(f"Warning: New graphs for xF bin {i_x} not found. Skipping.")
            continue

        g_new_stat = ROOT.TGraphErrors()
        g_new_sys = ROOT.TGraphErrors()
        
        y_max_data = 1e-3
        y_min_data = 10.0
        
        for i in range(g_new_stat_orig.GetN()):
            x_val = g_new_stat_orig.GetX()[i]
            y_val = g_new_stat_orig.GetY()[i]
            ey_stat = g_new_stat_orig.GetErrorY(i)
            ex_sys = g_new_sys_orig.GetErrorX(i)
            ey_sys = g_new_sys_orig.GetErrorY(i)
            
            if y_val > 0:
                g_new_stat.SetPoint(i, x_val, y_val)
                g_new_stat.SetPointError(i, 0.0, ey_stat)
                g_new_sys.SetPoint(i, x_val, y_val)
                g_new_sys.SetPointError(i, ex_sys, ey_sys)
                
                y_max_data = max(y_max_data, float(y_val) + max(ey_stat, ey_sys))
                y_min_data = min(y_min_data, float(y_val) - max(ey_stat, ey_sys))

        # --- B. Load "Previous" TGraphs (Blue) ---
        old_filename = f"LH2_{i_x}_updatedPoT.root"
        if not os.path.exists(old_filename):
            old_filename = f"result_rootfiles/LH2_{i_x}_updatedPoT.root"
            
        g_old_stat = ROOT.TGraphErrors()
        g_old_sys = ROOT.TGraphErrors()
        
        if os.path.exists(old_filename):
            f_old = ROOT.TFile.Open(old_filename)
            g_orig_stat = f_old.Get("gAccCor_Stat_avgMass")
            g_orig_sys = f_old.Get("gAccCor_Syst_avgMass")
            
            if g_orig_stat and g_orig_sys:
                for i in range(g_orig_stat.GetN()):
                    x_geom = g_orig_sys.GetX()[i]
                    y = g_orig_stat.GetY()[i]
                    ey_stat = g_orig_stat.GetErrorY(i)
                    ey_sys = g_orig_sys.GetErrorY(i)
                    
                    # Fix: Dynamic width calculation for previous results
                    ex_geom = 0.15 
                    for b_idx in range(len(mass_bins_np)-1):
                        if mass_bins_np[b_idx] <= x_geom <= mass_bins_np[b_idx+1]:
                            ex_geom = (mass_bins_np[b_idx+1] - mass_bins_np[b_idx]) / 2.0
                            break
                    
                    if y > 0:
                        g_old_stat.SetPoint(i, x_geom, y)
                        g_old_stat.SetPointError(i, 0.0, ey_stat)
                        g_old_sys.SetPoint(i, x_geom, y)
                        g_old_sys.SetPointError(i, ex_geom, ey_sys)
                        
                        y_max_data = max(y_max_data, float(y) + max(ey_stat, ey_sys))
                        y_min_data = min(y_min_data, float(y) - max(ey_stat, ey_sys))
            f_old.Close()

        # --- C. Visual Formatting & Canvas ---
        c = ROOT.TCanvas(f"c_comp_{i_x}", f"Comparison xF {i_x}", 800, 600)
        c.SetLogy()
        c.SetTickx(1)
        c.SetTicky(1)
        
        mg = ROOT.TMultiGraph()
        mg.SetTitle(f";Mass [GeV];M^{{3}} d^{{2}}\\sigma / dM dx_{{F}} [nb GeV^{{2}}]")
        
        leg = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
        leg.SetBorderSize(0)

        # 1. Theoretical Curves
        gr_name = f"gr_xFbin{i_x}"
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
                g_nnpdf.SetLineColor(ROOT.kGray + 2)
                g_nnpdf.SetFillColorAlpha(ROOT.kGray + 1, 0.5)
                g_nnpdf.SetFillStyle(3002)
                mg.Add(g_nnpdf, "L3")
                leg.AddEntry(g_nnpdf, "NNPDF4.0 NLO", "lf")

        # 2. Previous Results (Blue)
        if g_old_stat.GetN() > 0:
            g_old_sys.SetMarkerSize(0)
            g_old_sys.SetLineColor(ROOT.kBlue)
            g_old_sys.SetFillColorAlpha(ROOT.kAzure - 9, 0.5)
            g_old_sys.SetFillStyle(1001)
            mg.Add(g_old_sys, "2")
            g_old_stat.SetMarkerStyle(21) 
            g_old_stat.SetMarkerColor(ROOT.kBlue)
            g_old_stat.SetLineColor(ROOT.kBlue)
            mg.Add(g_old_stat, "P")
            leg.AddEntry(g_old_stat, "Previous Results (LH2)", "lep")
            leg.AddEntry(g_old_sys, "Previous Syst.", "f")

        # 3. Current Results (Red)
        if g_new_stat.GetN() > 0:
            g_new_sys.SetMarkerSize(0)
            g_new_sys.SetLineColor(ROOT.kRed)
            g_new_sys.SetFillColorAlpha(ROOT.kPink - 9, 0.5)
            g_new_sys.SetFillStyle(1001)
            mg.Add(g_new_sys, "2")
            g_new_stat.SetMarkerStyle(20) 
            g_new_stat.SetMarkerColor(ROOT.kRed)
            g_new_stat.SetLineColor(ROOT.kRed)
            mg.Add(g_new_stat, "P")
            leg.AddEntry(g_new_stat, "Current Data (LH2)", "lep")
            leg.AddEntry(g_new_sys, "Current Syst.", "f")

        # Determine plot limits
        if y_min_data > 0 and y_min_data < y_max_data:
            y_min_plot = y_min_data * 0.2
            y_max_plot = y_max_data * 5.0
        else:
            y_min_plot = 1e-3
            y_max_plot = 3.0

        mg.SetMinimum(y_min_plot)
        mg.SetMaximum(y_max_plot)
        mg.Draw("A")
        mg.GetXaxis().CenterTitle()
        mg.GetYaxis().CenterTitle()
        mg.GetXaxis().SetLimits(4.2, 8.7)
        leg.Draw()

        # --- D. Annotations ---
        internal_title_text = f"pp cross-section within {xf_min:.2f} #leq x_{{F}} < {xf_max:.2f}"
        internal_title = ROOT.TLatex()
        internal_title.SetNDC(True)
        internal_title.SetTextFont(42)
        internal_title.SetTextSize(0.04)
        internal_title.SetTextAlign(13)
        internal_title.DrawLatex(0.14, 0.86, internal_title_text)

        lumi_note = ROOT.TLatex()
        lumi_note.SetNDC(False)
        lumi_note.SetTextColor(ROOT.kBlack)
        lumi_note.SetTextAlign(13)
        lumi_note.SetTextSize(0.025)
        # Use calculated plot minimum instead of GetMinimum()
        lumi_note.DrawLatex(4.3, y_min_plot * 1.5, "10% global uncertainty due to the integrated luminosity is not included in the error bands")

        # Save Output
        out_pdf = f"Comparison_LH2_xF_{xf_min:.2f}_{xf_max:.2f}.pdf"
        c.SaveAs(out_pdf)
        c.Close()
        print(f"Saved comparison: {out_pdf}")

    f_kin.Close()
    if f_ct18: f_ct18.Close()
    if f_nnpdf: f_nnpdf.Close()

if __name__ == "__main__":
    main()