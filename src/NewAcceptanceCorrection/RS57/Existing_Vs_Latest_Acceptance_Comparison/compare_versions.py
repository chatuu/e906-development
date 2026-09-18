import os
import ROOT
import numpy as np

# ==============================================================================
# Helper Functions
# ==============================================================================
def get_hist_max_with_error(h):
    """Calculates the absolute maximum of a TH1 including error bars."""
    max_val = 0.0
    for i in range(1, h.GetNbinsX() + 1):
        val = h.GetBinContent(i) + h.GetBinError(i)
        if val > max_val:
            max_val = val
    return max_val

def get_hist_min_with_error(h):
    """Calculates the absolute minimum of a TH1 including error bars, ignoring 0 bins."""
    min_val = 1e9
    valid = False
    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(i) > 0: 
            val = h.GetBinContent(i) - h.GetBinError(i)
            if val < min_val:
                min_val = val
            valid = True
    return min_val if valid else 0.0

def get_optimal_legend_position(hists, y_min=0.75, y_max=0.88):
    """
    Determines optimal TLegend X-coordinates to avoid masking data.
    """
    if not isinstance(hists, list):
        hists = [hists]

    left_max = -np.inf
    right_max = -np.inf

    for h in hists:
        n_bins = h.GetNbinsX()
        mid_bin = n_bins // 2
        
        def get_side_max(start_bin, end_bin):
            side_max = -np.inf
            for i in range(start_bin, end_bin + 1):
                c = h.GetBinContent(i)
                e = h.GetBinError(i)
                if (c == 0 and e == 0) or (e >= 1.5 * c and c > 0):
                    continue
                val = c + e
                if val > side_max:
                    side_max = val
            return side_max
        
        curr_left_max = get_side_max(1, mid_bin)
        curr_right_max = get_side_max(mid_bin + 1, n_bins)
        
        if curr_left_max > left_max: left_max = curr_left_max
        if curr_right_max > right_max: right_max = curr_right_max

    if left_max == -np.inf and right_max == -np.inf:
        return (0.65, y_min, 0.88, y_max)

    if left_max > right_max:
        return (0.65, y_min, 0.88, y_max) # Top Right
    else:
        return (0.15, y_min, 0.38, y_max) # Top Left

def add_fit_and_band(h_ratio, pad):
    pad.cd()
    h_ratio.SetStats(0) 
    
    # Fit with a constant line
    h_ratio.Fit("pol0", "QS")
    pad.Update()
    
    fit_func = h_ratio.GetFunction("pol0")
    if fit_func:
        fit_func.SetLineColor(ROOT.kRed)
        fit_func.SetLineWidth(2)
        p0 = fit_func.GetParameter(0)
        p0_err = fit_func.GetParError(0)
        x_min = h_ratio.GetXaxis().GetXmin()
        x_max = h_ratio.GetXaxis().GetXmax()
        
        # Error band
        error_box = ROOT.TBox(x_min, p0 - p0_err, x_max, p0 + p0_err)
        error_box.SetFillColorAlpha(ROOT.kRed, 0.3) 
        error_box.Draw("SAME")
        pad._error_band = error_box 
        
        fit_func.Draw("SAME")
        h_ratio.Draw("E1 SAME")
        
        # TLatex label
        x_pos = x_min + (x_max - x_min) * 0.02 
        y_range = h_ratio.GetMaximum() - h_ratio.GetMinimum()
        y_pos = p0 + (y_range * 0.05) 
        
        latex = ROOT.TLatex()
        label_size = h_ratio.GetYaxis().GetLabelSize()
        if label_size == 0: label_size = 0.05 
        latex.SetTextSize(label_size * 0.9)
        latex.SetTextColor(ROOT.kRed)
        latex.SetTextAlign(11) 
        latex.DrawLatex(x_pos, y_pos, f"Fit: {p0:.3f} #pm {p0_err:.3f}")
        
        # Dashed line at ratio = 1.0
        line = ROOT.TLine(x_min, 1.0, x_max, 1.0)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw("SAME")
        pad._fit_line = line

# ==============================================================================
# Main Plotting Logic
# ==============================================================================
def create_comparison_plot(file_old, file_new, plot_base_name, title, x_title, target):
    """
    Extracts the old and new histograms for a specific target, divides them, 
    and draws the standard split canvas.
    """
    hist_path = f"Acceptance_Ratios_SplitCanvas/h_{target.lower()}_acc_{plot_base_name}"
    
    h_old = file_old.Get(hist_path)
    h_new = file_new.Get(hist_path)
    
    if not h_old or not h_new:
        print(f"Warning: Could not find {hist_path} in one or both files. Skipping.")
        return

    # Clone to detach from file memory
    h_old = h_old.Clone(f"old_{target}_{plot_base_name}")
    h_new = h_new.Clone(f"new_{target}_{plot_base_name}")
    
    # Calculate Ratio (Latest / Existing). ROOT automatically propagates independent errors here.
    h_ratio = h_new.Clone(f"ratio_{target}_{plot_base_name}")
    h_ratio.Divide(h_old)

    # Styling
    h_old.SetLineColor(ROOT.kBlue)
    h_old.SetMarkerColor(ROOT.kBlue)
    h_old.SetTitle("")
    
    h_new.SetLineColor(ROOT.kRed)
    h_new.SetMarkerColor(ROOT.kRed)
    h_new.SetTitle("")
    
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetMarkerColor(ROOT.kBlack)
    h_ratio.SetTitle("")

    # Create Canvas and Pads
    c_name = f"Compare_{target}_{plot_base_name}"
    c = ROOT.TCanvas(c_name, title, 800, 800)
    
    pad1 = ROOT.TPad(f"pad1_{c_name}", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.15)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.Draw()
    
    pad2 = ROOT.TPad(f"pad2_{c_name}", "pad2", 0, 0.0, 1, 0.35)
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.3)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.Draw()

    # --- TOP PAD (Distributions) ---
    pad1.cd()
    max_y = max(get_hist_max_with_error(h_old), get_hist_max_with_error(h_new))
    h_old.SetMaximum(max_y * 1.35) 
    h_old.SetMinimum(0)
    
    # Format Top Axes
    h_old.GetXaxis().SetLabelSize(0.04)
    h_old.GetXaxis().SetTitleSize(0.045)
    h_old.GetYaxis().SetLabelSize(0.04)
    h_old.GetYaxis().SetTitleSize(0.045)
    h_old.GetYaxis().SetTitleOffset(1.2)
    h_old.GetYaxis().SetTitle(f"{target} Acceptance")
    h_old.GetXaxis().SetTitle(x_title)
    
    h_old.Draw("E1")
    h_new.Draw("E1 SAME")
    
    # Label the overall title manually inside the pad
    latex_title = ROOT.TLatex()
    latex_title.SetNDC()
    latex_title.SetTextSize(0.045)
    latex_title.SetTextAlign(21)
    latex_title.DrawLatex(0.5, 0.92, f"{target} Existing vs Latest: {title}")

    # Legend
    leg_coords = get_optimal_legend_position([h_old, h_new], y_min=0.75, y_max=0.88)
    leg = ROOT.TLegend(*leg_coords)
    leg.SetBorderSize(0)
    leg.SetFillColor(ROOT.kWhite)
    leg.AddEntry(h_old, "Existing Acceptance", "lep")
    leg.AddEntry(h_new, "Latest Acceptance", "lep")
    leg.Draw()

    # --- BOTTOM PAD (Ratio) ---
    pad2.cd()
    
    # Intelligent Scaling for Ratio
    valid_bins = []
    for i in range(1, h_ratio.GetNbinsX() + 1):
        if h_ratio.GetBinContent(i) > 0 and h_ratio.GetBinError(i) < h_ratio.GetBinContent(i):
            valid_bins.append(h_ratio.GetBinContent(i) + h_ratio.GetBinError(i))
            valid_bins.append(h_ratio.GetBinContent(i) - h_ratio.GetBinError(i))
            
    if len(valid_bins) > 0:
        max_ratio_val = np.max(valid_bins)
        min_ratio_val = np.min(valid_bins)
        range_padding = (max_ratio_val - min_ratio_val) * 0.3
        if range_padding == 0: range_padding = 0.2
        h_ratio.SetMaximum(max_ratio_val + range_padding)
        h_ratio.SetMinimum(max(0.0, min_ratio_val - range_padding))
    else:
        h_ratio.SetMaximum(2.0)
        h_ratio.SetMinimum(0.0)
        
    # Format Bottom Axes
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetLabelSize(0.08)
    h_ratio.GetYaxis().SetTitleSize(0.08)
    h_ratio.GetYaxis().SetTitleOffset(0.6)
    h_ratio.GetYaxis().SetTitle("Latest / Existing")
    
    h_ratio.GetXaxis().SetLabelSize(0.1)
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetTitleOffset(1.0)
    h_ratio.GetXaxis().SetTitle(x_title)
    
    h_ratio.Draw("E1")
    add_fit_and_band(h_ratio, pad2)
    
    # Save the canvas
    output_filename = f"Compare_Latest_Vs_Existing_{target}_{plot_base_name}.pdf"
    c.SaveAs(output_filename)
    print(f"Generated {output_filename}")


def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gROOT.SetBatch(True) 

    # Expand user paths
    old_file_path = os.path.expanduser("~/github/e906-development/src/AcceptanceCorrection/acceptance_no_unfolding_bins/acceptance_mass_xF_unfolding.root")
    new_file_path = os.path.expanduser("~/github/e906-development/src/NewAcceptanceCorrection/RS57/acceptance_mass_xF_unfolding.root")

    if not os.path.exists(old_file_path):
        print(f"Error: Old file not found at {old_file_path}")
        return
    if not os.path.exists(new_file_path):
        print(f"Error: New file not found at {new_file_path}")
        return

    print("Opening ROOT files...")
    file_old = ROOT.TFile.Open(old_file_path, "READ")
    file_new = ROOT.TFile.Open(new_file_path, "READ")
    
    targets = ["LH2", "LD2"]

    # 1. Fully Integrated Plots
    print("\n--- Generating Fully Integrated Comparisons ---")
    for t in targets:
        create_comparison_plot(file_old, file_new, "Acceptance_Mass_All_xF_pT", "Integrated Mass", "Mass [GeV]", t)
        create_comparison_plot(file_old, file_new, "Acceptance_xF_All_Mass_pT", "Integrated x_{F}", "x_{F}", t)
        create_comparison_plot(file_old, file_new, "Acceptance_pT_All_Mass_xF", "Integrated p_{T}", "p_{T} [GeV/c]", t)

    # 2. Sliced Binned Plots (Mass in xF bins 0 to 15)
    print("\n--- Generating Sliced xF Bin Comparisons ---")
    for t in targets:
        for i in range(16):
            plot_base = f"Acceptance_Mass_All_pT_xF_bin{i}"
            title = f"Mass in x_{{F}} Bin {i}"
            create_comparison_plot(file_old, file_new, plot_base, title, "Mass [GeV]", t)

    file_old.Close()
    file_new.Close()
    print("\nAll comparisons completed successfully.")

if __name__ == "__main__":
    main()