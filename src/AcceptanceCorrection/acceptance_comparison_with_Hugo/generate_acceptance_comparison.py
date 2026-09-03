import ROOT
import os

def get_my_hists(tfile, bin_idx):
    """
    Dynamically looks for the LH2 and LD2 mass histograms, checking both 
    potential directory structures and naming conventions.
    """
    possible_dirs = [
        f"mass_sliced_by_xF_bin{bin_idx}",  # No-D1 structure
        f"xF_bin{bin_idx}"                  # With-D1 structure
    ]
    
    possible_hists = [
        f"h_ratio_LH2_mass_sliced_by_xF_bin{bin_idx}",
        f"h_ratio_LH2_mass_sliced_by_xF_bin"
    ]
    
    for d in possible_dirs:
        tdir = tfile.Get(d)
        if tdir:
            for h in possible_hists:
                h_lh2 = tdir.Get(h)
                if h_lh2:
                    h_ld2_name = h.replace("LH2", "LD2")
                    h_ld2 = tdir.Get(h_ld2_name)
                    if h_ld2:
                        h_lh2.SetDirectory(0)
                        h_ld2.SetDirectory(0)
                        return h_lh2, h_ld2
                        
    return None, None

def style_hist(hist, color, marker):
    """Helper to apply consistent styling to histograms"""
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker)
    hist.SetLineWidth(2)

def get_optimal_y_range(hists):
    """
    Scans the actual data points across all histograms to find the optimal Y-axis range,
    minimizing whitespace while leaving room for the legend.
    """
    min_y = float('inf')
    max_y = float('-inf')
    
    for h in hists:
        for i in range(1, h.GetNbinsX() + 1):
            val = h.GetBinContent(i)
            err = h.GetBinError(i)
            if val > 0:  # Ignore completely empty/zero bins
                if (val - err) < min_y:
                    min_y = val - err
                if (val + err) > max_y:
                    max_y = val + err
                    
    # Fallback if all bins were somehow empty
    if min_y == float('inf'):
        return 0.0, 1.0
        
    diff = max_y - min_y
    if diff == 0:
        diff = max_y * 0.1  # Prevent divide-by-zero or flatline issues
        
    # Calculate padding: 10% below (capped at 0) and 45% above for the legend
    pad_bottom = max(0.0, min_y - (diff * 0.1))
    pad_top = max_y + (diff * 0.45)
    
    return pad_bottom, pad_top

def generate_comparison_plots():
    # Keep plots from popping up during generation
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # Define absolute paths
    path_hugo = "/root/github/e906-development/ROOTFiles/Hugo/acceptance_mass_xF_67.root"
    path_nod1 = "/root/github/e906-development/src/AcceptanceCorrection/NoD1Cut/acceptance_mass_xF.root"
    path_withd1 = "/root/github/e906-development/src/AcceptanceCorrection/acceptance_mass_xF_unfolding.root"

    # Open ROOT files
    f_hugo = ROOT.TFile.Open(path_hugo, "READ")
    f_nod1 = ROOT.TFile.Open(path_nod1, "READ")
    f_withd1 = ROOT.TFile.Open(path_withd1, "READ")

    if not f_hugo or not f_nod1 or not f_withd1:
        print("Error: Could not open one or more ROOT files. Check paths.")
        return

    # Output directory
    out_dir = "acceptance_comparisons"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Loop over Hugo's 16 bins (0 to 15)
    for i in range(16):
        hugo_bin = i
        nod1_bin = i       # Starts at 0 (No unfolding)
        withd1_bin = i + 1 # Starts at 1 (Unfolding offset)
        
        print(f"Processing Hugo bin {hugo_bin} / No-D1 bin {nod1_bin} / With-D1 bin {withd1_bin}...")

        # 1. Get Hugo's LH2 and LD2 acceptances
        h_hugo_lh2 = f_hugo.Get(f"h_ratio_LH2_xF_bin{hugo_bin}")
        h_hugo_ld2 = f_hugo.Get(f"h_ratio_LD2_xF_bin{hugo_bin}")
        
        if not h_hugo_lh2 or not h_hugo_ld2:
            print(f"  -> Skipping bin {hugo_bin} for Hugo (histograms not found).")
            continue

        # 2. Get Your acceptances WITHOUT D1 cut
        h_lh2_nod1, h_ld2_nod1 = get_my_hists(f_nod1, nod1_bin)
        if not h_lh2_nod1:
            print(f"  -> ERROR: Could not find No-D1 mass histograms for bin {nod1_bin}.")
            continue

        # 3. Get Your acceptances WITH D1 cut
        h_lh2_withd1, h_ld2_withd1 = get_my_hists(f_withd1, withd1_bin)
        if not h_lh2_withd1:
            print(f"  -> ERROR: Could not find With-D1 mass histograms for bin {withd1_bin}.")
            continue

        # --- Plot LH2 ---
        style_hist(h_hugo_lh2, ROOT.kBlack, 20)
        style_hist(h_lh2_nod1, ROOT.kRed, 21)
        style_hist(h_lh2_withd1, ROOT.kBlue, 22)

        # Dynamic Range & Center Labels
        lh2_min, lh2_max = get_optimal_y_range([h_hugo_lh2, h_lh2_nod1, h_lh2_withd1])
        h_hugo_lh2.SetTitle(f"LH2 Acceptance: xF Bin {hugo_bin}")
        
        h_hugo_lh2.GetYaxis().SetRangeUser(lh2_min, lh2_max)
        h_hugo_lh2.GetYaxis().SetTitle("Acceptance")
        h_hugo_lh2.GetYaxis().CenterTitle(True)
        
        h_hugo_lh2.GetXaxis().SetTitle("Mass (GeV)")
        h_hugo_lh2.GetXaxis().CenterTitle(True)

        c_lh2 = ROOT.TCanvas(f"c_lh2_bin{i}", f"LH2 Acceptance Bin {i}", 800, 600)
        c_lh2.SetGrid()
        h_hugo_lh2.Draw("E1 P")
        h_lh2_nod1.Draw("E1 P SAME")
        h_lh2_withd1.Draw("E1 P SAME")

        leg_lh2 = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        leg_lh2.SetBorderSize(0)
        leg_lh2.SetFillStyle(0)
        leg_lh2.AddEntry(h_hugo_lh2, "Hugo (Existing)", "lp")
        leg_lh2.AddEntry(h_lh2_nod1, "Existing (No D1 Cut)", "lp")
        leg_lh2.AddEntry(h_lh2_withd1, "Existing (With D1 Cut)", "lp")
        leg_lh2.Draw()

        c_lh2.SaveAs(f"{out_dir}/compare_LH2_xFbin{i}.pdf")
        c_lh2.SaveAs(f"{out_dir}/compare_LH2_xFbin{i}.png")
        c_lh2.Close()


        # --- Plot LD2 ---
        style_hist(h_hugo_ld2, ROOT.kBlack, 20)
        style_hist(h_ld2_nod1, ROOT.kRed, 21)
        style_hist(h_ld2_withd1, ROOT.kBlue, 22)

        # Dynamic Range & Center Labels
        ld2_min, ld2_max = get_optimal_y_range([h_hugo_ld2, h_ld2_nod1, h_ld2_withd1])
        h_hugo_ld2.SetTitle(f"LD2 Acceptance: xF Bin {hugo_bin}")
        
        h_hugo_ld2.GetYaxis().SetRangeUser(ld2_min, ld2_max)
        h_hugo_ld2.GetYaxis().SetTitle("Acceptance")
        h_hugo_ld2.GetYaxis().CenterTitle(True)
        
        h_hugo_ld2.GetXaxis().SetTitle("Mass (GeV)")
        h_hugo_ld2.GetXaxis().CenterTitle(True)

        c_ld2 = ROOT.TCanvas(f"c_ld2_bin{i}", f"LD2 Acceptance Bin {i}", 800, 600)
        c_ld2.SetGrid()
        h_hugo_ld2.Draw("E1 P")
        h_ld2_nod1.Draw("E1 P SAME")
        h_ld2_withd1.Draw("E1 P SAME")

        leg_ld2 = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        leg_ld2.SetBorderSize(0)
        leg_ld2.SetFillStyle(0)
        leg_ld2.AddEntry(h_hugo_ld2, "Hugo (Existing)", "lp")
        leg_ld2.AddEntry(h_ld2_nod1, "Yours (No D1 Cut)", "lp")
        leg_ld2.AddEntry(h_ld2_withd1, "Yours (With D1 Cut)", "lp")
        leg_ld2.Draw()

        c_ld2.SaveAs(f"{out_dir}/compare_LD2_xFbin{i}.pdf")
        c_ld2.SaveAs(f"{out_dir}/compare_LD2_xFbin{i}.png")
        c_ld2.Close()

    # Close files
    f_hugo.Close()
    f_nod1.Close()
    f_withd1.Close()
    print(f"\nFinished! Individual LH2 and LD2 plots generated in the '{out_dir}' directory.")

if __name__ == "__main__":
    generate_comparison_plots()