import ROOT
import os

def main():
    # Run in batch mode to prevent X11 windows from opening during plotting
    ROOT.gROOT.SetBatch(True)
    
    # Hide ROOT statistics box and set plain style
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    # File paths (using os.path.expanduser to resolve the '~')
    path_no_d1 = os.path.expanduser("~/github/e906-development/src/xsec_pT/RS57-70_NoD1Cut/All_XSec_Objects.root")
    path_with_d1 = os.path.expanduser("~/github/e906-development/src/xsec_pT/RS57-70/All_XSec_Objects.root")

    # Open ROOT files
    f_no_d1 = ROOT.TFile.Open(path_no_d1, "READ")
    f_with_d1 = ROOT.TFile.Open(path_with_d1, "READ")

    if not f_no_d1 or f_no_d1.IsZombie():
        print(f"Error: Could not open {path_no_d1}")
        return
    if not f_with_d1 or f_with_d1.IsZombie():
        print(f"Error: Could not open {path_with_d1}")
        return

    targets = ["LH2", "LD2"]

    # Define the requested axis limits
    x_min, x_max = 0.0, 2.0
    ratio_y_min, ratio_y_max = 0.95, 1.05

    for target in targets:
        # --- Setup Canvas and Pads ---
        c1 = ROOT.TCanvas(f"c_{target}", f"{target} Cross Section Comparison", 800, 800)
        
        # Main Upper Pad
        pad1 = ROOT.TPad(f"pad1_{target}", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.12)
        pad1.SetTickx(1)
        pad1.SetTicky(1)
        pad1.Draw()
        
        # Lower Ratio Pad
        pad2 = ROOT.TPad(f"pad2_{target}", "pad2", 0, 0.0, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.3)
        pad2.SetLeftMargin(0.12)
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.Draw()

        dir_name = f"CrossSections_{target}"
        h_xsec_name = f"h1_xsec_{target}_geom"
        h_sys_name  = f"h1_sys_{target}_geom"

        # Extract "No D1 Cut" histograms
        h_xsec_no_d1 = f_no_d1.Get(f"{dir_name}/{h_xsec_name}")
        h_sys_no_d1  = f_no_d1.Get(f"{dir_name}/{h_sys_name}")

        # Extract "With D1 Cut" histograms
        h_xsec_with_d1 = f_with_d1.Get(f"{dir_name}/{h_xsec_name}")
        h_sys_with_d1  = f_with_d1.Get(f"{dir_name}/{h_sys_name}")

        # Safety check to ensure histograms exist
        if not all([h_xsec_no_d1, h_sys_no_d1, h_xsec_with_d1, h_sys_with_d1]):
            print(f"Warning: Missing histograms for target {target}. Skipping.")
            continue

        # ==========================================
        # UPPER PAD: Main Plot
        # ==========================================
        pad1.cd()

        # Formatting "No D1 Cut" (Black points, Gray band)
        h_sys_no_d1.SetFillColorAlpha(ROOT.kBlack, 0.3)
        h_sys_no_d1.SetMarkerStyle(0)
        h_sys_no_d1.SetLineColor(ROOT.kBlack)
        
        h_xsec_no_d1.SetMarkerStyle(20)
        h_xsec_no_d1.SetMarkerColor(ROOT.kBlack)
        h_xsec_no_d1.SetLineColor(ROOT.kBlack)

        # Formatting "With D1 Cut" (Red points, Light Red band)
        h_sys_with_d1.SetFillColorAlpha(ROOT.kRed, 0.3)
        h_sys_with_d1.SetMarkerStyle(0)
        h_sys_with_d1.SetLineColor(ROOT.kRed)

        h_xsec_with_d1.SetMarkerStyle(21)
        h_xsec_with_d1.SetMarkerColor(ROOT.kRed)
        h_xsec_with_d1.SetLineColor(ROOT.kRed)

        # Draw a strict frame to perfectly lock the X and Y axes
        max_y = max(h_sys_no_d1.GetMaximum(), h_sys_with_d1.GetMaximum())
        frame1 = pad1.DrawFrame(x_min, 0.0, x_max, max_y * 1.5)
        
        frame1.GetYaxis().SetTitle("d#sigma / dp_{T} (nb / GeV / Nucleus)")
        frame1.GetYaxis().SetTitleSize(0.04)
        frame1.GetYaxis().CenterTitle(True)
        
        # Hide X-axis labels on the top pad
        frame1.GetXaxis().SetLabelSize(0)
        frame1.GetXaxis().SetTitleSize(0)

        # Draw Logic (Over top of the strict frame)
        h_sys_no_d1.Draw("E2 SAME")
        h_sys_with_d1.Draw("E2 SAME")
        h_xsec_no_d1.Draw("PE1 SAME")
        h_xsec_with_d1.Draw("PE1 SAME")

        # Legend
        leg = ROOT.TLegend(0.45, 0.65, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.035)
        leg.SetHeader(f"Target: {target}", "C")
        
        leg.AddEntry(h_xsec_no_d1, "Without D1 Cut (Stat. Err.)", "pe")
        leg.AddEntry(h_sys_no_d1, "Without D1 Cut (Sys. Err.)", "f")
        leg.AddEntry(h_xsec_with_d1, "After D1 Cut (Stat. Err.)", "pe")
        leg.AddEntry(h_sys_with_d1, "After D1 Cut (Sys. Err.)", "f")
        leg.Draw()

        # ==========================================
        # LOWER PAD: Ratio Plot
        # ==========================================
        pad2.cd()
        
        # Create ratio histogram: (With D1) / (No D1)
        h_ratio = h_xsec_with_d1.Clone(f"h_ratio_{target}")
        h_ratio.Divide(h_xsec_no_d1)
        
        # Fit horizontal line (pol0) strictly within the x_min to x_max bounds
        h_ratio.Fit("pol0", "Q0", "", x_min, x_max)
        func = h_ratio.GetFunction("pol0")
        if func:
            p0 = func.GetParameter(0)
            err0 = func.GetParError(0)
        else:
            p0, err0 = 1.0, 0.0

        # Draw a strict frame to completely force the bounds [0.0, 2.0] and [0.95, 1.05]
        frame2 = pad2.DrawFrame(x_min, ratio_y_min, x_max, ratio_y_max)
        frame2.GetYaxis().SetTitle("Ratio (D1 / No D1)")
        frame2.GetXaxis().SetTitle("p_{T} (GeV)")

        # Text sizes 
        frame2.GetYaxis().SetTitleSize(0.06) 
        frame2.GetYaxis().SetLabelSize(0.06) 
        frame2.GetYaxis().SetTitleOffset(0.75) 
        frame2.GetYaxis().SetNdivisions(505)
        frame2.GetYaxis().CenterTitle(True)

        frame2.GetXaxis().SetTitleSize(0.12)
        frame2.GetXaxis().SetLabelSize(0.1)
        frame2.GetXaxis().SetTitleOffset(0.9)
        frame2.GetXaxis().CenterTitle(True)

        # Set ratio styling
        h_ratio.SetMarkerStyle(20)
        h_ratio.SetMarkerColor(ROOT.kBlack)
        h_ratio.SetLineColor(ROOT.kBlack)

        # 1. Draw Ratio Points First
        h_ratio.Draw("PE1 SAME")

        # 2. Draw Error Band using TBox 
        box = ROOT.TBox(x_min, p0 - err0, x_max, p0 + err0)
        box.SetFillColorAlpha(ROOT.kBlue, 0.3)
        box.Draw("SAME")

        # 3. Draw Best Fit Line
        line = ROOT.TLine(x_min, p0, x_max, p0)
        line.SetLineColor(ROOT.kBlue)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw("SAME")

        # 4. Redraw Ratio Points on top of the band
        h_ratio.Draw("PE1 SAME")

        # 5. Add TLatex for Fit Value
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        latex.SetTextSize(0.08)
        latex.SetTextColor(ROOT.kBlue)
        
        # Position label cleanly away from the axes
        latex.DrawLatex(0.65, 0.85, f"Fit: {p0:.3f} #pm {err0:.3f}")

        # Save to PDF
        c1.cd()
        pdf_name = f"Compare_XSec_{target}.pdf"
        c1.SaveAs(pdf_name)
        print(f"Saved: {pdf_name}")

    # Close files
    f_no_d1.Close()
    f_with_d1.Close()

if __name__ == "__main__":
    main()