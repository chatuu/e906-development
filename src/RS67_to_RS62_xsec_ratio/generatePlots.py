import ROOT
import os

def main():
    # Run in batch mode to prevent canvases from popping up and slowing down the script
    ROOT.gROOT.SetBatch(True)
    
    # Create an output directory for the plots
    output_dir = "comparison_plots"
    os.makedirs(output_dir, exist_ok=True)

    # 1. Open the ROOT files
    file_RS67 = ROOT.TFile.Open("./root_files/RS67/All_XSec_Objects.root", "READ")
    file_RS5770 = ROOT.TFile.Open("./root_files/RS57-70/All_XSec_Objects.root", "READ")

    if not file_RS67 or file_RS67.IsZombie() or not file_RS5770 or file_RS5770.IsZombie():
        print("Error: Could not open one or both ROOT files. Check the paths.")
        return

    targets = ["LH2", "LD2"]
    
    # DEFINE YOUR xF BIN EDGES HERE (16 bins = 17 edges)
    xf_edges = [
        0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
        0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80
    ]
    
    # 2. Loop over targets and xF bins
    for target in targets:
        dir_name = f"CrossSections_{target}"
        
        # Get directories
        dir_RS67 = file_RS67.Get(dir_name)
        dir_RS5770 = file_RS5770.Get(dir_name)
        
        if not dir_RS67 or not dir_RS5770:
            print(f"Error: Directory {dir_name} not found in one of the files.")
            continue

        for xf_bin in range(16):  # xF bins from 0 to 15
            hist_name = f"h1_xsec_cent_{target}_xF_{xf_bin}"
            
            # Fetch histograms
            h_RS67 = dir_RS67.Get(hist_name)
            h_RS5770 = dir_RS5770.Get(hist_name)
            
            if not h_RS67 or not h_RS5770:
                print(f"Warning: {hist_name} missing in one or both files. Skipping.")
                continue

            # Ensure sum of squares of weights is calculated for proper error propagation
            if not h_RS67.GetSumw2N(): h_RS67.Sumw2()
            if not h_RS5770.GetSumw2N(): h_RS5770.Sumw2()

            # 3. Create Canvas and Pads
            canvas_name = f"c_{target}_xF_{xf_bin}"
            c = ROOT.TCanvas(canvas_name, f"Cross Section {target} xF {xf_bin}", 800, 800)

            # --- Top Pad: Cross Sections ---
            pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.33, 1.0, 1.0)
            pad1.SetBottomMargin(0.15)  
            pad1.SetLeftMargin(0.15)
            pad1.Draw()
            pad1.cd()
            
            # Apply tick settings requested
            pad1.SetTickx(1)
            pad1.SetTicky(1)

            # Styling for RS67 (Blue)
            h_RS67.SetLineColor(ROOT.kBlue)
            h_RS67.SetMarkerColor(ROOT.kBlue)
            h_RS67.SetMarkerStyle(20)
            h_RS67.SetStats(0)
            
            # Title
            xf_min = xf_edges[xf_bin]
            xf_max = xf_edges[xf_bin + 1]
            h_RS67.SetTitle(f"{target} Cross Section, {xf_min} #leq x_{{F}} < {xf_max}")
            
            # Top Y-axis styling & centering
            h_RS67.GetYaxis().SetTitle("Cross Section")
            h_RS67.GetYaxis().CenterTitle(True)
            h_RS67.GetYaxis().SetTitleSize(20)
            h_RS67.GetYaxis().SetTitleFont(43)
            h_RS67.GetYaxis().SetTitleOffset(1.8)
            h_RS67.GetYaxis().SetLabelFont(43)
            h_RS67.GetYaxis().SetLabelSize(15)

            # Top X-axis styling & centering
            h_RS67.GetXaxis().SetTitle("Invariant Mass (GeV)")
            h_RS67.GetXaxis().CenterTitle(True)
            h_RS67.GetXaxis().SetTitleSize(20)
            h_RS67.GetXaxis().SetTitleFont(43)
            h_RS67.GetXaxis().SetTitleOffset(1.8)
            h_RS67.GetXaxis().SetLabelFont(43)
            h_RS67.GetXaxis().SetLabelSize(15)

            # Styling for RS57-70 (Red)
            h_RS5770.SetLineColor(ROOT.kRed)
            h_RS5770.SetMarkerColor(ROOT.kRed)
            h_RS5770.SetMarkerStyle(24) # Open circle
            h_RS5770.SetStats(0)

            # Adjust Y-axis maximum to fit both histograms
            max_y = max(h_RS67.GetMaximum(), h_RS5770.GetMaximum())
            h_RS67.SetMaximum(max_y * 1.3)

            # Draw top histograms
            h_RS67.Draw("E1 P")
            h_RS5770.Draw("E1 P SAME")

            # Add Legend
            leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.85)
            leg.SetBorderSize(0)
            leg.AddEntry(h_RS67, "RS67", "lep")
            leg.AddEntry(h_RS5770, "RS62", "lep")
            leg.Draw()

            c.cd()  # Go back to main canvas before drawing pad2

            # --- Bottom Pad: Ratio ---
            pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.33)
            pad2.SetTopMargin(0.05)
            pad2.SetBottomMargin(0.35)
            pad2.SetLeftMargin(0.15)
            pad2.Draw()
            pad2.cd()
            
            # Apply tick settings requested
            pad2.SetTickx(1)
            pad2.SetTicky(1)

            # Calculate Ratio (RS67 / RS57-70)
            h_ratio = h_RS67.Clone(f"ratio_{target}_xF_{xf_bin}")
            h_ratio.SetTitle("")
            h_ratio.Divide(h_RS5770)
            
            # Fit the ratio *before* calculating limits to include fit constraints
            h_ratio.Fit("pol0", "Q0")
            fit_func = h_ratio.GetFunction("pol0")
            
            # --- Dynamic Y-Axis Calculation ---
            min_y_val = float('inf')
            max_y_val = float('-inf')
            
            # 1. Loop through bins to check data bounds (value +/- error)
            for i in range(1, h_ratio.GetNbinsX() + 1):
                content = h_ratio.GetBinContent(i)
                error = h_ratio.GetBinError(i)
                # Ignore empty bins
                if content != 0 or error != 0:
                    min_y_val = min(min_y_val, content - error)
                    max_y_val = max(max_y_val, content + error)
            
            # 2. Factor in the fit line and error band
            if fit_func:
                f_val = fit_func.GetParameter(0)
                f_err = fit_func.GetParError(0)
                min_y_val = min(min_y_val, f_val - f_err)
                max_y_val = max(max_y_val, f_val + f_err)

            # 3. Ensure the Ratio=1 reference line is always well within bounds
            min_y_val = min(min_y_val, 0.95)
            max_y_val = max(max_y_val, 1.05)

            # 4. Calculate final padded ranges (adding 20% margin top and bottom)
            if min_y_val != float('inf') and max_y_val != float('-inf'):
                y_span = max_y_val - min_y_val
                if y_span == 0: y_span = 0.1 # Fallback for perfectly flat line
                final_min = min_y_val - (0.20 * y_span)
                final_max = max_y_val + (0.20 * y_span)
            else:
                final_min, final_max = 0.5, 1.5 # Safe fallback

            # Styling for Ratio
            h_ratio.SetLineColor(ROOT.kBlack)
            h_ratio.SetMarkerColor(ROOT.kBlack)
            h_ratio.SetMarkerStyle(20)
            h_ratio.SetStats(0)

            # Ratio Y-axis styling & centering
            h_ratio.GetYaxis().SetTitle("RS67 / RS62")
            h_ratio.GetYaxis().CenterTitle(True)
            
            # APPLIED DYNAMIC RANGE HERE
            h_ratio.GetYaxis().SetRangeUser(final_min, final_max)
            
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(20)
            h_ratio.GetYaxis().SetTitleFont(43)
            h_ratio.GetYaxis().SetTitleOffset(1.8)
            h_ratio.GetYaxis().SetLabelFont(43)
            h_ratio.GetYaxis().SetLabelSize(15)

            # Ratio X-axis styling & centering
            h_ratio.GetXaxis().SetTitle("Invariant Mass (GeV)") 
            h_ratio.GetXaxis().CenterTitle(True)
            h_ratio.GetXaxis().SetTitleSize(20)
            h_ratio.GetXaxis().SetTitleFont(43)
            h_ratio.GetXaxis().SetTitleOffset(3.2)
            h_ratio.GetXaxis().SetLabelFont(43)
            h_ratio.GetXaxis().SetLabelSize(15)

            # Draw ratio histogram first to establish axes
            h_ratio.Draw("E1 P")
            
            if fit_func:
                fit_val = fit_func.GetParameter(0)
                fit_err = fit_func.GetParError(0)
                
                x_min = h_ratio.GetXaxis().GetXmin()
                x_max = h_ratio.GetXaxis().GetXmax()
                
                # Create and draw the error band (TBox)
                error_box = ROOT.TBox(x_min, fit_val - fit_err, x_max, fit_val + fit_err)
                error_box.SetFillColorAlpha(ROOT.kPink, 0.4) # Transparent pink
                error_box.Draw("SAME")

                # Style and draw the central fit line
                fit_func.SetLineColor(ROOT.kRed)
                fit_func.SetLineStyle(1)
                fit_func.SetLineWidth(2)
                fit_func.Draw("SAME")
                
                # Re-draw data points so they are on top of the error band
                h_ratio.Draw("E1 P SAME")
                
                # Draw the TLatex text dynamically
                latex = ROOT.TLatex()
                latex.SetTextFont(43)
                latex.SetTextSize(16)
                latex.SetTextColor(ROOT.kRed)
                
                x_pos = x_min + 0.05 * (x_max - x_min)
                
                # If the fit is in the top 20% of the newly calculated range, draw text below it
                y_span_final = final_max - final_min
                if (final_max - fit_val) < (0.2 * y_span_final):
                    y_pos = fit_val - (0.15 * y_span_final)
                else:
                    y_pos = fit_val + (0.05 * y_span_final)
                    
                latex.DrawLatex(x_pos, y_pos, f"y = {fit_val:.3f} #pm {fit_err:.3f}")

            # Add a reference line at Ratio = 1
            line = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1.0, h_ratio.GetXaxis().GetXmax(), 1.0)
            line.SetLineColor(ROOT.kGray + 2)
            line.SetLineStyle(2)
            line.Draw("SAME")

            # 4. Save Canvas
            output_file = os.path.join(output_dir, f"XSec_Compare_{target}_xF_{xf_bin}.pdf")
            c.SaveAs(output_file)
            
            # Cleanup to prevent memory leaks
            c.Close()

    # Close files
    file_RS67.Close()
    file_RS5770.Close()
    print(f"\nProcessing complete. Plots saved in the '{output_dir}' directory.")

if __name__ == "__main__":
    main()