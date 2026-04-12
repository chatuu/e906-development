import pandas as pd
import numpy as np
import ROOT
import os

def main():
    # Force ROOT into batch mode so it doesn't pop up windows for every plot
    ROOT.gROOT.SetBatch(True)

    # File names
    file_red = "hodoscope_eff_RS57-70_final.tsv"
    #file_blue = "hodoscope_eff_RS67_final.tsv"
    file_blue = "hodoscope_eff_RS67_final_2026-01-31.tsv"  # Updated file with sample std deviation
    
    # Check if files exist
    if not os.path.exists(file_red) or not os.path.exists(file_blue):
        print("Error: One or both of the .tsv files are missing in the current directory.")
        return

    # Column layout based on your data
    columns = ['Plane', 'Paddle', 'Efficiency', 'ErrUp', 'ErrDown']

    print("Loading data...")
    # sep='\s+' handles spaces or tabs seamlessly
    df_red = pd.read_csv(file_red, sep='\s+', header=None, names=columns)
    df_blue = pd.read_csv(file_blue, sep='\s+', header=None, names=columns)

    planes = ["H1T", "H1B", "H2T", "H2B", "H3T", "H3B", "H4T", "H4B"]

    print("Generating ROOT PDF plots...")
    for plane in planes:
        data_red = df_red[df_red['Plane'] == plane]
        data_blue = df_blue[df_blue['Plane'] == plane]
        
        if data_red.empty and data_blue.empty:
            continue
            
        # 1. Create TCanvas and apply requested styling
        c1 = ROOT.TCanvas(f"c_{plane}", f"Efficiency {plane}", 800, 600)
        c1.SetTickx(1)  # Ticks on top and bottom
        c1.SetTicky(1)  # Ticks on left and right
        
        # We will use a TMultiGraph to draw multiple TGraphs on the same pad
        mg = ROOT.TMultiGraph()
        mg.SetTitle(f"Hodoscope Paddle Efficiency - {plane};Paddle ID;Efficiency")
        
        # Position at bottom center and remove border
        leg = ROOT.TLegend(0.3, 0.15, 0.7, 0.25)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)

        # 2. Add Blue Data (Old)
        if not data_blue.empty:
            n_b = len(data_blue)
            x_b = np.array(data_blue['Paddle'], dtype=float)
            y_b = np.array(data_blue['Efficiency'], dtype=float)
            eyl_b = np.array(data_blue['ErrDown'], dtype=float)
            eyh_b = np.array(data_blue['ErrUp'], dtype=float)
            ex_b = np.zeros(n_b, dtype=float)
            
            gr_blue = ROOT.TGraphAsymmErrors(n_b, x_b, y_b, ex_b, ex_b, eyl_b, eyh_b)
            gr_blue.SetMarkerColor(ROOT.kBlue)
            gr_blue.SetLineColor(ROOT.kBlue)
            gr_blue.SetMarkerStyle(ROOT.kFullCircle) # 20
            gr_blue.SetMarkerSize(1.0)
            
            mg.Add(gr_blue, "P")
            leg.AddEntry(gr_blue, "RS67 Final (Old)", "lp")

        # 3. Add Red Data (New)
        if not data_red.empty:
            n_r = len(data_red)
            x_r = np.array(data_red['Paddle'], dtype=float)
            y_r = np.array(data_red['Efficiency'], dtype=float)
            eyl_r = np.array(data_red['ErrDown'], dtype=float)
            eyh_r = np.array(data_red['ErrUp'], dtype=float)
            ex_r = np.zeros(n_r, dtype=float)
            
            gr_red = ROOT.TGraphAsymmErrors(n_r, x_r, y_r, ex_r, ex_r, eyl_r, eyh_r)
            gr_red.SetMarkerColor(ROOT.kRed)
            gr_red.SetLineColor(ROOT.kRed)
            gr_red.SetMarkerStyle(ROOT.kFullSquare) # 21
            gr_red.SetMarkerSize(1.0)
            
            mg.Add(gr_red, "P")
            leg.AddEntry(gr_red, "RS57-70 Final (New)", "lp")

        # 4. Draw the MultiGraph 
        mg.Draw("AP")
        
        # ROOT idiosyncrasy: Axes on a TMultiGraph only exist AFTER Draw() is called. 
        ROOT.gPad.Update()
        
        # 5. Apply Axis Formatting (Center titles & adjust limits)
        if mg.GetXaxis():
            mg.GetXaxis().CenterTitle(True)
            max_x = max(data_red['Paddle'].max() if not data_red.empty else 0,
                        data_blue['Paddle'].max() if not data_blue.empty else 0)
            mg.GetXaxis().SetLimits(0, max_x + 1)
            
        if mg.GetYaxis():
            mg.GetYaxis().CenterTitle(True)
            
            # Apply custom Y-axis ranges based on the specific plane
            if plane in ["H3T", "H3B", "H4T", "H4B"]:
                mg.SetMinimum(0.8)
                mg.SetMaximum(1.02)
            else:
                mg.SetMinimum(-0.05)
                mg.SetMaximum(1.1)

        leg.Draw()
        
        # 6. Save Canvas as PDF
        output_filename = f"eff_plot_{plane}_root.pdf"
        c1.SaveAs(output_filename)
        
        print(f" -> Saved {output_filename}")

if __name__ == "__main__":
    main()