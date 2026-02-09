import ROOT
import csv

def create_hodoscope_histograms(filename):
    """
    Reads TSV file and creates ROOT histograms and TGraphs.
    Generates PDF plots with dynamic Y-axis ranges based on error bars.
    """
    
    # Dictionary to store data: { 'Label': [ {id, eff, up, low}, ... ] }
    data_map = {}

    print(f"Reading file: {filename}")
    
    try:
        with open(filename, 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            
            for row in reader:
                if not row: continue # Skip empty lines
                
                # Clean empty strings
                parts = [p for p in row if p.strip()]
                
                if len(parts) < 5:
                    continue

                # Skip header line if it exists
                try:
                    elem_id = int(parts[1])
                except ValueError:
                    continue

                label = parts[0]
                eff = float(parts[2])
                err_up = float(parts[3])
                err_low = float(parts[4])

                if label not in data_map:
                    data_map[label] = []

                data_map[label].append({
                    'id': elem_id,
                    'eff': eff,
                    'up': err_up,
                    'low': err_low
                })
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return

    # Create output ROOT file
    out_file = ROOT.TFile("hodoscope_results.root", "RECREATE")
    
    print(f"Found labels: {list(data_map.keys())}")

    # Create Histograms and Graphs
    for label, points in data_map.items():
        
        # 1. Determine Binning and Y-Axis Limits
        ids = [p['id'] for p in points]
        max_id = max(ids)
        
        # Calculate the absolute min and max including error bars
        # min_y = lowest (efficiency - low_error)
        # max_y = highest (efficiency + up_error)
        y_vals_low = [p['eff'] - p['low'] for p in points]
        y_vals_high = [p['eff'] + p['up'] for p in points]
        
        min_y_limit = min(y_vals_low)
        max_y_limit = max(y_vals_high)

        # Calculate a small margin (5%) so error bars don't touch the frame edge
        y_span = max_y_limit - min_y_limit
        margin = y_span * 0.05
        if margin == 0: margin = 0.01  # Safety if all points are identical

        final_y_min = min_y_limit - margin
        final_y_max = max_y_limit + margin

        # 2. Create TH1D (Histogram)
        h_name = f"{label}"
        h_title = f"{label};Element ID;Efficiency"
        h1 = ROOT.TH1D(h_name, h_title, max_id, 0.5, max_id + 0.5)
        
        # 3. Create TGraphAsymmErrors (Graph)
        g_name = f"{label}_graph"
        g_title = f"{label} Efficiency;Element ID;Efficiency"
        tg = ROOT.TGraphAsymmErrors()
        tg.SetName(g_name)
        tg.SetTitle(g_title)

        # --- Visual Styling ---
        tg.SetMarkerStyle(20)   # 20 = Full Circle
        tg.SetMarkerSize(1.0)
        
        point_index = 0
        
        for p in points:
            eid = p['id']
            eff = p['eff']
            e_up = p['up']
            e_low = p['low']
            
            # Fill Histogram
            h1.SetBinContent(eid, eff)
            avg_err = (e_up + e_low) / 2.0
            h1.SetBinError(eid, avg_err)
            
            # Fill Graph
            tg.SetPoint(point_index, float(eid), eff)
            tg.SetPointError(point_index, 0, 0, e_low, e_up)
            
            point_index += 1

        # Save Histogram
        h1.Write()

        # --- Plotting and Saving Canvas ---
        c = ROOT.TCanvas(f"c_{label}", f"Canvas for {label}", 800, 600)
        c.SetTickx(1)
        c.SetTicky(1)
        
        # Draw first so the axis object is created
        tg.Draw("AP")
        
        # --- AXIS MODIFICATIONS ---
        # 1. Force X-axis to start at 0
        tg.GetXaxis().SetLimits(0, max_id + 1)
        
        # 2. Dynamic Y-axis based on Error Bars
        tg.GetYaxis().SetRangeUser(final_y_min, final_y_max)
        
        tg.GetXaxis().CenterTitle(1)
        tg.GetYaxis().CenterTitle(1)
        
        # Update canvas to apply changes
        c.Update()
        
        tg.Write()
        c.SaveAs(f"{label}_efficiency.pdf")
        c.Write()
        
        print(f"Created Plot: {label} (Y-axis: {final_y_min:.3f} to {final_y_max:.3f})")

    out_file.Close()
    print("------------------------------------------------")
    print("Processing complete.")
    print("Output saved to: hodoscope_results.root")

if __name__ == "__main__":
    # Ensure you have your TSV file named correctly here
    create_hodoscope_histograms("hodoscope_eff_RS67_final.tsv")