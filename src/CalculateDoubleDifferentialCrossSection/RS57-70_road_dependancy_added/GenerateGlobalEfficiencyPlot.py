import numpy as np
import ROOT

def main():
    # Set ROOT to batch mode so it doesn't try to open X11 windows
    ROOT.gROOT.SetBatch(True)
    
    # 1. Load the .npz file (replace with your actual file name)
    filename = "interpolation_data_d1.npz"
    
    try:
        data = np.load(filename)
    except FileNotFoundError:
        print(f"Error: Could not find '{filename}'. Please check the path.")
        return

    # Print the keys in the .npz file so you know what names to use below
    print("Keys found in the .npz file:", data.files)

    # 2. Extract the arrays
    # IMPORTANT: Change the string keys ('d1', 'efficiency', etc.) to exactly match 
    # the keys printed by data.files above!
    x       = np.ascontiguousarray(data['x'], dtype=np.float64)
    y       = np.ascontiguousarray(data['y'], dtype=np.float64)
    ey_low  = np.ascontiguousarray(data['y_error_low'], dtype=np.float64)
    ey_high = np.ascontiguousarray(data['y_error_high'], dtype=np.float64)

    # Get number of points
    n = len(x)

    # Since no x-errors were provided, we create arrays of zeros for them
    ex = np.zeros(n, dtype=np.float64)

    # 3. Create the TGraphAsymmErrors object
    # Constructor: n_points, x, y, ex_low, ex_high, ey_low, ey_high
    graph = ROOT.TGraphAsymmErrors(n, x, y, ex, ex, ey_low, ey_high)

    # 4. Format the plot
    graph.SetTitle("Efficiency vs D1;D1;Efficiency")
    graph.SetMarkerStyle(20)        # Filled circle
    graph.SetMarkerSize(1.2)
    graph.SetMarkerColor(ROOT.kBlue)
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(2)

    # 5. Draw and save
    canvas = ROOT.TCanvas("c1", "Efficiency Plot", 800, 600)
    
    # Enable tick marks on top and right sides
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    
    # "A" draws the axes, "P" draws the markers (data points)
    graph.Draw("AP")

    # Center the axis labels/titles (Done after Draw() to ensure axes are fully initialized)
    graph.GetXaxis().CenterTitle(True)
    graph.GetYaxis().CenterTitle(True)

    # Optional: adjust the y-axis range if you want specific boundaries
    # graph.GetYaxis().SetRangeUser(0.0, 1.1)
    
    # Update the canvas to apply axis configurations before saving
    canvas.Update()

    output_file = "efficiency_d1.pdf"
    canvas.SaveAs(output_file)
    print(f"Successfully created: {output_file}")

if __name__ == "__main__":
    main()