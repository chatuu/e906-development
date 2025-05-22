import ROOT
from array import array

def create_dy_plot(yield_info):
    """
    Generates a TCanvas plotting DY candidate yields (data) and DY MC prediction
    with their respective uncertainties.

    Args:
        yield_info (dict): A dictionary containing yield and error information.
                           Expected keys: "data_yield", "data_error",
                                          "mc_yield", "mc_error".
    """

    # Correcting the potential typo in the input dictionary key
    if "mc_error:" in yield_info and "mc_error" not in yield_info:
        yield_info["mc_error"] = yield_info.pop("mc_error:")

    # --- 1. Extract data from the dictionary ---
    data_yield = yield_info["data_yield"]
    data_error = yield_info["data_error"]
    mc_yield = yield_info["mc_yield"]
    mc_error = yield_info["mc_error"]

    # --- 2. Create TCanvas ---
    canvas = ROOT.TCanvas("c_dy_yields", "DY Yields: Data vs MC", 800, 600)
    canvas.SetGrid() # Add a grid for better readability
    canvas.SetLeftMargin(0.15)
    canvas.SetTickx(1) # Enable ticks on X-axis
    canvas.SetTicky(1) # Enable ticks on Y-axis
    canvas.SetBottomMargin(0.12)

    # --- 3. Prepare data for TGraphErrors ---
    # We'll place "Data" at x=1 and "MC" at x=2 on the categorical axis
    x_points = array('d', [1.0, 2.0])
    y_points = array('d', [data_yield, mc_yield])
    x_errors = array('d', [0.0, 0.0]) # No horizontal error for categorical points
    y_errors = array('d', [data_error, mc_error])

    # --- 4. Create TGraphErrors for Data ---
    # For TGraphErrors, it's often clearer to make separate graphs if styling differs significantly
    # or if they represent distinct categories as here.
    # However, for this specific case with two points, one TGraphErrors can also work.
    # Let's make two separate TGraphErrors for clarity in legend and styling.

    # Data point
    x_data = array('d', [1.0])
    y_data = array('d', [data_yield])
    ex_data = array('d', [0.0]) # No horizontal error
    ey_data = array('d', [data_error])
    gr_data = ROOT.TGraphErrors(1, x_data, y_data, ex_data, ey_data)

    gr_data.SetMarkerStyle(ROOT.kFullCircle)
    gr_data.SetMarkerSize(1.5)
    gr_data.SetMarkerColor(ROOT.kBlack)
    gr_data.SetLineColor(ROOT.kBlack)
    gr_data.SetLineWidth(2)

    # MC point
    x_mc = array('d', [2.0])
    y_mc = array('d', [mc_yield])
    ex_mc = array('d', [0.0]) # No horizontal error
    ey_mc = array('d', [mc_error])
    gr_mc = ROOT.TGraphErrors(1, x_mc, y_mc, ex_mc, ey_mc)

    gr_mc.SetMarkerStyle(ROOT.kFullSquare)
    gr_mc.SetMarkerSize(1.5)
    gr_mc.SetMarkerColor(ROOT.kBlue)
    gr_mc.SetLineColor(ROOT.kBlue)
    gr_mc.SetLineWidth(2)
    # To show MC uncertainty as a band, you could use fill style:
    # gr_mc.SetFillColorAlpha(ROOT.kBlue, 0.35)
    # gr_mc.SetFillStyle(3001) # Example fill style
    # And draw with "E2" or "a3" option later if using a band.
    # For now, we'll use error bars similar to data.

    # --- 5. Create a dummy TH1F to set up axes ---
    # This histogram will define the axis ranges and labels.
    # X-axis: "Type" (Categorical: Data, MC)
    # Y-axis: "Yield"

    # Determine Y-axis range dynamically
    min_y = min(data_yield - data_error, mc_yield - mc_error)
    max_y = max(data_yield + data_error, mc_yield + mc_error)
    padding = (max_y - min_y) * 0.2 # 20% padding

    # Create a histogram to define the axes
    # Title: "Plot Title;X-axis Title;Y-axis Title"
    # 2 bins on x-axis, from x=0.5 to x=2.5
    hist_axis = ROOT.TH1F("hist_axis", "Drell-Yan Yields;Type;Yield", 2, 0.5, 2.5)

    # Set Y-axis range
    hist_axis.SetMinimum(min_y - padding)
    hist_axis.SetMaximum(max_y + padding)

    # Set X-axis bin labels for categorical data
    hist_axis.GetXaxis().SetBinLabel(1, "Corrected Data") # Bin 1 (center x=1)
    hist_axis.GetXaxis().SetBinLabel(2, "DY MC")   # Bin 2 (center x=2)
    hist_axis.GetXaxis().SetLabelSize(0.05)
    hist_axis.GetXaxis().SetTitleSize(0.05)
    hist_axis.GetXaxis().CenterTitle(True)

    hist_axis.GetYaxis().SetTitleOffset(1.4)
    hist_axis.GetYaxis().SetLabelSize(0.04)
    hist_axis.GetYaxis().SetTitleSize(0.05)
    hist_axis.GetYaxis().CenterTitle(True)

    # --- 6. Draw the histogram (for axes) and graphs ---
    hist_axis.SetStats(ROOT.kFALSE) # Ensure no statistics box is drawn for the axis histogram
    hist_axis.Draw("AXIS") # Draw only the axes
    gr_data.Draw("P SAME") # "P" for markers, "E" for error bars (P includes E by default for TGraphErrors)
    gr_mc.Draw("P SAME")   # Draw on the same canvas

    # --- 7. Create TLegend ---
    # Adjusted coordinates for a slightly larger legend box & increased text size
    # Shifted legend left to ensure text fits within the plot region
    legend = ROOT.TLegend(0.32, 0.70, 0.75, 0.88) # x1, y1, x2, y2 (NDC coordinates)
    legend.AddEntry(gr_data, "Corrected Data (DY Candidates)", "pe") # "pe" for marker with error bar
    legend.AddEntry(gr_mc, "DY MC Prediction", "pe")
    legend.SetBorderSize(0) # Remove the border
    legend.SetTextSize(0.04) # Increase text size (adjust as needed)
    legend.SetFillStyle(0) # Transparent background
    legend.Draw()

    # --- 8. Update canvas and save ---
    canvas.Update()
    
    # You might want to save the canvas to a file
    canvas.SaveAs("dy_yields_plot.root")
    canvas.SaveAs("dy_yields_plot.pdf")
    print("Plot saved as dy_yields_plot.root and dy_yields_plot.pdf")

    # To keep the canvas open if running interactively, you might need:
    # ROOT.gApplication.Run() # or input("Press Enter to continue...")
    
    return canvas # Return the canvas object if needed elsewhere

if __name__ == '__main__':
    # Provided yield information
    yield_info = {
        "data_yield": 16140,
        "data_error": 220,
        "mc_yield": 16054,
        "mc_error": 182.8, # Corrected key
    }

    # Generate the plot
    # To prevent the canvas from disappearing immediately in batch mode,
    # ROOT.gROOT.SetBatch(True) can be set at the beginning.
    # However, for interactive display or saving, it's fine.
    
    # If you want to run in batch mode (e.g. script finishes without graphics)
    # ROOT.gROOT.SetBatch(True) 
    
    dy_plot_canvas = create_dy_plot(yield_info)

    # If not running in batch mode and you want to see the plot interactively:
    # For example, if you run `python your_script_name.py`
    # you might need to keep the script alive to see the window.
    # If ROOT is set up for interactive display, this might not be needed,
    # or you might use:
    # if not ROOT.gROOT.IsBatch():
    #     ROOT.gApplication.Run()
    # Or simply:
    # input("Press Enter to exit and close the plot...")
