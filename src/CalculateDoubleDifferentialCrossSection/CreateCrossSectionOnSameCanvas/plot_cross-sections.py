import ROOT
import os

# --- Configuration ---
# Directory containing your ROOT files
input_directory = "result_rootfiles"
# Number of xF bins to process
num_xF_bins = 16
# Output file name for the plot
output_filename = "cross_section_comparison.pdf"
output_filename1 = "cross_section_comparison.root"
# Define the xF bin ranges to be displayed as labels
# You can adjust these values to match your specific analysis bins
xF_bin_ranges = [
    (0.10, 0.15), (0.15, 0.20), (0.20, 0.25), (0.25, 0.30),
    (0.30, 0.35), (0.35, 0.40), (0.40, 0.45), (0.45, 0.50),
    (0.50, 0.55), (0.55, 0.60), (0.60, 0.65), (0.65, 0.70),
    (0.70, 0.75), (0.75, 0.80), (0.80, 0.85), (0.85, 0.90)
]

def plot_cross_sections():
    """
    This script reads histograms from multiple ROOT files, scales them by an offset,
    and plots them on a single TCanvas with a logarithmic y-axis.
    """
    # Set a clean plot style
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    
    # --- Style settings for a centered title ---
    ROOT.gStyle.SetTitleAlign(23) # Horizontally center (2) and vertically top (3)
    ROOT.gStyle.SetTitleX(0.5)    # Set horizontal position of title to center of pad
    ROOT.gStyle.SetTitleY(0.99)   # Set vertical position of title at the top of the pad
    ROOT.gStyle.SetTitleH(0.04)   # Set title height
    ROOT.gStyle.SetTitleBorderSize(0) # Remove the border around the title

    # Create a canvas to draw on
    canvas = ROOT.TCanvas("canvas", "Cross-Section Comparison", 1200, 900)
    canvas.SetLogy()
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)

    # Create a legend to identify the different xF bins. Made taller for better spacing.
    legend = ROOT.TLegend(0.75, 0.45, 0.9, 0.9)
    legend.SetHeader("x_{F} Bins")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.025) # Adjust text size for readability

    # Lists to store objects to prevent them from being garbage collected
    processed_hists = []
    latex_labels = []

    # Define a list of colors to use for different bins
    colors = [
        ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kCyan,
        ROOT.kOrange + 7, ROOT.kSpring + 5, ROOT.kTeal + 5, ROOT.kAzure + 1,
        ROOT.kGray + 2, ROOT.kViolet - 5, ROOT.kYellow + 2, ROOT.kPink + 1,
        ROOT.kGreen - 9, ROOT.kRed - 9
    ]

    # --- First pass: Loop over all xF bins to load and scale ---
    for i in range(num_xF_bins):
        file_path = os.path.join(input_directory, f"histograms_xF_{i}.root")

        if not os.path.exists(file_path):
            print(f"Warning: File not found, skipping: {file_path}")
            continue

        root_file = ROOT.TFile.Open(file_path, "READ")
        if not root_file or root_file.IsZombie():
            print(f"Error: Could not open file: {file_path}")
            continue

        hData = root_file.Get("hData")
        hDataOverlay = root_file.Get("hDataOverlay")

        if not hData or not hDataOverlay:
            print(f"Error: Histograms not found in file: {file_path}")
            root_file.Close()
            continue

        hData.SetDirectory(0)
        hDataOverlay.SetDirectory(0)
        root_file.Close()

        scale_factor = 5 * (10**(2*i))
        hData.Scale(scale_factor)
        hDataOverlay.Scale(scale_factor)
        
        # Store processed histograms, their original index, and the scale factor
        processed_hists.append({'data': hData, 'sys': hDataOverlay, 'index': i, 'scale': scale_factor})

    # --- Drawing logic ---
    if not processed_hists:
        print("Error: No histograms were loaded. Cannot create plot.")
        return

    # Create a frame histogram to define the axes and their ranges
    h_frame = canvas.DrawFrame(2.8, 1e-3, 9.0, 1e35)

    # Set titles and labels on the frame
    h_frame.SetTitle("DY Absolute Cross-Section Vs Mass for x_{F} bins")
    h_frame.GetXaxis().SetTitle("Invariant Mass (GeV)")
    h_frame.GetXaxis().CenterTitle()
    h_frame.GetXaxis().SetTitleOffset(1.2)
    h_frame.GetYaxis().SetTitle("M^{3} #frac{d^{2}#sigma}{dMdx_{F}} (nb GeV^{2}) #times 5 #times 10^{x_{F} bin index}")
    h_frame.GetYaxis().CenterTitle()
    h_frame.GetYaxis().SetTitleOffset(1.4)
    
    # Loop through and draw all processed histograms
    for hist_info in processed_hists:
        hData = hist_info['data']
        hDataOverlay = hist_info['sys']
        bin_index = hist_info['index']
        scale_factor = hist_info['scale']
        color = colors[bin_index % len(colors)]

        # --- Set visual styles ---
        hDataOverlay.SetLineColor(color)
        hDataOverlay.SetFillColorAlpha(color, 0.35)
        hDataOverlay.SetFillStyle(1001)
        hDataOverlay.SetMarkerSize(0)

        hData.SetLineColor(color)
        hData.SetMarkerColor(color)
        hData.SetMarkerStyle(ROOT.kFullCircle)
        hData.SetMarkerSize(1.0)

        # Draw on top of the frame using the "SAME" option
        hDataOverlay.Draw("E2 SAME")
        hData.Draw("E1 P SAME")

        # Add TLatex label for the xF bin range
        if bin_index < len(xF_bin_ranges):
            low_edge, high_edge = xF_bin_ranges[bin_index]
            y_pos = 1.0 * scale_factor
            latex = ROOT.TLatex(3.1, y_pos, f"{low_edge:.2f} #leq x_{{F}} < {high_edge:.2f}")
            latex.SetTextSize(0.025)
            latex.SetTextColor(color)
            latex.Draw()
            latex_labels.append(latex) # Store to prevent garbage collection

        legend.AddEntry(hData, f"x_{{F}} bin {bin_index}", "pl")

    #legend.Draw()
    canvas.Update()
    canvas.SaveAs(output_filename)
    canvas.SaveAs(output_filename1)

    print(f"Plot saved to {output_filename}")
    print(f"Plot saved to {output_filename1}")

if __name__ == "__main__":
    plot_cross_sections()

