import sys
import os
import ROOT
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle

def get_integral_above_cut(hist, cut_value):
    """Return the integral above the given cut value."""
    if not hist:
        return 0
    bin_cut = hist.GetXaxis().FindBin(cut_value)
    return hist.Integral(bin_cut, hist.GetNbinsX())

def main():
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <root_file> <cut_value>")
        sys.exit(1)

    root_file_path = sys.argv[1]
    cut_value = float(sys.argv[2])
    file_name = os.path.basename(root_file_path)

    # Open ROOT file
    f = ROOT.TFile.Open(root_file_path)
    if not f or f.IsZombie():
        print(f"Error: Cannot open file {root_file_path}")
        sys.exit(1)

    hist_names = ["data", "dy", "jpsi", "psip", "mix", "flask"]
    integrals_above = {}

    # Get integrals for each histogram
    for name in hist_names:
        hist = f.Get(name)
        if not hist:
            print(f"Warning: Histogram '{name}' not found in file.")
            integrals_above[name] = 0
        else:
            integrals_above[name] = get_integral_above_cut(hist, cut_value)

    f.Close()

    # Calculate percentage relative to data
    data_integral_above = integrals_above["data"]
    percentages = {
        name: (val / data_integral_above * 100 if data_integral_above != 0 else 0)
        for name, val in integrals_above.items()
    }

    # Prepare table data
    table_data = [["Component", f"Integral > {cut_value}", "Percentage to data (%)"]]
    for name in hist_names:
        table_data.append([
            name,
            f"{integrals_above[name]:.3f}",
            f"{percentages[name]:.3f}"
        ])

    # Save as PDF
    pdf_filename = f"integrals_table_{os.path.splitext(file_name)[0]}.pdf"
    doc = SimpleDocTemplate(pdf_filename, pagesize=letter)
    table = Table(table_data)
    style = TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.grey),
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('BOTTOMPADDING', (0,0), (-1,0), 12),
        ('BACKGROUND', (0,1), (-1,-1), colors.beige),
        ('GRID', (0,0), (-1,-1), 1, colors.black),
    ])
    table.setStyle(style)
    doc.build([table])

    print(f"PDF saved as: {pdf_filename}")

if __name__ == "__main__":
    main()
