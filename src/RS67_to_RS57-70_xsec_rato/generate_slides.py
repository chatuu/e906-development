import os
import csv

def format_tex_sci(val, bold=False):
    """Formats a float into a LaTeX scientific notation string."""
    if val == 0:
        return r"$\mathbf{0}$" if bold else "$0$"
    s = f"{val:.2e}"
    base, exp = s.split('e')
    exp = int(exp) # Convert to int to remove leading zeros (e.g., +016 -> 16)
    
    if bold:
        return f"$\\mathbf{{{base} \\times 10^{{{exp}}}}}$"
    return f"${base} \\times 10^{{{exp}}}$"

def parse_val(v):
    """Parses scientific notation, fixing missing 'e' (e.g., 5.65+16 -> 5.65e+16)."""
    v = v.strip()
    if '+' in v and 'e' not in v.lower():
        v = v.replace('+', 'e+')
    return float(v)

def main():
    # Configuration
    plot_dir = "comparison_plots"
    tex_filename = "CrossSection_Comparison.tex"
    csv_filename = "POT_list.csv"
    
    # Define your xF bin edges again to match the slide titles with the plot titles
    xf_edges = [
        0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
        0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80
    ]

    # Process the CSV Data
    pot_data = []
    totals = {'LH2': 0.0, 'Flask': 0.0, 'LD2': 0.0, 'Total': 0.0}
    
    with open(csv_filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            rs = row['Roadset']
            
            lh2 = parse_val(row['LH2'])
            ld2 = parse_val(row['LD2'])
            flask = parse_val(row['Flask'])
            row_total = lh2 + ld2 + flask
            
            # Add to combined totals
            totals['LH2'] += lh2
            totals['LD2'] += ld2
            totals['Flask'] += flask
            totals['Total'] += row_total
            
            pot_data.append({
                'RS': rs,
                'LH2': lh2,
                'Flask': flask,
                'LD2': ld2,
                'Total': row_total
            })

    # Start writing the LaTeX content
    tex_content = []
    
    # 1. Preamble
    tex_content.append(r"\documentclass[aspectratio=169]{beamer}")
    tex_content.append(r"\usetheme{Madrid}")
    tex_content.append(r"\usecolortheme{whale}")
    tex_content.append(r"\usepackage{graphicx}")
    tex_content.append(r"\usepackage{booktabs} % For better looking tables")
    tex_content.append(r"\title{Cross-Section Comparison: RS67 vs RS57-70}")
    tex_content.append(r"\author{Analysis Report}")
    tex_content.append(r"\date{\today}")
    
    tex_content.append(r"\begin{document}")
    
    # 2. Title Slide
    tex_content.append(r"\begin{frame}")
    tex_content.append(r"  \titlepage")
    tex_content.append(r"\end{frame}")

    # --- NEW SLIDE: POT Table ---
    tex_content.append(r"\begin{frame}{Live Protons on Target (POT) Summary}")
    tex_content.append(r"  \centering")
    tex_content.append(r"  \small")
    tex_content.append(r"  \begin{tabular}{lcccc}")
    tex_content.append(r"    \toprule")
    tex_content.append(r"    \textbf{Roadset} & \textbf{LH2 POT} & \textbf{Flask POT} & \textbf{LD2 POT} & \textbf{Total (H2+D2+F)} \\")
    tex_content.append(r"    \midrule")
    
    # Loop through CSV data for rows
    for d in pot_data:
        lh2_str = format_tex_sci(d['LH2'])
        flask_str = format_tex_sci(d['Flask'])
        ld2_str = format_tex_sci(d['LD2'])
        total_str = format_tex_sci(d['Total'])
        
        tex_content.append(f"    {d['RS']} & {lh2_str} & {flask_str} & {ld2_str} & {total_str} \\\\")

    tex_content.append(r"    \midrule")
    
    # Combined Totals Row
    lh2_tot = format_tex_sci(totals['LH2'], bold=True)
    flask_tot = format_tex_sci(totals['Flask'], bold=True)
    ld2_tot = format_tex_sci(totals['LD2'], bold=True)
    comb_tot = format_tex_sci(totals['Total'], bold=True)
    
    tex_content.append(f"    \\textbf{{57--70 Combined}} & {lh2_tot} & {flask_tot} & {ld2_tot} & {comb_tot} \\\\")
    tex_content.append(r"    \bottomrule")
    tex_content.append(r"  \end{tabular}")
    tex_content.append(r"  \vspace{0.5cm}")
    tex_content.append(r"  \begin{itemize}")
    tex_content.append(r"    \item POT values extracted from analysis header file.")
    tex_content.append(r"    \item RS57--70 used as the primary integrated dataset for comparison.")
    tex_content.append(r"  \end{itemize}")
    tex_content.append(r"\end{frame}")

    # 3. Loop through xF bins (0 to 15)
    for i in range(16):
        xf_min = xf_edges[i]
        xf_max = xf_edges[i+1]
        
        # Slide Header
        tex_content.append(r"\begin{frame}{Comparison for $" + str(xf_min) + r" \leq x_F < " + str(xf_max) + r"$}")
        
        # Use columns to put LH2 and LD2 side-by-side
        tex_content.append(r"  \begin{columns}")
        
        # --- Left Column: LH2 ---
        lh2_path = os.path.join(plot_dir, f"XSec_Compare_LH2_xF_{i}.pdf")
        tex_content.append(r"    \begin{column}{0.5\textwidth}")
        if os.path.exists(lh2_path):
            tex_content.append(r"      \centering \textbf{Target: LH2}")
            tex_content.append(r"      \includegraphics[width=\textwidth]{" + lh2_path + r"}")
        else:
            tex_content.append(r"      \centering LH2 Plot missing for bin " + str(i))
        tex_content.append(r"    \end{column}")

        # --- Right Column: LD2 ---
        ld2_path = os.path.join(plot_dir, f"XSec_Compare_LD2_xF_{i}.pdf")
        tex_content.append(r"    \begin{column}{0.5\textwidth}")
        if os.path.exists(ld2_path):
            tex_content.append(r"      \centering \textbf{Target: LD2}")
            tex_content.append(r"      \includegraphics[width=\textwidth]{" + ld2_path + r"}")
        else:
            tex_content.append(r"      \centering LD2 Plot missing for bin " + str(i))
        tex_content.append(r"    \end{column}")
        
        tex_content.append(r"  \end{columns}")
        tex_content.append(r"\end{frame}")

    # 4. End Document
    tex_content.append(r"\end{document}")

    # Write to file
    with open(tex_filename, "w") as f:
        f.write("\n".join(tex_content))

    print(f"LaTeX file '{tex_filename}' has been generated.")
    print("To compile it, run: pdflatex " + tex_filename)

if __name__ == "__main__":
    main()