import os

def main():
    # Configuration
    plot_dir = "comparison_plots"
    tex_filename = "CrossSection_Comparison.tex"
    
    # Define your xF bin edges again to match the slide titles with the plot titles
    xf_edges = [
        0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
        0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80
    ]

    # Start writing the LaTeX content
    tex_content = []
    
    # 1. Preamble
    tex_content.append(r"\documentclass[aspectratio=169]{beamer}")
    tex_content.append(r"\usetheme{Madrid}")
    tex_content.append(r"\usecolortheme{whale}")
    tex_content.append(r"\usepackage{graphicx}")
    tex_content.append(r"\usepackage{booktabs} % For better looking tables")
    tex_content.append(r"\title{Cross-Section Comparison: RS67 vs RS57}")
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
    tex_content.append(r"    57 & $3.79 \times 10^{16}$ & $4.23 \times 10^{15}$ & $1.90 \times 10^{16}$ & $6.11 \times 10^{16}$ \\")
    tex_content.append(r"    59 & $8.44 \times 10^{15}$ & $9.18 \times 10^{14}$ & $3.84 \times 10^{15}$ & $1.32 \times 10^{16}$ \\")
    tex_content.append(r"    62 & $5.51 \times 10^{16}$ & $1.13 \times 10^{16}$ & $2.55 \times 10^{16}$ & $9.19 \times 10^{16}$ \\")
    tex_content.append(r"    67 & $1.57 \times 10^{17}$ & $3.58 \times 10^{16}$ & $7.51 \times 10^{16}$ & $2.68 \times 10^{17}$ \\")
    tex_content.append(r"    70 & $1.76 \times 10^{16}$ & $3.71 \times 10^{15}$ & $8.54 \times 10^{15}$ & $2.99 \times 10^{16}$ \\")
    tex_content.append(r"    \midrule")
    tex_content.append(r"    \textbf{57--70 Combined} & $\mathbf{2.76 \times 10^{17}}$ & $\mathbf{5.59 \times 10^{16}}$ & $\mathbf{1.32 \times 10^{17}}$ & $\mathbf{4.64 \times 10^{17}}$ \\")
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