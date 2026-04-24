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
    #tex_content.append(r"\usepackage{siunitx}")
    tex_content.append(r"\title{Cross-Section Comparison: RS67 vs RS57-70}")
    tex_content.append(r"\author{Analysis Report}")
    tex_content.append(r"\date{\today}")
    
    tex_content.append(r"\begin{document}")
    
    # 2. Title Slide
    tex_content.append(r"\begin{frame}")
    tex_content.append(r"  \titlepage")
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