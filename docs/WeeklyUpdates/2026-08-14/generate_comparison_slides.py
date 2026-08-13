#!/usr/bin/env python3
import csv
import os

def main():
    # --- Configuration ---
    targets = ["LH2", "LD2"]
    xf_bins = [
        "-0.05_0.00", "0.00_0.05", "0.05_0.10", "0.10_0.15", "0.15_0.20", 
        "0.20_0.25", "0.25_0.30", "0.30_0.35", "0.35_0.40", "0.40_0.45", 
        "0.45_0.50", "0.50_0.55", "0.55_0.60", "0.60_0.65", "0.65_0.70", 
        "0.70_0.75", "0.75_0.80", "0.80_0.85"
    ]
    
    # Using absolute paths as requested
    path_iter3 = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70_with_unfolding/num_iter_3/"
    path_iter4 = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70_with_unfolding/num_iter_4/"
    csv_path = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70_with_unfolding/Unfolding_Systematics_Comparison.csv"
    
    tex_filename = "Unfolding_Comparison.tex"

    with open(tex_filename, "w") as f:
        # Preamble
        f.write(r"""\documentclass[aspectratio=169]{beamer}
\usetheme{Madrid}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{amsmath}

\title{Unfolding Stability Study}
\subtitle{Iteration 3 vs. Iteration 4 Comparison}
\author{Chatura Kuruppu}
\date{\today}

\begin{document}
\begin{frame}
    \titlepage
\end{frame}
""")

        # Generate comparative slides for LH2 and LD2
        for tgt in targets:
            f.write(f"\\section{{{tgt} Comparison}}\n")
            for bin_range in xf_bins:
                # Safely split and format with math mode ($)
                parts = bin_range.split("_")
                if len(parts) == 2:
                    clean_range = f"${parts[0]} \\leq x_F < {parts[1]}$"
                else:
                    clean_range = bin_range
                    
                f.write(r"""
\begin{frame}{""" + tgt + f": {clean_range}" + r"""}
    \begin{columns}
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Iterations = 3}\\
            \includegraphics[width=\textwidth,height=0.7\textheight,keepaspectratio]{""" + path_iter3 + f"CrossSection_{tgt}_xF_{bin_range}_GeoCenter.pdf" + r"""}
        \end{column}
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Iterations = 4}\\
            \includegraphics[width=\textwidth,height=0.7\textheight,keepaspectratio]{""" + path_iter4 + f"CrossSection_{tgt}_xF_{bin_range}_GeoCenter.pdf" + r"""}
        \end{column}
    \end{columns}
\end{frame}
""")

        # Generate Systematics Table Slides (Paginated)
        f.write("\n\\section{Systematics Table}\n")

        try:
            if os.path.exists(csv_path):
                # Read all rows into memory
                all_rows = []
                with open(csv_path, 'r') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        all_rows.append(row)
                
                # Define how many rows fit on one slide
                rows_per_slide = 12
                total_chunks = (len(all_rows) + rows_per_slide - 1) // rows_per_slide
                
                for i in range(total_chunks):
                    chunk = all_rows[i * rows_per_slide : (i + 1) * rows_per_slide]
                    
                    f.write(r"""
\begin{frame}{Unfolding Systematics Summary (Iter 3 vs 4) - Part """ + str(i + 1) + r"""}
    \centering
    \resizebox{0.95\textwidth}{!}{
    \begin{tabular}{llccc}
        \toprule
        \textbf{Tgt} & \textbf{xF Bin} & \textbf{Mass Bin} & \textbf{Syst (Clean) [\%]} & \textbf{Syst (Messy) [\%]} \\
        \midrule
""")
                    
                    # Inject rows for this chunk
                    for row in chunk:
                        # Clean strings and format them into safe LaTeX math mode
                        xf = row['xF bin'].replace("<=", r"\leq").replace("<", r"<").replace("xF", "x_F")
                        m = row['Mass bin'].replace("<=", r"\leq").replace("<", r"<").replace("Mass", r"\text{Mass}")
                        
                        syst_clean = row[r'Unfolding syst (clean) [%]']
                        syst_messy = row[r'Unfolding syst (messy) [%]']
                        
                        f.write(f"        {row['Target']} & ${xf}$ & ${m}$ & {syst_clean} & {syst_messy} \\\\\n")
                    
                    f.write(r"""        \bottomrule
    \end{tabular}
    }
\end{frame}
""")
            else:
                # Fallback if CSV is missing
                f.write(r"""
\begin{frame}{Unfolding Systematics Summary}
    \centering
    CSV file not found.
\end{frame}
""")
        except Exception as e:
            f.write(r"""
\begin{frame}{Unfolding Systematics Summary}
    \centering
    Error loading CSV: """ + str(e) + r"""
\end{frame}
""")

        f.write(r"\end{document}" + "\n")
        
    print(f"Successfully generated {tex_filename}")

if __name__ == "__main__":
    main()