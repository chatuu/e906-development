#!/usr/bin/env python3
import csv
import os
import numpy as np
from collections import defaultdict

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

\title{Weekly Analysis Update: August 14, 2026}
\author{Chatura Kuruppu}
\institute{New Mexico State University\\ \vspace{0.1cm} SeaQuest Experiment (E906)}
\date{\today}

\begin{document}

% --- TITLE SLIDE ---
\begin{frame}
    \titlepage
\end{frame}

% --- OVERVIEW SLIDE ---
\begin{frame}{Overview}
    \begin{itemize}
        \setlength\itemsep{1em}
        \item Review of comments from the recent E906 presentation
        \item Updated cross-section results
        \item Evaluation of unfolding systematic uncertainties
    \end{itemize}
\end{frame}

% --- COMMENTS SLIDE (PART 1) ---
\begin{frame}{Review of Comments on Unfolded Cross-Section Results (Part 1)}
    \begin{enumerate}
        \item \textbf{Were Messy/Clean corrections applied during unfolding?} \\
        Messy/Clean corrections are accounted for during the acceptance correction phase, which is handled separately from the unfolding process.
        
        \vspace{0.3cm}
        \item \textbf{How should the large standard deviations in the ratio above 7 GeV be interpreted?} \\
        Error bars were calculated by rigorously propagating uncertainties, including fully correlated systematic effects (which were previously underestimated).
        
        \vspace{0.3cm}
        \item \textbf{We need to define the systematic uncertainty introduced by unfolding. (Proposed method: calculating the systematic uncertainty by varying the number of iterations by one step.)} \\
        \textit{Why is Clean unfolding preferred over Messy unfolding?} \\
        Using Messy MC for unfolding risks double-counting corrections, as the acceptance correction already incorporates Messy MC effects.
    \end{enumerate}
\end{frame}

% --- COMMENTS SLIDE (PART 2) ---
\begin{frame}{Review of Comments on Unfolded Cross-Section Results (Part 2)}
    \begin{enumerate}
        \setcounter{enumi}{3}
        \item \textbf{Is it possible to use a different number of iterations for Messy versus Clean unfolding?} \\
        Currently, we determine the optimal number of iterations by minimizing the divergence between the mean statistical error and the $\Delta\chi^2$ convergence curves.
        
        \vspace{0.3cm}
        \item \textbf{To determine unfolding systematics, could we adopt the closure test methodology used by Harsha?} \\
        This requires further discussion with the working group.
        
        \vspace{0.3cm}
        \item \textbf{Combine overflow and underflow bins, and investigate the required bin widths (they should exceed the detector resolution).} \\
        We plan to broaden the peripheral bins and evaluate the impact on the final cross-section results.
    \end{enumerate}
\end{frame}

% --- SECTION TRANSITION SLIDE ---
\begin{frame}[plain,c]
    \begin{center}
        \Huge \textbf{Updated Cross-Section Results}
    \end{center}
\end{frame}
""")

        # Generate comparative Plot slides for LH2 and LD2
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
\begin{frame}{""" + tgt + f" Target: {clean_range}" + r"""}
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

        # --- SYSTEMATICS EQUATION SLIDE ---
        f.write(r"""
\begin{frame}{Evaluation of Systematic Uncertainty Due to Unfolding}
    We determine the systematic uncertainty introduced by the unfolding algorithm for both the Messy and Clean configurations.
    
    \vspace{0.5cm}
    \begin{block}{Systematic Uncertainty Calculation}
        \begin{equation*}
            \text{Systematic Uncertainty (\%)} = \left( \frac{| \sigma_{\text{unf}}^{N=4} - \sigma_{\text{unf}}^{N=3} |}{\sigma_{\text{unf}}^{N=3}} \right) \times 100
        \end{equation*}
    \end{block}
    \vspace{0.2cm}
    \textit{Note: The deviation between consecutive algorithmic iterations provides a robust estimate of the structural bias introduced by the unfolding matrix.}
\end{frame}
""")

        # Generate Systematics Table Slides (Organized by xF bin)
        f.write("\n\\section{Systematics Tables}\n")

        try:
            if os.path.exists(csv_path):
                # Read all rows into a nested dictionary grouped by xF bin and Target
                # Structure: data_by_xf["-0.05 <= xF < 0.00"]["LH2"] = [row1, row2...]
                data_by_xf = defaultdict(lambda: {"LH2": [], "LD2": []})
                
                with open(csv_path, 'r') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        xf_str = row['xF bin']
                        tgt_str = row['Target']
                        data_by_xf[xf_str][tgt_str].append(row)
                
                # Iterate through the extracted xF bins and create one slide per xF bin
                for xf_str, target_data in data_by_xf.items():
                    # Format the xF string for LaTeX Math Mode
                    xf_math = xf_str.replace("<=", r"\leq").replace("<", r"<").replace("xF", "x_F")
                    
                    f.write(r"""
\begin{frame}{Unfolding Systematics: $""" + xf_math + r"""$}
    \begin{columns}[t]
""")
                    # Create a side-by-side table for LH2 and LD2
                    for tgt in ["LH2", "LD2"]:
                        f.write(r"""
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{""" + tgt + r""" Target} \vspace{0.1cm} \\
            \resizebox{0.95\textwidth}{!}{
            \begin{tabular}{lcc}
                \toprule
                \textbf{Mass Bin} & \textbf{Syst. Error (Clean) [\%]} & \textbf{Syst. Error (Messy) [\%]} \\
                \midrule
""")
                        # Inject the rows for this specific target and xF bin
                        if target_data[tgt]:
                            for row in target_data[tgt]:
                                m_math = row['Mass bin'].replace("<=", r"\leq").replace("<", r"<").replace("Mass", r"\text{Mass}")
                                syst_clean = row[r'Unfolding syst (clean) [%]']
                                syst_messy = row[r'Unfolding syst (messy) [%]']
                                f.write(f"                ${m_math}$ & {syst_clean} & {syst_messy} \\\\\n")
                        else:
                            f.write(r"                \multicolumn{3}{c}{\textit{No data available for this bin}} \\" + "\n")
                        
                        f.write(r"""                \bottomrule
            \end{tabular}
            }
        \end{column}
""")
                    f.write(r"""    \end{columns}
\end{frame}
""")
            else:
                # Fallback if CSV is missing
                f.write(r"""
\begin{frame}{Unfolding Systematics Summary}
    \centering
    CSV file not found at the specified path.
\end{frame}
""")
        except Exception as e:
            f.write(r"""
\begin{frame}{Unfolding Systematics Summary}
    \centering
    Error loading CSV: """ + str(e) + r"""
\end{frame}
""")

        # --- CONCLUSIONS SLIDE ---
        f.write(r"""
\begin{frame}{Conclusions and Next Steps}
    \begin{itemize}
        \setlength\itemsep{1em}
        \item Successfully implemented multi-dimensional unfolding for the double differential cross-section.
        \item Recalculation of acceptance corrections for Runs 2 and 3 independently is currently underway (Kenichi and Harsha).
        \item Finalize the publication draft and initiate the collaboration review process.
    \end{itemize}
\end{frame}
""")

        f.write(r"\end{document}" + "\n")
        
    print(f"Successfully generated {tex_filename}")

if __name__ == "__main__":
    main()