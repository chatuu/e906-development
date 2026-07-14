"""
generate_beamer.py
Generates a LaTeX Beamer presentation for the Roadset Systematic Procedure,
dynamically pulling table values from CSV files and omitting the final pT bin.
"""

import os
import csv

# Define absolute paths to the generated CSV files
LH2_CSV_PATH = "/root/github/e906-development/src/xsec_pT/RS57-70_weighted_average/Roadset_Sys_StdDev_LH2_geom.csv"
LD2_CSV_PATH = "/root/github/e906-development/src/xsec_pT/RS57-70_weighted_average/Roadset_Sys_StdDev_LD2_geom.csv"

# Global pT Bins to map bin_idx back to the kinematic range
PT_BINS = [0.0, 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8, 2.5]

def generate_table_rows(csv_filepath):
    """Reads the CSV file and returns formatted LaTeX table rows, excluding the final bin."""
    if not os.path.exists(csv_filepath):
        print(f"Warning: Could not find {csv_filepath}")
        return r"        \multicolumn{5}{c}{\textit{Data file not found. Run variance calculation first.}} \\"
    
    latex_rows = []
    with open(csv_filepath, "r") as f:
        reader = csv.reader(f)
        
        # Skip the first two header rows
        next(reader, None) 
        next(reader, None) 
        
        for row in reader:
            if not row or len(row) < 4: 
                continue
                
            try:
                bin_idx = int(row[0])
                mean_val = float(row[1])
                std_dev = float(row[2])
                rel_err = float(row[3])
                
                pt_min = PT_BINS[bin_idx]
                pt_max = PT_BINS[bin_idx + 1]
                
                # Format into a valid LaTeX table row
                latex_row = f"        {bin_idx} & $[{pt_min:.2f}, {pt_max:.2f})$ & {mean_val:.6f} & {std_dev:.6f} & {rel_err:.2f}\\% \\\\"
                latex_rows.append(latex_row)
            except ValueError:
                # Catch any unexpected parsing errors gracefully
                continue
                
    # --- Omitting the final pT bin row ---
    if latex_rows:
        latex_rows = latex_rows[:-1]
            
    return "\n".join(latex_rows)

def generate_tex_file(filename="roadset_systematics.tex"):
    
    # 1. Dynamically read CSV data
    lh2_table_body = generate_table_rows(LH2_CSV_PATH)
    ld2_table_body = generate_table_rows(LD2_CSV_PATH)
    
    # 2. LaTeX Template (using a raw string)
    latex_template = r'''\documentclass[aspectratio=169]{beamer}
\usetheme{Madrid}
\usecolortheme{default}

\usepackage{amsmath}
\usepackage{booktabs}
\usepackage{graphicx}

% Title Information
\title[Roadset Systematics]{Evaluation of Tracking Systematic Uncertainties via Roadset Variance}
\subtitle{Step-by-Step Methodology for the Drell-Yan Absolute Cross-Section}
\author{Chatura Kuruppu}
\institute[Fermilab / NMSU]{New Mexico State University \\ Fermilab SpinQuest/SeaQuest}
\date{July 20, 2026}

\begin{document}

% Slide 1: Title
\begin{frame}
    \titlepage
\end{frame}

% Slide 2: Motivation
\begin{frame}{Motivation: Tracking Systematics}
    \textbf{Objective:} Quantify the systematic uncertainty introduced by different trigger road definitions (RS57, RS59, RS62, RS67, RS70).
    \vspace{0.5cm}
    \begin{itemize}
        \item \textbf{Challenge:} Different roadsets have drastically different Protons on Target (POT) and differing trigger acceptances.
        \item \textbf{Approach:} Calculate the single differential cross-section independently for each roadset, then evaluate the bin-by-bin variance.
        \item \textbf{Constraint:} Because the reliability of high-statistics roadsets (e.g., RS67) is under evaluation, we must treat roadsets equally to prevent biased pulls.
    \end{itemize}
\end{frame}

% Slide 3: Step 1 & 2
\begin{frame}{Step 1 \& 2: Independent Execution \& Dynamic Normalization}
    To prevent normalization distortions (the "Global POT Trap"), the analysis pipeline was refactored to execute \textit{per roadset}. 
    \vspace{0.3cm}
    
    Each execution dynamically overrides the global POT with the roadset-specific POT to yield the true independent cross-section:
    \vspace{0.3cm}
    
    \begin{table}[]
        \centering
        \begin{tabular}{@{}lccc@{}}
        \toprule
        \textbf{Roadset} & \textbf{LH2 POT} & \textbf{LD2 POT} & \textbf{Flask POT} \\ \midrule
        \textbf{RS57} & $3.53 \times 10^{16}$ & $1.76 \times 10^{16}$ & $3.92 \times 10^{15}$ \\
        \textbf{RS59} & $9.37 \times 10^{15}$ & $4.32 \times 10^{15}$ & $1.01 \times 10^{15}$ \\
        \textbf{RS62} & $5.28 \times 10^{16}$ & $2.38 \times 10^{16}$ & $1.10 \times 10^{16}$ \\
        \textbf{RS67} & $1.61 \times 10^{17}$ & $7.69 \times 10^{16}$ & $3.66 \times 10^{16}$ \\
        \textbf{RS70} & $1.79 \times 10^{16}$ & $8.75 \times 10^{15}$ & $3.84 \times 10^{15}$ \\ \bottomrule
        \end{tabular}
    \end{table}
\end{frame}

% Slide 4: Step 3 - Mathematics
\begin{frame}{Step 3: Variance Calculation Mathematics}
    For a given $p_T$ bin, let $x_i$ be the absolute cross-section for roadset $i$, and $e_i$ be its statistical error.
    \vspace{0.5cm}
    
    \textbf{Inverse-Variance Weighting:} \\
    To suppress massive Poisson fluctuations from low-POT datasets (e.g., RS59), we scale by statistical power ($w_i = 1/e_i^2$):
    \vspace{0.3cm}
    \[
    \mu_w = \frac{\sum w_i x_i}{\sum w_i}
    \]
    \vspace{0.3cm}
    \[
    \sigma_{sys, w} = \sqrt{ \frac{\sum w_i}{\left(\sum w_i\right)^2 - \sum w_i^2} \sum w_i (x_i - \mu_w)^2 }
    \]
\end{frame}

% Slide 5: Results Table (LH2)
\begin{frame}{Step 4: Resulting Uncertainties (LH2)}
    Using the inverse-variance weighted strategy suppresses low-statistics noise, yielding the following systematic errors for Liquid Hydrogen.
    \vspace{0.3cm}
    
    \begin{table}[]
        \centering
        \begin{tabular}{@{}ccccc@{}}
        \toprule
        \textbf{$p_T$ Bin} & \textbf{Range [GeV]} & \textbf{Mean $\sigma_w$} & \textbf{Std Dev $\sigma_{sys,w}$} & \textbf{Relative Error} \\ \midrule
__LH2_TABLE_BODY__
        \bottomrule
        \end{tabular}
    \end{table}
\end{frame}

% Slide 6: Results Table (LD2)
\begin{frame}{Step 4: Resulting Uncertainties (LD2)}
    Applying the same inverse-variance weighted procedure to Liquid Deuterium targets yields comparable systematic margins.
    \vspace{0.3cm}
    
    \begin{table}[]
        \centering
        \begin{tabular}{@{}ccccc@{}}
        \toprule
        \textbf{$p_T$ Bin} & \textbf{Range [GeV]} & \textbf{Mean $\sigma_w$} & \textbf{Std Dev $\sigma_{sys,w}$} & \textbf{Relative Error} \\ \midrule
__LD2_TABLE_BODY__
        \bottomrule
        \end{tabular}
    \end{table}
\end{frame}

% Slide 7: Conclusions
\begin{frame}{Discussion \& Next Steps}
    \textbf{Why is the relative error still $\sim$13--18\%?}
    \vspace{0.3cm}
    \begin{enumerate}
        \item \textbf{Statistical Noise Limitations:} While inverse-variance weighting heavily suppresses the noise of low-POT datasets (RS59, RS70), it cannot completely eliminate the pull from extreme outliers.
        \item \textbf{Global Acceptance Mismatch:} A single \texttt{acceptance\_mass\_xF.root} file is currently used across all roadsets. Because different triggers have distinct kinematic acceptances, dividing by a global average creates artificial variance.
    \end{enumerate}
    \vspace{0.4cm}
    \textbf{Path Forward:} Implement roadset-specific Monte Carlo acceptance maps and evaluate the necessity of excluding extremely low-POT datasets from the overall variance calculation.
\end{frame}

\end{document}
'''
    
    # 3. Inject the dynamic table rows into the LaTeX template
    final_latex = latex_template.replace("__LH2_TABLE_BODY__", lh2_table_body)
    final_latex = final_latex.replace("__LD2_TABLE_BODY__", ld2_table_body)

    # 4. Write to disk
    with open(filename, "w") as f:
        f.write(final_latex)
        
    print(f"[*] Successfully generated {filename} using dynamic CSV data.")

if __name__ == "__main__":
    generate_tex_file()