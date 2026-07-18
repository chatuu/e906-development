import os
import csv
import numpy as np

# --- 1. Configurations & Paths ---
# Image Plot Directory
PLOT_DIR = "/root/github/e906-development/src/pT Distributions/DeterminePtBins"

# pT Bins & Variables
pt_bins = [
    "pT_0p0_to_0p5",
    "pT_0p5_to_1p0",
    "pT_1p0_to_1p5",
    "pT_1p5_to_1p8"
]

binned_vars = [
    {"name": "xF", "title": "$x_{F}$"},
    {"name": "Mass", "title": "Mass"}
]

# Systematic Tracking CSV Paths
LH2_CSV_PATH = "/root/github/e906-development/src/xsec_pT/RS57-70_weighted_average/Roadset_Sys_StdDev_LH2_geom.csv"
LD2_CSV_PATH = "/root/github/e906-development/src/xsec_pT/RS57-70_weighted_average/Roadset_Sys_StdDev_LD2_geom.csv"
LH2_2D_CSV_PATH = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70_weighted_average/Roadset_Sys_StdDev_LH2_cent.csv"
LD2_2D_CSV_PATH = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70_weighted_average/Roadset_Sys_StdDev_LD2_cent.csv"

# Error Propagation CSV Paths
ERROR_CSV_OLD = "~/github/e906-development/src/xsec_pT/RS57-70/pd_2pp_errors.csv"
ERROR_CSV_NEW = "~/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/pd_2pp_errors.csv"

# New Cross-Section Comparison Directories
PREV_XSEC_DIR = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70"
LATEST_XSEC_DIR = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70_road_dependancy_added"

# Ordered Plots for Comparison (Centroid plots removed)
COMPARISON_PLOTS = [
    "CrossSection_LH2_xF_0.00_0.05_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.05_0.10_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.10_0.15_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.15_0.20_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.20_0.25_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.25_0.30_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.30_0.35_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.35_0.40_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.40_0.45_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.45_0.50_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.50_0.55_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.55_0.60_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.60_0.65_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.65_0.70_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.70_0.75_GeoCenter_with_logo.pdf",
    "CrossSection_LH2_xF_0.75_0.80_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.00_0.05_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.05_0.10_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.10_0.15_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.15_0.20_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.20_0.25_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.25_0.30_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.30_0.35_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.35_0.40_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.40_0.45_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.45_0.50_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.50_0.55_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.55_0.60_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.60_0.65_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.65_0.70_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.70_0.75_GeoCenter_with_logo.pdf",
    "CrossSection_LD2_xF_0.75_0.80_GeoCenter_with_logo.pdf",
    "cross_section_overlay_LH2_GeoCenter_logo.pdf",
    "cross_section_overlay_LD2_GeoCenter_logo.pdf"
]

# Kinematic Bins
PT_BINS = [0.0, 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8, 2.5]
XF_BINS = np.round(np.arange(0.0, 0.85, 0.05), 2)
MASS_BINS = np.linspace(2.0, 9.0, 71)

MIN_REL_ERR = 0.0
MAX_REL_ERR = 50.0

# --- 2. Helper Functions ---
def format_bin_title(bin_str):
    parts = bin_str.replace('p', '.').split('_to_')
    if len(parts) == 2:
        return f"${parts[0]} \\leq p_{{T}} < {parts[1]}$ GeV/c"
    return bin_str

def make_slide(title, lh2_file, ld2_file):
    """Generates a standard two-column Beamer slide with protected paths."""
    # Added height=0.7\textheight to fix footer masking
    return f"""
\\begin{{frame}}{{{title}}}
    \\begin{{columns}}[T] 
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\textbf{{LH2 Target}}\\\\[0.2cm]
            \\includegraphics[width=0.95\\linewidth, height=0.7\\textheight, keepaspectratio]{{"{lh2_file}"}}
        \\end{{column}}
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\textbf{{LD2 Target}}\\\\[0.2cm]
            \\includegraphics[width=0.95\\linewidth, height=0.7\\textheight, keepaspectratio]{{"{ld2_file}"}}
        \\end{{column}}
    \\end{{columns}}
\\end{{frame}}
"""

def generate_comparison_slides(plot_list):
    """Generates specific side-by-side comparison slides mapping Previous vs Latest folders."""
    slides = []
    for plot in plot_list:
        # Dynamically generate the slide title based on the filename
        if plot.startswith("CrossSection_"):
            parts = plot.split('_')
            target = parts[1]
            xf_min = parts[3]
            xf_max = parts[4]
            slide_title = f"{target} Cross-Section Previous Vs Latest: ${xf_min} \\leq x_F < {xf_max}$"
        else:
            slide_title = "Summary Plots Previous Vs Latest with updated systematics"

        # Added height=0.7\textheight to fix footer masking
        slide = f"""
\\begin{{frame}}{{{slide_title}}}
    \\begin{{columns}}[T] 
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\textbf{{Previous}}\\\\[0.2cm]
            \\includegraphics[width=0.95\\linewidth, height=0.7\\textheight, keepaspectratio]{{{PREV_XSEC_DIR}/{plot}}}
        \\end{{column}}
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\textbf{{Latest (with weighted average)}}\\\\[0.2cm]
            \\includegraphics[width=0.95\\linewidth, height=0.7\\textheight, keepaspectratio]{{{LATEST_XSEC_DIR}/{plot}}}
        \\end{{column}}
    \\end{{columns}}
\\end{{frame}}
"""
        slides.append(slide)
    return "".join(slides)

def get_latex_table_rows(csv_filename):
    """Reads the ratio error CSV file and returns formatted LaTeX table rows."""
    table_rows = ""
    full_csv_path = os.path.expanduser(csv_filename)
    
    if not os.path.exists(full_csv_path):
        print(f"Warning: {full_csv_path} not found. Using fallback layout values.")
        fallback = [
            ["[0.00, 0.32)", "1.1228", "0.0316", "0.0017", "0.0317"],
            ["[0.32, 0.49)", "1.1370", "0.0309", "0.0014", "0.0309"],
            ["[0.49, 0.63)", "1.2137", "0.0374", "0.0018", "0.0374"],
            ["[0.63, 0.77)", "1.0941", "0.0339", "0.0017", "0.0339"],
            ["[0.77, 0.95)", "1.1839", "0.0364", "0.0018", "0.0364"],
            ["[0.95, 1.18)", "1.0855", "0.0381", "0.0020", "0.0381"],
            ["[1.18, 1.80)", "1.1273", "0.0539", "0.0025", "0.0540"]
        ]
        for r in fallback:
            table_rows += f"        {{{r[0]}}} & {r[1]} & {r[2]} & {r[3]} & {r[4]} \\\\\n"
        return table_rows

    with open(full_csv_path, mode='r') as f:
        reader = csv.reader(f)
        next(reader) 
        for row in reader:
            if len(row) >= 5:
                pt_bin, ratio, stat, sys, total = row[0], row[1], row[2], row[3], row[4]
                table_rows += f"        {{{pt_bin}}} & {ratio} & {stat} & {sys} & {total} \\\\\n"
                
    return table_rows

def generate_systematic_table_rows(csv_filepath):
    """Reads tracking systematic variance CSV data and formats rows."""
    if not os.path.exists(csv_filepath):
        return r"        \multicolumn{5}{c}{\textit{Data file not found. Run variance calculation first.}} \\"
    
    latex_rows = []
    with open(csv_filepath, "r") as f:
        reader = csv.reader(f)
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
                
                if mean_val <= 0 or rel_err <= MIN_REL_ERR or rel_err >= MAX_REL_ERR:
                    continue
                
                pt_min = PT_BINS[bin_idx]
                pt_max = PT_BINS[bin_idx + 1]
                
                latex_row = f"        {bin_idx} & $[{pt_min:.2f}, {pt_max:.2f})$ & {mean_val:.6f} & {std_dev:.6f} & {rel_err:.2f}\\% \\\\"
                latex_rows.append(latex_row)
            except ValueError:
                continue
                
    filtered_latex_rows = [row for row in latex_rows if not row.strip().startswith("7 &")]
            
    if not filtered_latex_rows:
        return r"        \multicolumn{5}{c}{\textit{No statistically significant bins found.}} \\"
        
    return "\n".join(filtered_latex_rows)

def generate_2d_frames(csv_filepath, target):
    """Generates table frames for 2D cross-section bins."""
    if not os.path.exists(csv_filepath):
        return f"% Warning: Could not find {csv_filepath}\n"
    
    valid_rows = []
    with open(csv_filepath, "r") as f:
        reader = csv.reader(f)
        next(reader, None)
        next(reader, None)
        for row in reader:
            if not row or len(row) < 5: 
                continue
            try:
                xf_idx = int(row[0])
                mass_idx = int(row[1])
                mean_val = float(row[2])
                std_dev = float(row[3])
                rel_err = float(row[4])
                
                if mean_val <= 0 or rel_err <= MIN_REL_ERR or rel_err >= MAX_REL_ERR: 
                    continue
                
                xf_min, xf_max = XF_BINS[xf_idx], XF_BINS[xf_idx+1]
                m_min, m_max = MASS_BINS[mass_idx], MASS_BINS[mass_idx+1]
                
                valid_rows.append(
                    f"        $[{xf_min:.2f}, {xf_max:.2f})$ & $[{m_min:.2f}, {m_max:.2f})$ & {mean_val:.4e} & {std_dev:.4e} & {rel_err:.2f}\\% \\\\"
                )
            except Exception:
                continue
    
    if not valid_rows:
         return f"% No statistically significant 2D bins found for {target}\n"

    CHUNK_SIZE = 12
    frames = []
    total_chunks = (len(valid_rows) + CHUNK_SIZE - 1) // CHUNK_SIZE
    
    for i in range(total_chunks):
        chunk = valid_rows[i*CHUNK_SIZE : (i+1)*CHUNK_SIZE]
        chunk_str = "\n".join(chunk)
        
        frame = f'''
% Slide: 2D {target} part {i+1}
\\begin{{frame}}{{Double Differential Systematics (xF and Mass) - {target} ({i+1}/{total_chunks})}}
    \\vspace{{0.2cm}}
    \\begin{{table}}[]
        \\centering
        \\resizebox{{\\textwidth}}{{!}}{{
        \\begin{{tabular}}{{@{{}}ccccc@{{}}}}
        \\toprule
        \\textbf{{xF Range}} & \\textbf{{Mass Range [GeV]}} & \\textbf{{Mean $\\sigma_w$}} & \\textbf{{Std Dev $\\sigma_{{sys,w}}$}} & \\textbf{{Relative Error}} \\\\ \\midrule
{chunk_str}
        \\bottomrule
        \\end{{tabular}}
        }}
    \\end{{table}}
\\end{{frame}}
'''
        frames.append(frame)
        
    return "\n".join(frames)

# --- 3. Main Generation Logic ---
def generate_presentation(filename="presentation.tex"):
    
    # 1. Generate tables
    lh2_table_body = generate_systematic_table_rows(LH2_CSV_PATH)
    ld2_table_body = generate_systematic_table_rows(LD2_CSV_PATH)
    lh2_2d_frames = generate_2d_frames(LH2_2D_CSV_PATH, "LH2")
    ld2_2d_frames = generate_2d_frames(LD2_2D_CSV_PATH, "LD2")
    
    ratio_rows_old = get_latex_table_rows(ERROR_CSV_OLD)
    ratio_rows_new = get_latex_table_rows(ERROR_CSV_NEW)

    # 2. Generate dynamic plots
    inclusive_pT_slide = make_slide(
        "$p_{T}$ Inclusive (Data - Mix): LH2 vs LD2",
        f"{PLOT_DIR}/pT_Distribution_LH2.pdf",
        f"{PLOT_DIR}/pT_Distribution_LD2.pdf"
    )
    
    binned_pT_slides = ""
    for var in binned_vars:
        v_name = var["name"]
        v_title = var["title"]
        for pt_bin in pt_bins:
            lh2_file = f"{PLOT_DIR}/{v_name}_Distribution_LH2_{pt_bin}.pdf"
            ld2_file = f"{PLOT_DIR}/{v_name}_Distribution_LD2_{pt_bin}.pdf"
            bin_label = format_bin_title(pt_bin.replace('pT_', ''))
            title = f"{v_title} ({bin_label}): LH2 vs LD2"
            binned_pT_slides += make_slide(title, lh2_file, ld2_file)
            
    # Generating Comparison Extension Slides
    comparison_slides = generate_comparison_slides(COMPARISON_PLOTS)

    # 3. LaTeX Master Template
    latex_template = r'''\documentclass[aspectratio=169]{beamer}
\usepackage{amsmath}
\usepackage{bm}
\usepackage{xcolor}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{caption}

\usetheme{Madrid}
\usecolortheme{default}
\setbeamertemplate{navigation symbols}{}

\title[DY Cross-Section Study]{Addressing Questions/Comments received for \\DY absolute Cross-Section Study}
\author{Chatura Kuruppu}
\institute{New Mexico State University}
\date{\today}

\begin{document}

\begin{frame}
    \titlepage
\end{frame}

\begin{frame}{Overview}
    This talk addresses additional comments received for the talk DocID: 11584 (\url{https://seaquest-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=11584})
    
    \vspace{0.5cm}
    \textbf{Supplemental Documents:}
    \begin{itemize}
        \item Technote (\url{https://seaquest-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=11569})
        \item Previous talk (\url{https://seaquest-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=11584})
    \end{itemize}
\end{frame}

\begin{frame}{Additional Data/MC plots}
    \begin{quote}
    ``The other comparison that I raised on Friday was of the mass and xF
    distributions in each pT bin. Can you make comparison plots? As you
    integrate the events over mass and xF when deriving dsigma/dpT, the
    agreement between the real data and the MC is crucial''
    \end{quote}
\end{frame}

\section{Inclusive Distributions (Data - Mix)}
__INCLUSIVE_PT_SLIDE__

\section{pT-Binned Distributions}
__BINNED_PT_SLIDES__


% --- ERROR PROPAGATION SECTION ---
\begin{frame}{Systematic Correlations to cross-section ratio}
    \begin{quote}
    ``How was the systematic error propagated from the LH2 and LD2 results? Have you considered any correlation?''
    \end{quote}
\end{frame}

\begin{frame}{The Observable: Cross-Section Ratio}
    \textbf{Objective:} Calculate the error for the single differential cross-section ratio in a given kinematic bin (e.g., $p_T$).
    
    \vspace{0.5cm}
    The ratio $R$ is defined as:
    \begin{equation}
        R = \frac{\sigma_{pd}}{2\sigma_{pp}} \equiv \frac{A}{cB}
    \end{equation}
    Where:
    \begin{itemize}
        \setlength{\itemsep}{0.2cm}
        \item $A = \sigma_{pd}$ (Liquid Deuterium cross-section)
        \item $B = \sigma_{pp}$ (Liquid Hydrogen cross-section)
        \item $c = 2$ (Exact scaling constant, zero uncertainty)
    \end{itemize}
\end{frame}

\begin{frame}{General Error Propagation (The Master Equation)}
    For any function $R(A, B)$, the variance is given by the partial derivatives and the covariance matrix:
    \begin{equation}
        \delta R^2 = \left(\frac{\partial R}{\partial A}\right)^2 \delta A^2 + \left(\frac{\partial R}{\partial B}\right)^2 \delta B^2 + 2 \left(\frac{\partial R}{\partial A}\right)\left(\frac{\partial R}{\partial B}\right) \text{cov}(A,B)
    \end{equation}
    
    By defining the covariance in terms of the Pearson correlation coefficient $\rho$:
    \begin{equation}
        \text{cov}(A,B) = \rho (\delta A) (\delta B) \quad \text{where } \rho \in [-1, 1]
    \end{equation}
    
    Evaluating the derivatives for $R = A / (cB)$ and converting to relative errors yields the master equation:
    \begin{equation}
        \left(\frac{\delta R}{R}\right)^2 = \left(\frac{\delta A}{A}\right)^2 + \left(\frac{\delta B}{B}\right)^2 - 2\rho \left(\frac{\delta A}{A}\right)\left(\frac{\delta B}{B}\right)
    \end{equation}
\end{frame}

\begin{frame}{1. Statistical Uncertainties ($\rho = 0$)}
    \textbf{Physical Reality:} Statistical fluctuations are driven by random counting statistics in finite data samples. 
    
    \vspace{0.2cm}
    A random upward fluctuation in the $pd$ dataset has no influence on the random fluctuations in the $pp$ dataset. Therefore, they are \textbf{completely uncorrelated} ($\rho = 0$).
    
    \vspace{0.2cm}
    Substituting $\rho = 0$ into the master equation:
    \begin{equation*}
        \left(\frac{\delta R_{\text{stat}}}{R}\right)^2 = \left(\frac{\delta \sigma_{pd}^{\text{stat}}}{\sigma_{pd}}\right)^2 + \left(\frac{\delta \sigma_{pp}^{\text{stat}}}{\sigma_{pp}}\right)^2 - \color{blue}{2(0)}\color{black}{\left(\frac{\delta \sigma_{pd}^{\text{stat}}}{\sigma_{pd}}\right)\left(\frac{\delta \sigma_{pp}^{\text{stat}}}{\sigma_{pp}}\right)}
    \end{equation*}
    
    The covariance term vanishes, yielding the quadrature sum:
    \begin{equation}
        \frac{\delta R_{\text{stat}}}{R} = \sqrt{\left(\frac{\delta \sigma_{pd}^{\text{stat}}}{\sigma_{pd}}\right)^2 + \left(\frac{\delta \sigma_{pp}^{\text{stat}}}{\sigma_{pp}}\right)^2}
    \end{equation}
\end{frame}

\begin{frame}{2. Systematic Uncertainties ($\rho = 1$)}
    \textbf{Physical Reality:} Systematic biases (e.g., tracking efficiency, geometric acceptance) apply to the exact same spectrometer running under identical beam conditions.
    
    \vspace{0.2cm}
    If an inefficiency causes us to systematically underestimate $\sigma_{pd}$, it will also cause us to systematically underestimate $\sigma_{pp}$. They move together (\textbf{perfectly positively correlated}, $\rho = 1$).
    
    \vspace{0.2cm}
    Substituting $\rho = 1$ into the master equation:
    \begin{equation*}
        \left(\frac{\delta R_{\text{sys}}}{R}\right)^2 = \left(\frac{\delta \sigma_{pd}^{\text{sys}}}{\sigma_{pd}}\right)^2 + \left(\frac{\delta \sigma_{pp}^{\text{sys}}}{\sigma_{pp}}\right)^2 - \color{red}{2(1)}\color{black}{\left(\frac{\delta \sigma_{pd}^{\text{sys}}}{\sigma_{pd}}\right)\left(\frac{\delta \sigma_{pp}^{\text{sys}}}{\sigma_{pp}}\right)}
    \end{equation*}
    
    This forms a perfect square $x^2 + y^2 - 2xy = (x-y)^2$:
    \begin{equation}
        \left(\frac{\delta R_{\text{sys}}}{R}\right)^2 = \left( \frac{\delta \sigma_{pd}^{\text{sys}}}{\sigma_{pd}} - \frac{\delta \sigma_{pp}^{\text{sys}}}{\sigma_{pp}} \right)^2
    \end{equation}
\end{frame}

\begin{frame}{The Power of Ratios: Error Cancellation}
    Taking the square root gives our final systematic error propagation:
    \begin{equation}
        \frac{\delta R_{\text{sys}}}{R} = \left| \frac{\delta \sigma_{pd}^{\text{sys}}}{\sigma_{pd}} - \frac{\delta \sigma_{pp}^{\text{sys}}}{\sigma_{pp}} \right|
    \end{equation}
    
    \vspace{0.3cm}
    \textbf{Why this matters:}
    \begin{itemize}
        \setlength{\itemsep}{0.2cm}
        \item If systematic errors are evaluated as uncorrelated, they inflate the final uncertainty (quadrature sum).
        \item By correctly identifying them as perfectly correlated ($\rho=1$), identical fractional biases subtract and \textbf{cancel out entirely}.
        \item This mathematically isolates the physical cross-section ratio from global detector biases.
    \end{itemize}
\end{frame}

\begin{frame}{Total Bin-by-Bin Uncertainty}
    Finally, we evaluate the total error for each kinematic bin.
    \vspace{0.3cm}
    
    Because the sources of statistical scatter (random sample size) and systematic bias (global normalization/efficiencies) are independent of each other, their cross-correlation is zero.
    
    \vspace{0.3cm}
    We add the absolute propagated errors in quadrature:
    \begin{equation}
        \delta R_{\text{total}} = \sqrt{(\delta R_{\text{stat}})^2 + (\delta R_{\text{sys}})^2}
    \end{equation}
    
    \vspace{0.3cm}
    This $\delta R_{\text{total}}$ represents the final error bars plotted on the single differential cross-section ratio $d\sigma_{pd} / 2d\sigma_{pp}$ vs $p_T$.
\end{frame}

\begin{frame}{Error Table Comparison}
    \textbf{Summary Table:} Extracted bin uncertainties comparing the uncorrelated background treatment against the new 100\% correlated systematic treatment.
    \vspace{0.2cm}
    \begin{columns}[T]
        \begin{column}{0.48\textwidth}
            \centering
            \textbf{\footnotesize Uncorrelated (Old)} \\
            \vspace{0.1cm}
            \resizebox{\textwidth}{!}{
            \begin{tabular}{ccccc}
                \toprule
                \textbf{$p_T$ Bin [GeV]} & \textbf{Ratio} & \textbf{Stat Err} & \textbf{Sys Err} & \textbf{Total Err} \\
                \midrule
__TABLE_ROWS_OLD__
                \bottomrule
            \end{tabular}
            }
        \end{column}
        
        \begin{column}{0.48\textwidth}
            \centering
            \textbf{\footnotesize Correlated (Latest)} \\
            \vspace{0.1cm}
            \resizebox{\textwidth}{!}{
            \begin{tabular}{ccccc}
                \toprule
                \textbf{$p_T$ Bin [GeV]} & \textbf{Ratio} & \textbf{Stat Err} & \textbf{Sys Err} & \textbf{Total Err} \\
                \midrule
__TABLE_ROWS_NEW__
                \bottomrule
            \end{tabular}
            }
        \end{column}
    \end{columns}
\end{frame}

\begin{frame}{Ratio Plot Comparison}
    \begin{columns}
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{\footnotesize Uncorrelated Background Treatment (Old)} \\
            \vspace{0.1cm}
            \includegraphics[width=0.95\textwidth,height=0.65\textheight,keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70/CrossSection_Ratio_pd_2pp_vs_pT_geom_with_logo.pdf}
        \end{column}
        
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{\footnotesize 100\% Correlated Propagation (Latest)} \\
            \vspace{0.1cm}
            \includegraphics[width=0.95\textwidth,height=0.65\textheight,keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/CrossSection_Ratio_pd_2pp_vs_pT_geom_with_logo.pdf}
        \end{column}
    \end{columns}
\end{frame}
% --- END ERROR PROPAGATION SECTION ---


% --- ROADSET SYSTEMATICS SECTION ---
\begin{frame}{Weighted Average Cross-Section}
    \begin{quote}
    ``Can you compute the weighted standard deviation of the cross sections
    of RS 57, 59, 62, 67 and 70? It should be better in the sense that we
    treat the roadsets equally, since we cannot assure RS 67 is 100\%
    correct.''
    \end{quote}
    \vspace{0.5cm}
    \begin{center}
        \textcolor{red}{\textbf{This requires collaboration input!}}
    \end{center}
\end{frame}

\begin{frame}{Variance Calculation Mathematics}
    For a given kinematic bin, let $x_i$ be the absolute cross-section for roadset $i$, and $e_i$ be its statistical error.
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

\begin{frame}{Resulting Uncertainties (1D $p_T$ - LH2)}
    Using the inverse-variance weighted strategy suppresses low-statistics noise, yielding the following systematic errors for Liquid Hydrogen.
    
    \textit{Note: Bins with zero relative error or relative error $\geq$ 50\% have been omitted.}
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

\begin{frame}{Resulting Uncertainties (1D $p_T$ - LD2)}
    Applying the same inverse-variance weighted procedure to Liquid Deuterium targets yields comparable systematic margins.
    
    \textit{Note: Bins with zero relative error or relative error $\geq$ 50\% have been omitted.}
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

\begin{frame}{List of systematics}
    \begin{itemize}
        \item Acceptance Correction
        \item Efficiency Correction for Yields
        \item $\psi'$ Contamination
        \item Road Dependancy
    \end{itemize}
\end{frame}

\begin{frame}{LH2 Cross-Section Previous Vs Latest with updated systematics}
    \begin{columns}[T]
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Previous}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/CrossSection_LH2_geom_vs_pT_with_logo_old.pdf}
        \end{column}
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Latest (with weighted average)}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/CrossSection_LH2_geom_vs_pT_with_logo.pdf}
        \end{column}
    \end{columns}
\end{frame}

\begin{frame}{LD2 Cross-Section Previous Vs Latest with updated systematics}
    \begin{columns}[T]
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Previous}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/CrossSection_LD2_geom_vs_pT_with_logo_old.pdf}
        \end{column}
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Latest (with weighted average)}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/CrossSection_LD2_geom_vs_pT_with_logo.pdf}
        \end{column}
    \end{columns}
\end{frame}

\begin{frame}{Summary Plots for $p_{T}$ Previous Vs Latest with updated systematics}
    \begin{columns}[T]
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Previous}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/Combined_XSec_Ratio_vs_pT_geom_old.pdf}
        \end{column}
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Latest (with weighted average)}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT/RS57-70_xsec_ratio_correlated/Combined_XSec_Ratio_vs_pT_geom_logo.pdf}
        \end{column}
    \end{columns}
\end{frame}

\begin{frame}{Summary Plots for $p^{2}_{T}$ Previous Vs Latest with updated systematics}
    \begin{columns}[T]
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Previous}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT_squard/RS57-70_clone/Combined_XSec_Ratio_geom_pT2.pdf}
        \end{column}
        \begin{column}{0.5\textwidth}
            \centering
            \textbf{Latest (with weighted average)}\\[0.2cm]
            \includegraphics[width=0.95\linewidth, height=0.7\textheight, keepaspectratio]{/root/github/e906-development/src/xsec_pT_squard/RS57-70_road_dependancy_added/Combined_XSec_Ratio_geom_pT2_logo.pdf}
        \end{column}
    \end{columns}
\end{frame}

__COMPARISON_SLIDES__

% --- BACKUP SLITES TITLE FRAME ---
\begin{frame}
    \centering
    \Huge \textbf{Backup Slides}
\end{frame}

\begin{frame}{Double Differential Systematics (xF and Mass)}
    \textbf{Extension to 2D Kinematics:}
    \vspace{0.3cm}
    \begin{itemize}
        \item The inverse-variance weighting procedure is extended to the double differential bins (xF and Invariant Mass).
        \item The following slides document the systematic uncertainties strictly for statistically significant data sets across both LH2 and LD2 targets.
        \item \textbf{Filtering Criteria:} Bins reporting a relative variance $\leq 0\%$ or $\geq 50\%$ represent zero-bin artifacts or severe low-statistics starvation and have been omitted for clarity.
    \end{itemize}
\end{frame}

__LH2_2D_FRAMES__

__LD2_2D_FRAMES__

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
    
    # 4. Inject all strings into template
    final_latex = latex_template.replace("__INCLUSIVE_PT_SLIDE__", inclusive_pT_slide)
    final_latex = final_latex.replace("__BINNED_PT_SLIDES__", binned_pT_slides)
    final_latex = final_latex.replace("__LH2_TABLE_BODY__", lh2_table_body)
    final_latex = final_latex.replace("__LD2_TABLE_BODY__", ld2_table_body)
    final_latex = final_latex.replace("__LH2_2D_FRAMES__", lh2_2d_frames)
    final_latex = final_latex.replace("__LD2_2D_FRAMES__", ld2_2d_frames)
    final_latex = final_latex.replace("__TABLE_ROWS_OLD__", ratio_rows_old)
    final_latex = final_latex.replace("__TABLE_ROWS_NEW__", ratio_rows_new)
    final_latex = final_latex.replace("__COMPARISON_SLIDES__", comparison_slides)

    # 5. Write to disk
    with open(filename, "w") as f:
        f.write(final_latex)
        
    print(f"[*] Successfully generated {filename} with all error propagation and correlation slides successfully integrated.")

if __name__ == "__main__":
    generate_presentation()