import os

# ==========================================
# Configuration
# ==========================================
old_dir = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS67/StepByStepPlots"
new_dir = "."

targets = ["LH2", "LD2", "Flask"]
plot_types = ["total", "mix"]
output_tex_file = "efficiency_comparison.tex"

# Bin counts derived from the GenerateROOTFiles.py arrays
num_mass_bins = 11  # 12 edges = 11 bins
num_xf_bins = 16    # 17 edges = 16 bins

# ==========================================
# LaTeX Document Generation
# ==========================================
def generate_tex():
    # LaTeX Preamble
    tex_content = r"""\documentclass[12pt, a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage{geometry}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{float}
\usepackage{hyperref}

% Configure hyperref for clean, clickable links
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,      
    urlcolor=cyan,
    pdftitle={Efficiency Comparison},
}

% Set margins to fit large plots side-by-side comfortably
\geometry{top=1in, bottom=1in, left=0.5in, right=0.5in}

\title{Reconstruction Efficiency Comparison: Correlated vs. Previous Results}
\author{Chatura Kuruppu}
\date{\today}

\begin{document}
\maketitle

\tableofcontents
\clearpage

\listoffigures
\clearpage
"""

    # ---------------------------------------------------------
    # PART 1: Efficiency Comparisons
    # ---------------------------------------------------------
    tex_content += "\\section{Efficiency Comparisons}\n"
    
    # Loop through targets and types to generate side-by-side figures
    for target in targets:
        tex_content += f"\\subsection{{{target} Target}}\n"
        
        for p_type in plot_types:
            # Construct exact file paths based on your prompt
            old_file = f"{old_dir}/E_{p_type}_reco_{target}.pdf"
            new_file = f"{new_dir}/E_{p_type}_reco_{target}_correlated.pdf"

            # Create the figure block
            tex_content += r"""
\begin{figure}[H]
    \centering
    \begin{subfigure}{0.49\textwidth}
        \centering
        \includegraphics[width=\linewidth]{""" + old_file + r"""}
        \caption{Old (Uncorrelated)}
    \end{subfigure}\hfill
    \begin{subfigure}{0.49\textwidth}
        \centering
        \includegraphics[width=\linewidth]{""" + new_file + r"""}
        \caption{New (Correlated)}
    \end{subfigure}
    \caption{Comparison of """ + p_type.capitalize() + r""" Reconstruction Efficiency for """ + target + r"""}
\end{figure}
"""
        # Put each target on its own page to keep it organized
        tex_content += "\\clearpage\n"

    # ---------------------------------------------------------
    # PART 2: Covariance Matrices
    # ---------------------------------------------------------
    tex_content += "\\section{Covariance Matrices}\n"
    
    for target in targets:
        tex_content += f"\\subsection{{{target} Target}}\n"
        
        for p_type in plot_types:
            tex_content += f"\\subsubsection{{{p_type.capitalize()} Yield}}\n"
            
            covar_dir = f"CovarianceMatrices_{target}_{p_type}"
            
            figs_added = 0
            for i_m in range(num_mass_bins):
                for i_x in range(num_xf_bins):
                    file_path = f"{covar_dir}/CovarianceMatrix_{target}_{p_type}_XfBin_{i_x}_MassBin_{i_m}.pdf"
                    
                    # Only include the plot if it was actually generated (not empty)
                    if os.path.exists(file_path):
                        tex_content += r"""
\begin{figure}[H]
    \centering
    \includegraphics[width=0.65\linewidth]{""" + file_path + r"""}
    \caption{Covariance Matrix for """ + target + f" ({p_type.capitalize()}) - Mass Bin {i_m}, $x_F$ Bin {i_x}" + r"""}
\end{figure}
"""
                        figs_added += 1
                        
                        # Add a clearpage every 2 figures to keep the document clean
                        if figs_added % 2 == 0:
                            tex_content += "\\clearpage\n"
            
            # Catch any dangling figure at the end of a section
            if figs_added % 2 != 0:
                tex_content += "\\clearpage\n"
            elif figs_added == 0:
                tex_content += "No non-empty covariance matrices were generated for this selection.\n\n"

    # End document
    tex_content += "\\end{document}\n"

    # Write to file
    with open(output_tex_file, "w") as f:
        f.write(tex_content)

    print(f"LaTeX script '{output_tex_file}' generated successfully.")

if __name__ == "__main__":
    generate_tex()