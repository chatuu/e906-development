import os

# ==========================================
# Configuration
# ==========================================
old_dir = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS67/StepByStepPlots"
new_dir = "."

targets = ["LH2", "LD2", "Flask"]
plot_types = ["total", "mix"]
output_tex_file = "efficiency_comparison.tex"

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

% Set margins to fit large plots side-by-side comfortably
\geometry{top=1in, bottom=1in, left=0.5in, right=0.5in}

\title{Reconstruction Efficiency Comparison: Correlated vs. Previous Results}
\author{Chatura Kuruppu}
\date{\today}

\begin{document}
\maketitle
\tableofcontents
\clearpage
"""

    # Loop through targets and types to generate side-by-side figures
    for target in targets:
        tex_content += f"\\section{{{target} Target}}\n"
        
        for p_type in plot_types:
            # Construct exact file paths based on your prompt
            old_file = f"{old_dir}/E_{p_type}_reco_{target}.pdf"
            new_file = f"{new_dir}/E_{p_type}_reco_{target}_correlated.pdf"

            # Create the figure block
            tex_content += r"""
\begin{figure}[htbp]
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

    # End document
    tex_content += "\\end{document}\n"

    # Write to file
    with open(output_tex_file, "w") as f:
        f.write(tex_content)

    print(f"LaTeX script '{output_tex_file}' generated successfully.")

if __name__ == "__main__":
    generate_tex()
