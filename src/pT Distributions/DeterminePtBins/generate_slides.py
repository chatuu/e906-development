import os

def format_bin_title(bin_str):
    """Converts a string like '0p0_to_0p5' into a LaTeX math string '$0.0 \\leq p_{T} < 0.5$'."""
    parts = bin_str.replace('p', '.').split('_to_')
    if len(parts) == 2:
        return f"${parts[0]} \\leq p_{{T}} < {parts[1]}$ GeV/c"
    return bin_str

def make_slide(title, lh2_file, ld2_file):
    """Generates a standard two-column Beamer slide."""
    return f"""
\\begin{{frame}}{{{title}}}
    \\begin{{columns}}[T] 
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\textbf{{LH2 Target}}\\\\[0.2cm]
            \\includegraphics[width=\\linewidth, keepaspectratio]{{{lh2_file}}}
        \\end{{column}}
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\textbf{{LD2 Target}}\\\\[0.2cm]
            \\includegraphics[width=\\linewidth, keepaspectratio]{{{ld2_file}}}
        \\end{{column}}
    \\end{{columns}}
\\end{{frame}}
"""

def generate_beamer_tex(output_filename="comparison_slides.tex"):
    
    # pT Bins specified in your directory
    pt_bins = [
        "pT_0p0_to_0p5",
        "pT_0p5_to_1p0",
        "pT_1p0_to_1p5",
        "pT_1p5_to_1p8"
    ]

    # Variables requiring binned plots
    binned_vars = [
        {"name": "xF", "title": "$x_{F}$"},
        {"name": "Mass", "title": "Mass"}
    ]

    # LaTeX document header
    tex_content = r"""\documentclass[aspectratio=169]{beamer}
\usepackage{graphicx}
\usepackage{caption}

% Theme settings
\usetheme{Madrid}
\usecolortheme{default}
\setbeamertemplate{navigation symbols}{}

\title{Data-Mix vs Messy MC Comparisons}
\subtitle{LH2 and LD2 Targets}
\author{Target Comparison Plots}
\date{\today}

\begin{document}

\begin{frame}
    \titlepage
\end{frame}
"""

    # --- 1. Inclusive Distribution (pT Only) ---
    tex_content += r"\section{Inclusive Distributions (Data - Mix)}" + "\n"
    title_pT = "$p_{T}$ Inclusive (Data - Mix): LH2 vs LD2"
    lh2_file_pT = "pT_Distribution_LH2.pdf"
    ld2_file_pT = "pT_Distribution_LD2.pdf"
    tex_content += make_slide(title_pT, lh2_file_pT, ld2_file_pT)

    # --- 2. pT-Binned Distributions (xF and Mass) ---
    tex_content += r"\section{$p_{T}$-Binned Distributions}" + "\n"
    
    for var in binned_vars:
        v_name = var["name"]
        v_title = var["title"]
        for pt_bin in pt_bins:
            lh2_file = f"{v_name}_Distribution_LH2_{pt_bin}.pdf"
            ld2_file = f"{v_name}_Distribution_LD2_{pt_bin}.pdf"
            bin_label = format_bin_title(pt_bin.replace('pT_', ''))
            title = f"{v_title} ({bin_label}): LH2 vs LD2"
            tex_content += make_slide(title, lh2_file, ld2_file)

    # Close document
    tex_content += r"\end{document}" + "\n"

    # Write to file
    with open(output_filename, "w") as f:
        f.write(tex_content)
        
    print(f"Successfully generated LaTeX file: {output_filename}")

if __name__ == "__main__":
    generate_beamer_tex()