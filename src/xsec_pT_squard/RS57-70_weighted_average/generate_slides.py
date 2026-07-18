"""
generate_slides.py
Generates a LaTeX Beamer presentation for Drell-Yan Cross-Section plots.
"""

import os

# Define the xF bins used in your analysis config
XF_BINS = [0.1, 0.16, 0.22, 0.28, 0.34, 0.40, 0.46, 0.58, 0.70, 0.95]
n_bins = len(XF_BINS) - 1

def generate_tex():
    # Standard Beamer Header
    tex_code = r"""\documentclass[aspectratio=169]{beamer}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{caption}

\usetheme{Madrid}
\usecolortheme{default}

% Title Page Information
\title[DY Cross-Sections]{Drell-Yan Absolute Cross-Sections and Ratios}
\subtitle{Binned in $x_F$ and $p_T$}
\author{Chatura Kuruppu}
\institute[Fermilab / NMSU]{New Mexico State University \\ Fermi National Accelerator Laboratory}
\date{\today}

\begin{document}

\begin{frame}
  \titlepage
\end{frame}

\begin{frame}{Outline}
  \tableofcontents
\end{frame}
"""

    # Generate slides for each xF bin
    for i in range(n_bins):
        xf_low = XF_BINS[i]
        xf_high = XF_BINS[i+1]
        
        section_title = f"$x_F \\in [{xf_low:.2f}, {xf_high:.2f})$"
        tex_code += f"\n\\section{{{section_title}}}\n"
        
        # Slide 1: Individual LH2 and LD2 Cross-Sections
        tex_code += f"""
\\begin{{frame}}{{LH2 and LD2 Single Differential Cross-Sections}}
  \\framesubtitle{{{section_title}}}
  \\begin{{columns}}
    \\begin{{column}}{{0.5\\textwidth}}
      \\centering
      \\includegraphics[width=\\textwidth]{{CrossSection_LH2_vs_pT_geom_xF{i}_with_logo.pdf}}
    \\end{{column}}
    \\begin{{column}}{{0.5\\textwidth}}
      \\centering
      \\includegraphics[width=\\textwidth]{{CrossSection_LD2_vs_pT_geom_xF{i}_with_logo.pdf}}
    \\end{{column}}
  \\end{{columns}}
\\end{{frame}}
"""
        
        # Slide 2: Ratio Plot
        tex_code += f"""
\\begin{{frame}}{{Cross-Section Ratio $\\sigma_{{pd}} / 2\\sigma_{{pp}}$}}
  \\framesubtitle{{{section_title}}}
  \\centering
  \\includegraphics[height=0.8\\textheight]{{CrossSection_Ratio_pd_2pp_vs_pT_geom_xF{i}_with_logo.pdf}}
\\end{{frame}}
"""

        # Slide 3: 3-in-1 Combined Canvas Plot
        tex_code += f"""
\\begin{{frame}}{{Combined Overlay and Ratio}}
  \\framesubtitle{{{section_title}}}
  \\centering
  \\includegraphics[height=0.85\\textheight]{{Combined_XSec_Ratio_vs_pT_geom_xF{i}.pdf}}
\\end{{frame}}
"""

    tex_code += r"\end{document}"

    output_filename = "dy_cross_section_slides.tex"
    with open(output_filename, "w") as f:
        f.write(tex_code)

    print(f"Successfully generated {output_filename}")

if __name__ == "__main__":
    generate_tex()