import os

def generate_tex():
    # Enforced ordering based on your requirements
    roadsets = ["RS57", "RS59", "RS62", "RS67", "RS70", "RS57-70"]
    targets = ["LH2", "LD2", "Flask"]
    plot_types = ["total", "mix", "corrected", "corrected_Subtracted"]

    # LaTeX document header with Beamer templates
    tex_code = r"""\documentclass[aspectratio=169]{beamer}
\usetheme{Madrid}        % Professional Beamer template
\usecolortheme{whale}    % Clean blue color palette
\usepackage{graphicx}

% Customize the footline to be a bit cleaner if desired
\setbeamertemplate{navigation symbols}{}

\title{Yield Distributions for $p_T$ Cross-Section Study}
\author{Chatura Kuruppu}
\date{\today}

\begin{document}

\begin{frame}
    \titlepage
\end{frame}
"""

    # Generate a slide for each plot
    for rs in roadsets:
        for tgt in targets:
            for ptype in plot_types:
                # The Flask target does not have a subtracted plot
                if ptype == "corrected_Subtracted" and tgt == "Flask":
                    continue
                
                filename = f"/root/github/e906-development/src/xsec_pT/{rs}/Y_{ptype}_{tgt}.pdf"
                
                # Check if file exists to prevent LaTeX compilation errors just in case
                if not os.path.exists(filename):
                    print(f"Warning: {filename} not found. Skipping slide.")
                    continue

                # Escape underscores for the LaTeX frame title
                display_type = ptype.replace('_', r'\_')

                tex_code += f"\n\\begin{{frame}}{{{rs} | Target: {tgt} | Type: {display_type}}}\n"
                tex_code += r"    \centering" + "\n"
                tex_code += f"    \\includegraphics[width=\\textwidth, height=0.85\\textheight, keepaspectratio]{{{filename}}}\n"
                tex_code += r"\end{frame}" + "\n"

    # Close document
    tex_code += r"\end{document}" + "\n"

    # Write output
    output_file = "yield_slides.tex"
    with open(output_file, "w") as f:
        f.write(tex_code)
    
    print(f"Successfully generated {output_file}")

if __name__ == "__main__":
    generate_tex()