import os

def generate_beamer_tex(output_filename="comparison_slides.tex"):
    # Variables and their exact string matches to the PDF filenames generated earlier
    variables = [
        {"name": "pT", "title": "p_{T}"},
        {"name": "xF", "title": "x_{F}"},
        {"name": "Mass", "title": "Mass"}
    ]

    # LaTeX document header (Not an f-string, so single braces are fine here)
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

    # Generate a slide for each variable
    for var in variables:
        v_name = var["name"]
        v_title = var["title"]
        
        # Notice we use the exact naming convention from your previous Python script
        lh2_file = f"{v_name}_Distribution_LH2.pdf"
        ld2_file = f"{v_name}_Distribution_LD2.pdf"
        
        # This IS an f-string, so all literal LaTeX braces must be doubled (e.g., {{frame}})
        slide = f"""
\\begin{{frame}}{{{v_title} Distribution: LH2 vs LD2}}
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
        tex_content += slide

    # Close document
    tex_content += r"\end{document}" + "\n"

    # Write to file
    with open(output_filename, "w") as f:
        f.write(tex_content)
        
    print(f"Successfully generated LaTeX file: {output_filename}")

if __name__ == "__main__":
    generate_beamer_tex()