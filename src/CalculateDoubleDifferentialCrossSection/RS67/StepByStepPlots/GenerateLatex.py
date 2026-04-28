import os

def create_latex_script(output_filename="Analysis_Report.tex"):
    # Define the file pairs according to your requested order
    # Format: (Left Image, Right Image, Caption)
    
    # --- LH2 Section ---
    lh2_plots = [
        ("Y_total_LH2.pdf", "Y_mix_LH2.pdf", "LH2: Total Yield vs Mix Yield"),
        ("E_total_reco_LH2.pdf", "E_mix_reco_LH2.pdf", "LH2: Reconstruction Efficiency (Total vs Mix)"),
        ("E_total_hodo_LH2.pdf", "E_mix_hodo_LH2.pdf", "LH2: Hodoscope Efficiency (Total vs Mix)"),
        ("E_total_final_LH2.pdf", "E_mix_final_LH2.pdf", "LH2: Final Efficiency (Total vs Mix)"),
        ("Y_corrected_LH2.pdf", "E_final_signal_LH2.pdf", "LH2: Corrected Yield and Signal Efficiency")
    ]

    # --- Flask Section ---
    flask_plots = [
        ("Y_total_Flask.pdf", "Y_mix_Flask.pdf", "Flask: Total Yield vs Mix Yield"),
        ("E_total_reco_Flask.pdf", "E_mix_reco_Flask.pdf", "Flask: Reconstruction Efficiency (Total vs Mix)"),
        ("E_total_hodo_Flask.pdf", "E_mix_hodo_Flask.pdf", "Flask: Hodoscope Efficiency (Total vs Mix)"),
        ("E_total_final_Flask.pdf", "E_mix_final_Flask.pdf", "Flask: Final Efficiency (Total vs Mix)"),
        ("Y_corrected_Flask.pdf", "E_final_signal_Flask.pdf", "Flask: Corrected Yield and Signal Efficiency")
    ]

    # --- 1. LaTeX Header (Raw String) ---
    latex_content = r"""\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage{subcaption}
\usepackage[margin=0.5in]{geometry}
\usepackage{float}

\title{Kinematic Bin Analysis Report}
\author{Automated Script}
\date{\today}

\begin{document}

\maketitle

\section*{LH2 Data Analysis}
"""

    # --- Helper Function to Create Figure Blocks ---
    # We use %s formatting here. It is much safer for LaTeX than f-strings 
    # because it doesn't conflict with LaTeX's {curly braces}.
    def get_figure_block(left_file, right_file, caption_text):
        safe_caption = caption_text.replace("_", r"\_")
        safe_left_label = left_file.replace("_", r"\_")
        safe_right_label = right_file.replace("_", r"\_")

        # The %s placeholders will be replaced by the variables at the end
        block = r"""
\begin{figure}[H]
    \centering
    \begin{subfigure}[b]{0.49\textwidth}
        \centering
        \includegraphics[width=\textwidth]{%s}
        \caption{%s}
    \end{subfigure}
    \hfill
    \begin{subfigure}[b]{0.49\textwidth}
        \centering
        \includegraphics[width=\textwidth]{%s}
        \caption{%s}
    \end{subfigure}
    \caption{%s}
\end{figure}
""" % (left_file, safe_left_label, right_file, safe_right_label, safe_caption)
        return block

    # --- 2. Add LH2 Plots ---
    for i, (left, right, caption) in enumerate(lh2_plots):
        latex_content += get_figure_block(left, right, caption)
        
        # Add page break and new section after the last LH2 plot
        if i == len(lh2_plots) - 1:
            latex_content += r"\newpage" + "\n" + r"\section*{Flask Data Analysis}" + "\n"

    # --- 3. Add Flask Plots ---
    for left, right, caption in flask_plots:
        latex_content += get_figure_block(left, right, caption)

    # --- 4. LaTeX Footer ---
    latex_content += r"\end{document}"

    # --- 5. Write to File ---
    with open(output_filename, "w") as f:
        f.write(latex_content)
    
    print(f"Successfully generated LaTeX script: {output_filename}")
    print("To compile (if you have pdflatex installed):")
    print(f"  pdflatex {output_filename}")

if __name__ == "__main__":
    create_latex_script()