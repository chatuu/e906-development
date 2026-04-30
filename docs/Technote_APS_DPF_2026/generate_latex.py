import os

def generate_latex():
    # Configuration
    #base_dir = "/root/github/e906-development/src/xsec_pT_mass/RS57-70"
    base_dir = "/root/github/e906-development/src/xsec_pT/RS57-70"
    title = r"Measurement of Absolute Double Differential Cross-Section in Invariant Mass and $p_T$ Bins"
    author = "Chatura Kuruppu"
    
    # 1. Define the targets and plotting prefixes
    targets = ["LH2", "LD2", "Flask"]
    kinematic_prefixes = [
        "Y_total",
        "Y_mix",
        "E_total_reco",
        "E_mix_reco",
        "E_total_hodo",
        "E_mix_hodo",
        "E_total_final",
        "E_mix_final",
        "E_final_signal",
        "Y_corrected"
    ]

    # Build the ordered list of files based on user rules
    ordered_files = []

    # Rules 2-11: Kinematic and efficiency plots for all targets
    for target in targets:
        for prefix in kinematic_prefixes:
            ordered_files.append(f"{prefix}_{target}.pdf")

    # Rule 12: Subtracted Yields (Only for LH2 and LD2)
    ordered_files.append("Y_corrected_Subtracted_LH2.pdf")
    ordered_files.append("Y_corrected_Subtracted_LD2.pdf")

    # Rule 13: LH2 Cross-sections (New 1D pT plot)
    ordered_files.append("CrossSection_LH2_vs_pT_with_logo.pdf")

    # Rule 14: LD2 Cross-sections (New 1D pT plot)
    ordered_files.append("CrossSection_LD2_vs_pT_with_logo.pdf")

    # Function to escape underscores for LaTeX text (captions/titles)
    def tex_escape(text):
        return text.replace('_', r'\_')

    # Content Strings
    inputs_and_bins_latex = r"""
\textbf{Input Data Files:}
\begin{itemize}
    \item Merged ROOT files containing reconstruction and hodoscope efficiencies.
    \item \textbf{Roadsets:} RS57, RS59, RS62, RS67 (run 3089), and RS70.
    \item \textbf{Targets:} Liquid Hydrogen (LH2), Liquid Deuterium (LD2), and Empty Flask.
\end{itemize}
\vspace{0.3cm}
\textbf{Kinematic Binning Definition:}
\begin{itemize}
    \item \textbf{Mass Bins (GeV):} [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7]
    \item \textbf{$p_T$ Bins (GeV):} [0.0, 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8]
    \item \textbf{$x_F$ Bins:} [0.0, 0.85] with a fixed step width of 0.05.
\end{itemize}
"""

    event_selection_latex = r"""
\textbf{Dimuon Kinematics:}
\begin{itemize}
    \item $4.2 < \text{Mass} < 8.8$ GeV, $-0.1 < x_F < 0.95$, $0.05 < x_T \le 0.58$, and $|\cos\theta| < 0.5$.
    \item Vertex cuts: $-280 < dz < -5$, target/dump tracking separation criteria met.
    \item Transverse momentum bounds: $|dp_x| < 1.8$, $|dp_y| < 2.0$, and $dp_x^2 + dp_y^2 < 5$.
\end{itemize}
\vspace{0.2cm}
\textbf{Track Quality \& Acceptance:}
\begin{itemize}
    \item $\chi^2_{\text{target}} < 15$, $\chi^2_{\text{dimuon}} < 18$, and $\chi^2/(\text{nHits}-5) < 12$.
    \item Hit requirements: nHits $> 13$ per track, $\text{nHits}_1+\text{nHits}_2 > 29$, St1 total $> 8$.
    \item Track vertex longitudinal limits: $-320 < z_v < -5$.
    \item A dynamic beam offset correction is applied to y-coordinates: offset is $1.6$ for runID $\ge 11000$ and $0.4$ otherwise.
\end{itemize}
\vspace{0.2cm}
\textbf{Detector Occupancy Limits:}
\begin{itemize}
    \item Drift chamber hits: $D1 < 400$, $D2 < 400$, $D3 < 400$, and total $D1+D2+D3 < 1000$.
\end{itemize}
"""

    # ==========================================
    # Generate Beamer Slides
    # ==========================================
    beamer_file = "slides.tex"
    with open(beamer_file, "w") as f:
        # Added aspectratio=169 for widescreen, and the Madrid theme for professional formatting
        f.write(r"""\documentclass[aspectratio=169]{beamer}
\usetheme{Madrid}
\usecolortheme{default}
\usepackage{graphicx}
\usepackage{hyperref}
\setbeamertemplate{navigation symbols}{}
\graphicspath{{%s/}}

\title{%s}
\author{%s}
\institute{SeaQuest Experiment (E906) \\ Fermi National Accelerator Laboratory}
\date{\today}

\begin{document}

\begin{frame}
    \titlepage
\end{frame}

\begin{frame}{Input Files and Bin Ranges}
    %s
\end{frame}

\begin{frame}{Event Selection Criteria}
    \footnotesize
    %s
\end{frame}
""" % (base_dir, title, author, inputs_and_bins_latex, event_selection_latex))

        for file in ordered_files:
            # Skip missing files to prevent compilation crashes
            if not os.path.exists(os.path.join(base_dir, file)):
                print(f"Warning: {file} not found in {base_dir}. Skipping in presentation.")
                continue

            frame_title = tex_escape(file.replace('.pdf', ''))
            f.write(r"""
\begin{frame}{%s}
    \begin{center}
        \includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{%s}
    \end{center}
\end{frame}
""" % (frame_title, file))

        f.write("\n\\end{document}\n")

    # ==========================================
    # Generate Article Document
    # ==========================================
    article_file = "document.tex"
    with open(article_file, "w") as f:
        f.write(r"""\documentclass[12pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{float}
\usepackage{hyperref}
\graphicspath{{%s/}}

\title{%s}
\author{%s}
\date{\today}

\begin{document}

\maketitle

\section{Input Files and Bin Ranges}
%s

\section{Event Selection Criteria}
%s

\clearpage
\section{Plots and Distributions}
""" % (base_dir, title, author, inputs_and_bins_latex, event_selection_latex))

        for idx, file in enumerate(ordered_files):
            # Skip missing files to prevent compilation crashes
            if not os.path.exists(os.path.join(base_dir, file)):
                continue

            caption = tex_escape(file.replace('.pdf', ''))
            f.write(r"""
\begin{figure}[H]
    \centering
    \includegraphics[width=0.9\textwidth]{%s}
    \caption{%s}
\end{figure}
""" % (file, caption))
            
            # Add a clearpage every 2 images to keep the formatting neat
            if (idx + 1) % 2 == 0:
                f.write("\\clearpage\n")

        f.write("\n\\end{document}\n")

    print(f"Successfully generated '{beamer_file}' and '{article_file}'.")

if __name__ == "__main__":
    generate_latex()