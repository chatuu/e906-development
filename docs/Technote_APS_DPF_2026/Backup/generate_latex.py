import os

def generate_latex():
    # Configuration
    base_dir = "/root/github/e906-development/src/xsec_pT/RS57-70"
    title = r"Measurement of Absolute Double Differential Cross-Section in Invariant Mass and $p_T$ Bins"
    author = "Chatura Kuruppu"
    
    # Define exact plot order incorporating all 38 requested files logically grouped
    plot_names = [
        # Raw & Mixed Yields
        "Y_total_LH2",
        "Y_total_LD2",
        "Y_total_Flask",
        "Y_mix_LH2",
        "Y_mix_LD2",
        "Y_mix_Flask",
        
        # Reconstruction Efficiencies
        "E_total_reco_LH2",
        "E_total_reco_LD2",
        "E_total_reco_Flask",
        "E_mix_reco_LH2",
        "E_mix_reco_LD2",
        "E_mix_reco_Flask",
        
        # Hodoscope Efficiencies
        "E_total_hodo_LH2",
        "E_total_hodo_LD2",
        "E_total_hodo_Flask",
        "E_mix_hodo_LH2",
        "E_mix_hodo_LD2",
        "E_mix_hodo_Flask",
        
        # Final Efficiencies
        "E_total_final_LH2",
        "E_total_final_LD2",
        "E_total_final_Flask",
        "E_mix_final_LH2",
        "E_mix_final_LD2",
        "E_mix_final_Flask",
        
        # Signal Efficiencies
        "E_final_signal_LH2",
        "E_final_signal_LD2",
        "E_final_signal_Flask",
        
        # Corrected & Subtracted Yields
        "Y_corrected_LH2",
        "Y_corrected_LD2",
        "Y_corrected_Flask",
        "Y_corrected_Subtracted_LH2",
        "Y_corrected_Subtracted_LD2",
        
        # Cross-Sections (Geometric vs True pT)
        "CrossSection_LH2_geom_vs_pT_with_logo",
        "CrossSection_LH2_true_pt_vs_pT_with_logo",
        "CrossSection_LD2_geom_vs_pT_with_logo",
        "CrossSection_LD2_true_pt_vs_pT_with_logo",
        
        # Ratios (Geometric vs True pT)
        "CrossSection_Ratio_pd_2pp_vs_pT_geom_with_logo",
        "CrossSection_Ratio_pd_2pp_vs_pT_true_pt_with_logo"
    ]

    # Append .pdf to each plot name
    ordered_files = [f"{name}.pdf" for name in plot_names]

    # Function to escape underscores for LaTeX text (captions/titles)
    def tex_escape(text):
        return text.replace('_', r'\_')

    # Content Strings
    inputs_latex = r"""
\textbf{Note:} We do not use RS5a and RS5b files for this study.

\vspace{0.2cm}
\textbf{ROOT files used:}
\begin{multicols}{2}
\begin{itemize}
    \tiny
    \item \texttt{merged\_RS57\_Empty\_2\_1138\_FinalData.root}
    \item \texttt{merged\_RS57\_LD2\_3\_1138\_FinalData.root}
    \item \texttt{merged\_RS57\_LH2\_1\_1138\_FinalData.root}
    \item \texttt{merged\_RS59\_Empty\_2\_466\_FinalData.root}
    \item \texttt{merged\_RS59\_LD2\_3\_466\_FinalData.root}
    \item \texttt{merged\_RS59\_LH2\_1\_465\_FinalData.root}
    \item \texttt{merged\_RS5a\_Empty\_2\_1680\_FinalData.root}
    \item \texttt{merged\_RS5a\_LD2\_3\_1689\_FinalData.root}
    \item \texttt{merged\_RS5a\_LH2\_1\_1476\_FinalData.root}
    \item \texttt{merged\_RS5b\_Empty\_2\_917\_FinalData.root}
    \item \texttt{merged\_RS5b\_LD2\_3\_918\_FinalData.root}
    \item \texttt{merged\_RS5b\_LH2\_1\_770\_FinalData.root}
    \item \texttt{merged\_RS62\_Empty\_2\_1234\_FinalData.root}
    \item \texttt{merged\_RS62\_LD2\_3\_1237\_FinalData.root}
    \item \texttt{merged\_RS62\_LH2\_1\_1234\_FinalData.root}
    \item \texttt{merged\_RS67\_Empty\_2\_14\_FinalData.root}
    \item \texttt{merged\_RS67\_LD2\_3\_15\_FinalData.root}
    \item \texttt{merged\_RS67\_LH2\_1\_5\_FinalData.root}
    \item \texttt{merged\_RS70\_Empty\_2\_267\_FinalData.root}
    \item \texttt{merged\_RS70\_LD2\_3\_266\_FinalData.root}
    \item \texttt{merged\_RS70\_LH2\_1\_264\_FinalData.root}
\end{itemize}
\end{multicols}
"""

    kinematic_bins_latex = r"""
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
        f.write(r"""\documentclass[aspectratio=169]{beamer}
\usetheme{Madrid}
\usecolortheme{default}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{multicol}
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

\begin{frame}{Input Files Used}
    %s
\end{frame}

\begin{frame}{Kinematic Phase Space}
    %s
\end{frame}

\begin{frame}{Event Selection Criteria}
    \footnotesize
    %s
\end{frame}
""" % (base_dir, title, author, inputs_latex, kinematic_bins_latex, event_selection_latex))

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
\usepackage{multicol}
\graphicspath{{%s/}}

\title{%s}
\author{%s}
\date{\today}

\begin{document}

\maketitle

\section{Input Files Used}
%s

\section{Kinematic Phase Space}
%s

\section{Event Selection Criteria}
%s

\clearpage
\section{Plots and Distributions}
""" % (base_dir, title, author, inputs_latex, kinematic_bins_latex, event_selection_latex))

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