import os

def make_frame_2cols(title, img_left, img_right):
    """Helper to generate a two-column Beamer slide."""
    return f"""
\\begin{{frame}}{{{title}}}
    \\begin{{columns}}
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\includegraphics[width=\\textwidth, keepaspectratio]{{{img_left}}}
        \\end{{column}}
        \\begin{{column}}{{0.5\\textwidth}}
            \\centering
            \\includegraphics[width=\\textwidth, keepaspectratio]{{{img_right}}}
        \\end{{column}}
    \\end{{columns}}
\\end{{frame}}
"""

def make_frame_1col(title, img):
    """Helper to generate a single-image Beamer slide."""
    return f"""
\\begin{{frame}}{{{title}}}
    \\begin{{center}}
        \\includegraphics[height=0.85\\textheight, keepaspectratio]{{{img}}}
    \\end{{center}}
\\end{{frame}}
"""

def main():
    output_tex = "acceptance_slides.tex"

    # OLD_DIR remains constant across all roadsets
    OLD_DIR = os.path.expanduser("~/github/e906-development/src/AcceptanceCorrection/acceptance_no_unfolding_bins/").replace('\\', '/')

    # LaTeX document header using a 16:9 widescreen aspect ratio
    latex_content = r"""\documentclass[aspectratio=169]{beamer}
\usepackage{graphicx}
\usepackage{amsmath}

% Use a blue theme that natively includes slide numbers at the bottom right
\usetheme{Madrid}
\usecolortheme{whale} % Enforces a strong blue palette

% Remove navigation symbols for a cleaner look
\setbeamertemplate{navigation symbols}{}

\title{Detector Acceptance Study \\ (Existing Vs Latest)}
\author{Chatura Kuruppu}
\institute{New Mexico State University}
\date{\today}

\begin{document}

\begin{frame}
    \titlepage
\end{frame}

% ==========================================
% INTRO SLIDES
% ==========================================
\begin{frame}{Overview}
    \begin{itemize}
        \item Re-Calculate Acceptance Correction for each roadset (long term)
        \item Currently calculated acceptance corrections by using Kenichi's new root files
        \item Compare Latest Acceptance Correction Vs Existing Acceptance Correction
    \end{itemize}
\end{frame}

\begin{frame}{Files Used \& Kinematic Bins}
    \small
    \textbf{Latest Files Used (Kenichi):}
    \begin{itemize}
        \item \textbf{Thrown:} \texttt{rs<57/59/62/67/70>\_<lh2/ld2>\_acc.root:tree\_4pi:event}
        \item \textbf{Accepted:} \texttt{rs<57/59/62/67/70>\_<lh2/ld2>\_acc.root:tree\_acc:event}
    \end{itemize}
    \vspace{0.1cm}
    \textbf{Existing Files Used (Hugo):}
    \begin{itemize}
        \item \textbf{Thrown:} \texttt{mc\_drellyan\_<LH2/LD2>\_M027\_S001\_4pi\_pTxFweight\_v2.root}
        \item \textbf{Accepted:} \texttt{mc\_drellyan\_<LH2/LD2>\_M027\_S001\_clean\_occ\_pTxFweight\_v2.root}
    \end{itemize}
    \vspace{0.1cm}
    \textbf{Kinematic Bin Boundaries:}
    \begin{itemize}
        \item \textbf{Mass [GeV]:} [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.8]
        \item \textbf{$x_F$:} [0.0, 0.05, ..., 0.80]
        \item \textbf{$p_T$ [GeV/c]:} [0.0, 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8]
    \end{itemize}
\end{frame}

\begin{frame}{Event Selection Criteria (Chuck's Cuts)}
    \begin{columns}[T]
        \begin{column}{0.5\textwidth}
            \textbf{Dimuon Kinematics:}
            \begin{itemize}
                \item $4.2 < \text{Mass} < 8.8 \text{ GeV}$
                \item $-0.1 < x_F < 0.95$
                \item $0.05 < x_T \leq 0.58$
                \item $38 < dp_z < 116 \text{ GeV}$
                \item $dp_x^2 + dp_y^2 < 5.0 \text{ GeV}^2$
            \end{itemize}
            \vspace{0.2cm}
            \textbf{Target \& Vertex:}
            \begin{itemize}
                \item $-280 < dz < -5 \text{ cm}$
                \item $|dx| < 0.25 \text{ cm}$
                \item $|dy - 1.6| < 0.22 \text{ cm}$
            \end{itemize}
        \end{column}
        
        \begin{column}{0.5\textwidth}
            \textbf{Track Quality:}
            \begin{itemize}
                \item $\chi^2_{\text{dimuon}} < 18$
                \item Track separation $< 270 \text{ cm}$
                \item $p_{z,\text{st1}} \in [9, 75] \text{ GeV}$
                \item $\text{nHits}_1 > 13, \text{nHits}_2 > 13$
            \end{itemize}
            \vspace{0.2cm}
            \textbf{Generator-Level Fiducial Cuts (Thrown):}
            \begin{itemize}
                \item $4.2 < \text{Mass} < 8.8 \text{ GeV}$
                \item $-0.1 < x_F < 0.95$
                \item $0.0 < p_T \leq 3.0 \text{ GeV/c}$
            \end{itemize}
        \end{column}
    \end{columns}
\end{frame}
"""

    # List of all roadsets to loop through
    roadsets = [57, 59, 62, 67, 70]
    
    for rs in roadsets:
        rs_str = f"RS{rs}"
        
        # Ensure RS67 points to its original specific directories
        if rs == 67:
            NEW_DIR = os.path.expanduser("~/github/e906-development/src/NewAcceptanceCorrection/RS67").replace('\\', '/')
            COMP_DIR = os.path.expanduser("~/github/e906-development/src/Existing_Vs_Latest_Acceptance_Comparison").replace('\\', '/')
        else:
            # For all other roadsets, use the parameterized directories
            NEW_DIR = os.path.expanduser(f"~/github/e906-development/src/NewAcceptanceCorrection/{rs_str}").replace('\\', '/')
            COMP_DIR = os.path.expanduser(f"~/github/e906-development/src/NewAcceptanceCorrection/{rs_str}/Existing_Vs_Latest_Acceptance_Comparison").replace('\\', '/')

        # ==========================================
        # 1. Invariant Mass Studies
        # ==========================================
        latex_content += f"\\section{{{rs_str}: Invariant Mass Plots}}\n"
        latex_content += make_frame_2cols(f"{rs_str}: Existing (L) vs Latest (R): Mass (Integrated)", 
                                          f"{OLD_DIR}/acceptance_overlay_Integrated_mass.pdf", 
                                          f"{NEW_DIR}/acceptance_overlay_Integrated_mass.pdf")
        
        latex_content += make_frame_2cols(f"{rs_str}: Existing (L) vs Latest (R): Mass Acc. (All $x_F$, $p_T$)", 
                                          f"{OLD_DIR}/Acceptance_Mass_All_xF_pT.pdf", 
                                          f"{NEW_DIR}/Acceptance_Mass_All_xF_pT.pdf")
        
        latex_content += make_frame_2cols(f"{rs_str}: LH2 (L) vs LD2 (R): Latest/Existing Mass Acc.", 
                                          f"{COMP_DIR}/Compare_Latest_Vs_Existing_LH2_Acceptance_Mass_All_xF_pT.pdf", 
                                          f"{COMP_DIR}/Compare_Latest_Vs_Existing_LD2_Acceptance_Mass_All_xF_pT.pdf")

        # ==========================================
        # 2. xF Studies
        # ==========================================
        latex_content += f"\\section{{{rs_str}: $x_F$ Plots}}\n"
        latex_content += make_frame_2cols(f"{rs_str}: Existing (L) vs Latest (R): $x_F$ (Integrated)", 
                                          f"{OLD_DIR}/acceptance_overlay_Integrated_xF.pdf", 
                                          f"{NEW_DIR}/acceptance_overlay_Integrated_xF.pdf")
        
        latex_content += make_frame_2cols(f"{rs_str}: Existing (L) vs Latest (R): $x_F$ Acc. (All Mass, $p_T$)", 
                                          f"{OLD_DIR}/Acceptance_xF_All_Mass_pT.pdf", 
                                          f"{NEW_DIR}/Acceptance_xF_All_Mass_pT.pdf")

        latex_content += make_frame_2cols(f"{rs_str}: LH2 (L) vs LD2 (R): Latest/Existing $x_F$ Acc.", 
                                          f"{COMP_DIR}/Compare_Latest_Vs_Existing_LH2_Acceptance_xF_All_Mass_pT.pdf", 
                                          f"{COMP_DIR}/Compare_Latest_Vs_Existing_LD2_Acceptance_xF_All_Mass_pT.pdf")

        # ==========================================
        # 3. pT Studies
        # ==========================================
        latex_content += f"\\section{{{rs_str}: $p_T$ Plots}}\n"
        latex_content += make_frame_2cols(f"{rs_str}: Existing (L) vs Latest (R): $p_T$ (Integrated)", 
                                          f"{OLD_DIR}/acceptance_overlay_Integrated_pT.pdf", 
                                          f"{NEW_DIR}/acceptance_overlay_Integrated_pT.pdf")
        
        latex_content += make_frame_2cols(f"{rs_str}: Existing (L) vs Latest (R): $p_T$ Acc. (All Mass, $x_F$)", 
                                          f"{OLD_DIR}/Acceptance_pT_All_Mass_xF.pdf", 
                                          f"{NEW_DIR}/Acceptance_pT_All_Mass_xF.pdf")
                                          
        latex_content += make_frame_2cols(f"{rs_str}: LH2 (L) vs LD2 (R): Latest/Existing $p_T$ Acc.", 
                                          f"{COMP_DIR}/Compare_Latest_Vs_Existing_LH2_Acceptance_pT_All_Mass_xF.pdf", 
                                          f"{COMP_DIR}/Compare_Latest_Vs_Existing_LD2_Acceptance_pT_All_Mass_xF.pdf")

        # ==========================================
        # 4. Mass Sliced by xF Bin Overlays & Acceptances
        # ==========================================
        latex_content += f"\\section{{{rs_str}: Mass Sliced by $x_F$ Overlays}}\n"
        
        # 17 edges = 16 bins for xF
        xf_edges = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
                    0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80]
                    
        for i in range(16):
            xf_min = xf_edges[i]
            xf_max = xf_edges[i+1]
            bin_label = f"${xf_min:.2f} \\leq x_F < {xf_max:.2f}$"
            
            # 1. Overlay Slide
            title_overlay = f"{rs_str}: Existing vs Latest: Mass Overlay ({bin_label})"
            filename_overlay = f"acceptance_overlay_mass_sliced_by_xF_bin{i}.pdf"
            latex_content += make_frame_2cols(title_overlay, 
                                              f"{OLD_DIR}/{filename_overlay}", 
                                              f"{NEW_DIR}/{filename_overlay}")
            
            # 2. Acceptance Split-Canvas Slide
            title_acc = f"{rs_str}: Existing vs Latest: Mass Acc. ({bin_label})"
            filename_acc = f"Acceptance_Mass_All_pT_xF_bin{i}.pdf"
            latex_content += make_frame_2cols(title_acc, 
                                              f"{OLD_DIR}/{filename_acc}", 
                                              f"{NEW_DIR}/{filename_acc}")

            # 3. Latest/Existing Comparison Split-Canvas Slide
            title_comp = f"{rs_str}: LH2 (L) vs LD2 (R): Latest/Existing Mass Acc. ({bin_label})"
            filename_comp_lh2 = f"Compare_Latest_Vs_Existing_LH2_Acceptance_Mass_All_pT_xF_bin{i}.pdf"
            filename_comp_ld2 = f"Compare_Latest_Vs_Existing_LD2_Acceptance_Mass_All_pT_xF_bin{i}.pdf"
            latex_content += make_frame_2cols(title_comp, 
                                              f"{COMP_DIR}/{filename_comp_lh2}", 
                                              f"{COMP_DIR}/{filename_comp_ld2}")
    
    # ==========================================
    # 5. Conclusions
    # ==========================================
    latex_content += r"\section{Conclusions}" + "\n"
    latex_content += r"""\begin{frame}{Conclusions}
    \begin{itemize}
        \item Acceptance corrections calculated for the double differential cross-section study was minimally impacted.
        \vspace{0.3cm}
        \item Acceptance Vs $x_F$ and Acceptance Vs $p_T$ plots demonstrate larger differences.
    \end{itemize}
\end{frame}
"""

    # Close the document
    latex_content += r"\end{document}"

    # Write out to file
    with open(output_tex, "w") as f:
        f.write(latex_content)

    print(f"Successfully generated {output_tex} for all roadsets.")

if __name__ == "__main__":
    main()