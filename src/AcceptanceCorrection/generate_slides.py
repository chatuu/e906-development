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

    # Number of bins based on your numpy arrays
    num_xf_bins = 16 # np.arange(0.0, 0.85, 0.05) yields 17 edges = 16 bins
    num_pt_bins = 7  # [0., 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8] = 8 edges = 7 bins

    # LaTeX document header using a 16:9 widescreen aspect ratio
    latex_content = r"""\documentclass[aspectratio=169]{beamer}
\usepackage{graphicx}
\usepackage{amsmath}

% Use a blue theme that natively includes slide numbers at the bottom right
\usetheme{Madrid}
\usecolortheme{whale} % Enforces a strong blue palette

% Remove navigation symbols for a cleaner look
\setbeamertemplate{navigation symbols}{}

\title{Detector Acceptance Study}
\author{Chatura Kuruppu}
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
        \item Files Used \& Event Selection
        \item Invariant Mass Studies
        \item $x_F$ Studies (Yields \& Acceptances)
        \item $p_T$ Studies (Yields \& Acceptances)
        \item Track Kinematics ($dp_x, dp_y$)
        \item Binned Slice Studies ($x_F$ and $p_T$ bins)
    \end{itemize}
\end{frame}

\begin{frame}{Files Used \& Kinematic Bins}
    \textbf{Files Used:}
    \begin{itemize}
        \item LH2 and LD2 Drell-Yan Monte Carlo (M027\_S001)
        \item \textbf{Thrown:} \texttt{mc\_drellyan\_*\_4pi\_pTxFweight\_v2.root}
        \item \textbf{Accepted:} \texttt{mc\_drellyan\_*\_clean\_occ\_pTxFweight\_v2.root}
    \end{itemize}
    \vspace{0.2cm}
    \textbf{Acceptance Calculation:}
    \begin{equation*}
        \text{Acceptance} = \frac{\sum w_{\text{accepted}}}{\sum w_{\text{thrown}}} \quad \text{where } w = \text{ReWeight}
    \end{equation*}
    \vspace{0.1cm}
    \textbf{Kinematic Bin Boundaries:}
    \begin{itemize}
        \item \textbf{Mass [GeV]:} [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7]
        \item \textbf{$x_F$:} [0.0, 0.05, 0.10, ..., 0.80]
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

    # ==========================================
    # 1. Invariant Mass Studies
    # ==========================================
    latex_content += r"\section{Invariant Mass Plots}" + "\n"
    latex_content += make_frame_2cols("Overlay \\& Ratio: Mass (Integrated)", 
                                      "acceptance_overlay_Integrated_mass.pdf", 
                                      "acceptance_ratio_Integrated_mass.pdf")
    
    latex_content += make_frame_1col("Invariant Mass Yields \\& Ratio (All $x_F$, $p_T$)", "Split_Mass_All_xF_pT.pdf")
    latex_content += make_frame_1col("Invariant Mass Yields \\& Ratio ($0.0 < x_F \\leq 0.4$)", "Split_Mass_0.0_xF_0.4.pdf")
    latex_content += make_frame_1col("Invariant Mass Yields \\& Ratio ($0.4 < x_F \\leq 0.8$)", "Split_Mass_0.4_xF_0.8.pdf")


    # ==========================================
    # 2. xF Studies
    # ==========================================
    latex_content += r"\section{$x_F$ Plots}" + "\n"
    latex_content += make_frame_2cols("Overlay \\& Ratio: $x_F$ (Integrated)", 
                                      "acceptance_overlay_Integrated_xF.pdf", 
                                      "acceptance_ratio_Integrated_xF.pdf")
    
    # xF Yields
    latex_content += make_frame_1col("$x_F$ Yields \\& Ratio (All Mass, $p_T$)", "Split_xF_All_Mass_pT.pdf")
    latex_content += make_frame_1col("$x_F$ Yields \\& Ratio ($4.2 < \\text{Mass} \\leq 5.5$)", "Split_xF_4.2_Mass_5.5.pdf")
    latex_content += make_frame_1col("$x_F$ Yields \\& Ratio ($5.5 < \\text{Mass} \\leq 8.8$)", "Split_xF_5.5_Mass_8.8.pdf")

    # xF Acceptances (New Additions)
    latex_content += make_frame_1col("$x_F$ Acceptance \\& Ratio (All Mass, $p_T$)", "Acceptance_xF_All_Mass_pT.pdf")
    latex_content += make_frame_1col("$x_F$ Acceptance \\& Ratio ($4.2 < \\text{Mass} < 5.5$)", "Acceptance_xF_4.2_Mass_5.5.pdf")
    latex_content += make_frame_1col("$x_F$ Acceptance \\& Ratio ($5.5 < \\text{Mass} < 8.7$)", "Acceptance_xF_5.5_Mass_8.7.pdf")


    # ==========================================
    # 3. pT Studies (Yields & Acceptances)
    # ==========================================
    latex_content += r"\section{$p_T$ Plots}" + "\n"
    latex_content += make_frame_2cols("Overlay \\& Ratio: $p_T$ (Integrated)", 
                                      "acceptance_overlay_Integrated_pT.pdf", 
                                      "acceptance_ratio_Integrated_pT.pdf")
    
    # Fine Binning Yields
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (Fine Bins, All Mass, $x_F$)", "Split_pT_Fine_All_Mass_xF.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (Fine Bins, $0.0 < x_F < 0.4$)", "Split_pT_Fine_0.0_xF_0.4.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (Fine Bins, $0.4 < x_F < 0.8$)", "Split_pT_Fine_0.4_xF_0.8.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (Fine Bins, $4.2 < \\text{Mass} < 5.5$)", "Split_pT_Fine_4.2_Mass_5.5.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (Fine Bins, $5.5 < \\text{Mass} < 8.7$)", "Split_pT_Fine_5.5_Mass_8.7.pdf")

    # User Binning Yields
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (User Bins, All Mass, $x_F$)", "Split_pT_User_All_Mass_xF.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (User Bins, $0.0 < x_F < 0.4$)", "Split_pT_User_0.0_xF_0.4.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (User Bins, $0.4 < x_F < 0.8$)", "Split_pT_User_0.4_xF_0.8.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (User Bins, $4.2 < \\text{Mass} < 5.5$)", "Split_pT_User_4.2_Mass_5.5.pdf")
    latex_content += make_frame_1col("$p_T$ Yields \\& Ratio (User Bins, $5.5 < \\text{Mass} < 8.7$)", "Split_pT_User_5.5_Mass_8.7.pdf")

    # Acceptances
    latex_content += make_frame_1col("$p_T$ Acceptance \\& Ratio (User Bins, All Mass, $x_F$)", "Acceptance_pT_All_Mass_xF.pdf")
    latex_content += make_frame_1col("$p_T$ Acceptance \\& Ratio (User Bins, $0.0 < x_F < 0.4$)", "Acceptance_pT_0.0_xF_0.4.pdf")
    latex_content += make_frame_1col("$p_T$ Acceptance \\& Ratio (User Bins, $0.4 < x_F < 0.8$)", "Acceptance_pT_0.4_xF_0.8.pdf")
    latex_content += make_frame_1col("$p_T$ Acceptance \\& Ratio (User Bins, $4.2 < \\text{Mass} < 5.5$)", "Acceptance_pT_4.2_Mass_5.5.pdf")
    latex_content += make_frame_1col("$p_T$ Acceptance \\& Ratio (User Bins, $5.5 < \\text{Mass} < 8.7$)", "Acceptance_pT_5.5_Mass_8.7.pdf")


    # ==========================================
    # 4. Momentum Kinematics (dp_x, dp_y)
    # ==========================================
    latex_content += r"\section{Track Kinematics ($dp_x, dp_y$)}" + "\n"
    
    latex_content += make_frame_2cols("Kinematic Yields ($4\pi$ Thrown): Integrated", 
                                      "yield_th_dpx_integrated.pdf", 
                                      "yield_th_dpy_integrated.pdf")
    latex_content += make_frame_2cols("Kinematic Yields (Clean Accepted): Integrated", 
                                      "yield_ac_dpx_integrated.pdf", 
                                      "yield_ac_dpy_integrated.pdf")
                                      
    latex_content += make_frame_2cols("Kinematic Yields Squared ($4\pi$ Thrown): Integrated", 
                                      "yield_th_dpx2_integrated.pdf", 
                                      "yield_th_dpy2_integrated.pdf")
    latex_content += make_frame_2cols("Kinematic Yields Squared (Clean Accepted): Integrated", 
                                      "yield_ac_dpx2_integrated.pdf", 
                                      "yield_ac_dpy2_integrated.pdf")


    # ==========================================
    # 5. Binned Slice Studies (xF Slices)
    # ==========================================
    latex_content += r"\section{Binned Studies: $x_F$ Slices}" + "\n"
    for i in range(num_xf_bins):
        prefix = f"xF_bin{i}"
        title = f"$x_F$ Bin {i}"
        
        # 4pi and Clean Yield differences (Separated Canvases)
        latex_content += make_frame_2cols(f"Kinematic Yields ($4\pi$ Thrown): {title}", 
                                          f"yield_th_dpx_{prefix}.pdf", 
                                          f"yield_th_dpy_{prefix}.pdf")
        latex_content += make_frame_2cols(f"Kinematic Yields (Clean Accepted): {title}", 
                                          f"yield_ac_dpx_{prefix}.pdf", 
                                          f"yield_ac_dpy_{prefix}.pdf")
                                          
        latex_content += make_frame_2cols(f"Kinematic Yields Squared ($4\pi$ Thrown): {title}", 
                                          f"yield_th_dpx2_{prefix}.pdf", 
                                          f"yield_th_dpy2_{prefix}.pdf")
        latex_content += make_frame_2cols(f"Kinematic Yields Squared (Clean Accepted): {title}", 
                                          f"yield_ac_dpx2_{prefix}.pdf", 
                                          f"yield_ac_dpy2_{prefix}.pdf")
        
        latex_content += make_frame_2cols(f"Overlay \\& Ratio: {title}", 
                                          f"acceptance_overlay_{prefix}.pdf", 
                                          f"acceptance_ratio_{prefix}.pdf")

    # ==========================================
    # 6. Binned Slice Studies (pT Slices)
    # ==========================================
    latex_content += r"\section{Binned Studies: $p_T$ Slices}" + "\n"
    for i in range(num_pt_bins):
        prefix = f"pT_bin{i}"
        title = f"$p_T$ Bin {i}"
        latex_content += make_frame_2cols(f"Overlay \\& Ratio: {title}", 
                                          f"acceptance_overlay_{prefix}.pdf", 
                                          f"acceptance_ratio_{prefix}.pdf")

    # Close the document
    latex_content += r"\end{document}"

    # Write out to file
    with open(output_tex, "w") as f:
        f.write(latex_content)

    print(f"Successfully generated {output_tex}")

if __name__ == "__main__":
    main()