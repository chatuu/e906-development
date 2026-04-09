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

def main():
    output_tex = "acceptance_slides.tex"

    # Number of bins based on your numpy arrays
    num_xf_bins = 16 # np.arange(0.0, 0.85, 0.05) yields 17 edges = 16 bins
    num_pt_bins = 7  # [0., 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8] = 8 edges = 7 bins

    # LaTeX document header using a 16:9 widescreen aspect ratio
    latex_content = r"""\documentclass[aspectratio=169]{beamer}
\usepackage{graphicx}
\usepackage{amsmath}

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
        \item Files Used
        \item Event Selection Criteria
        \item Understanding the Acceptance Ratio Discrepancy
        \item Plots Generated
        \begin{itemize}
            \item Fully Integrated 1D Plots
            \item $x_F$ Binned Plots
            \item $p_T$ Binned Plots
        \end{itemize}
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
            \end{itemize}
        \end{column}
    \end{columns}
\end{frame}

% \begin{frame}{Why is the LH2/LD2 Ratio $\neq$ 1 for Integrated Bins?}
%     \textbf{The Phase-Space Convolution Trap:}
%     \vspace{0.2cm}
%   
%     A 1D projected acceptance is not just the detector geometry; it is mathematically weighted by the generated physics distributions ($N_{\text{gen}}$):
%     \begin{equation*}
%         \epsilon(p_T) = \frac{\int \int A(x_F, m, p_T) N_{\text{gen}}(x_F, m, p_T) \,dx_F \,dm}{\int \int N_{\text{gen}}(x_F, m, p_T) \,dx_F \,dm}
%     \end{equation*}
%
%     \vspace{0.2cm}
%     \begin{itemize}
%         \item \textbf{Different Physics:} Due to target isospin and PDFs, $N_{\text{gen}}^{\text{LH2}} \neq N_{\text{gen}}^{\text{LD2}}$. 
%         \item \textbf{Steep Acceptance Gradients:} The spectrometer acceptance changes drastically across $x_F$ and Mass.
%         \item \textbf{The Result:} When integrating out unbinned variables (like in a $p_T$ plot), the different LH2 and LD2 kinematic shapes sample the spectrometer's geometric acceptance differently, shifting the "average" ratio away from 1.0.
%     \end{itemize}
% \end{frame}
"""

    # ==========================================
    # 1. Fully Integrated 1D Plots
    # ==========================================
    latex_content += r"\section{Fully Integrated Plots}" + "\n"
    
    # 4pi and Clean Yield differences (Separated Canvases)
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

    integrated_vars = [
        ("pT", "$p_T$ (Integrated over $x_F$, Mass)"),
        ("xF", "$x_F$ (Integrated over $p_T$, Mass)"),
        ("mass", "Mass (Integrated over $x_F$, $p_T$)")
    ]

    for var_name, var_title in integrated_vars:
        prefix = f"Integrated_{var_name}"
        latex_content += make_frame_2cols(f"Overlay \\& Ratio: {var_title}", 
                                          f"acceptance_overlay_{prefix}.pdf", 
                                          f"acceptance_ratio_{prefix}.pdf")

    # ==========================================
    # 2. xF Binned Plots
    # ==========================================
    latex_content += r"\section{$x_F$ Binned Plots}" + "\n"
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
    # 3. pT Binned Plots
    # ==========================================
    latex_content += r"\section{$p_T$ Binned Plots}" + "\n"
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