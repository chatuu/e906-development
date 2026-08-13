import os
import re
from datetime import date

def generate_latex():
    # Configuration
    base_dir = "/root/github/e906-development/src/xsec_pT/RS57-70"
    backup_dir = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS57-70"
    acceptance_dir = "/root/github/e906-development/src/AcceptanceCorrection"
    pt2_dir = "/root/github/e906-development/src/xsec_pT_squard/RS57-70"
    xf_dir = "/root/github/e906-development/src/xsec_pT_xF/"
    ratio_comp_dir = "/root/github/e906-development/src/xsec_pT_xsec_comparison_to_RS67/ratio_plots/"

    author = "Chatura Kuruppu"
    footer_center_text = r"Addressing Comments Received for DocDB 11535"
    title = r"Measurement of Absolute Single Differential Cross-Section in $p_T$ Bins using Runs 2 and 3 Data\\ (Addressing comments received for DocDB 11535)"
    today_str = date.today().strftime("%Y-%m-%d")

    # Content LaTeX Strings
    overview_latex = r"""
\begin{itemize}
    \setlength{\itemsep}{0.5em}
    \item Events from runs 2 and 3
    \item Kinematic Phase Space
    \item Event Selection Criteria
    \item Comments received on previous talk:
    \begin{itemize} 
        \item $d\sigma/dp^{2}_T$ Distribution for different $x_F$ bins
        \item $d\sigma/dp^{2}_T$ Distribution and ratio plot
        \item Cross-section ratio as a function of $x_F$
    \end{itemize}
    \item Cross-section ratio from different roadsets to RS67
    \item Next Steps
\end{itemize}
"""

    inputs_latex = r"""
\textbf{Note:} Currently using all the runs saved in runs 2 and 3.
\vspace{0.3cm}
\tiny\\
\textbf{Harsha's ROOT files:} \\
\texttt{/seaquest/users/harshaka/e906\_project/e906-root-ana/work\_gpvm/scripts/results\_runFinalTree/merged\_results\_e906\_mixing/:}
\begin{multicols}{2}
\begin{itemize}
    \item \texttt{merged\_RS57\_LH2\_1\_1138.root}
    \item \texttt{merged\_RS57\_LD2\_3\_1138.root}
    \item \texttt{merged\_RS57\_Empty\_2\_1138.root}
    \item \texttt{merged\_RS59\_LH2\_1\_465.root}
    \item \texttt{merged\_RS59\_LD2\_3\_466.root}
    \item \texttt{merged\_RS59\_Empty\_2\_466.root}
    \item \texttt{merged\_RS62\_LH2\_1\_1234.root}
    \item \texttt{merged\_RS62\_LD2\_3\_1237.root}
    \item \texttt{merged\_RS62\_Empty\_2\_1234.root}
    \item \texttt{merged\_RS70\_LH2\_1\_264.root}
    \item \texttt{merged\_RS70\_LD2\_3\_266.root}
    \item \texttt{merged\_RS70\_Empty\_2\_267.root}
\end{itemize}
\end{multicols}
\textbf{Abi's ROOT files:} \\
\texttt{/seaquest/users/apun/e906\_projects/rs67\_merged\_files/:}
\begin{multicols}{1}
\begin{itemize}
    \item \texttt{merged\_RS67\_3089LH2.root}
    \item \texttt{merged\_RS67\_3089LD2.root}
    \item \texttt{merged\_RS67\_3089flask.root}
\end{itemize}
\end{multicols}

\normalsize
"""

    kinematic_bins_latex = r"""
\textbf{Kinematic Binning Definition:}
\begin{itemize}
    \item \textbf{Mass Bins (GeV):} [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.8]
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
\vspace{0.1cm}
\textbf{Track Quality \& Acceptance:}
\begin{itemize}
    \item $\chi^2_{\text{target}} < 15$, $\chi^2_{\text{dimuon}} < 18$, and $\chi^2/(\text{nHits}-5) < 12$.
    \item Hit requirements: nHits $> 13$ per track, $\text{nHits}_1+\text{nHits}_2 > 29$, St1 total $> 8$.
    \item Track vertex longitudinal limits: $-320 < z_v < -5$.
    \item A dynamic beam offset correction is applied to y-coordinates: offset is $1.6$ for runID $\ge 11000$ and $0.4$ otherwise.
\end{itemize}
\vspace{0.1cm}
\textbf{Detector Occupancy Limits:}
\begin{itemize}
    \item Drift chamber hits: $D1 < 400$, $D2 < 400$, $D3 < 400$, and total $D1+D2+D3 < 1000$.
    \item $20 < D1 < 385$ (global reconstruction efficiency curve defined in this range)
\end{itemize}
\textbf{ONLY RS62,  $runID > 11500$ applied to ensure correct target positions of dimuons}
\vspace{0.1cm}
\textbf{Note:} Currently using the same set of Chuck Cuts defined in \href{https://seaquest-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=2111}{DocDB 2111-V42}
"""

    pot_latex = r"""
\begin{table}[H]
    \centering
    \renewcommand{\arraystretch}{1.3}
    \begin{tabular}{|c|c|c|c|}
        \hline
        \textbf{Roadset} & \textbf{LH2} & \textbf{LD2} & \textbf{Flask} \\ \hline
        57 & $3.533324 \times 10^{16}$ & $1.768358 \times 10^{16}$ & $3.918550 \times 10^{15}$ \\ \hline
        59 & $9.365986 \times 10^{15}$ & $4.319952 \times 10^{15}$ & $1.010350 \times 10^{15}$ \\ \hline
        62 & $5.281659 \times 10^{16}$ & $2.382899 \times 10^{16}$ & $1.097456 \times 10^{16}$ \\ \hline
        67 & $1.611435 \times 10^{17}$ & $7.694541 \times 10^{16}$ & $3.662417 \times 10^{16}$ \\ \hline
        70 & $1.785745 \times 10^{16}$ & $8.752588 \times 10^{15}$ & $3.841280 \times 10^{15}$ \\ \hline
    \end{tabular}
\end{table}
\vspace{0.3cm}
\footnotesize
\textbf{Note:} POT Values from runs mixed by Harsha in \href{https://seaquest-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=11524}{DocDB 11524}.
\normalsize
"""

    docdb_comments_latex = r"""
\begin{itemize}
    \item Generate single differential cross-section as a function of $p^{2}_{T}$
    \item Calculate cross-section ratio
    \item Generate single differential cross-section as a function of $p_{T}$ for different $x_{F}$ bins
    \item Calculate cross-section ratio
    \item Generate cross-section ratio as a fuction of different $x_F$ bins
\end{itemize}
"""

    summary_table_latex = r"""
\begin{table}[H]
    \centering
    \renewcommand{\arraystretch}{1.5}
    \resizebox{\textwidth}{!}{
    \begin{tabular}{|l|c|c|c|c|c|}
        \hline
        \textbf{Roadset} & \textbf{57} & \textbf{59} & \textbf{62} & \textbf{70} & \textbf{57-70} \\ \hline
        Cross-Section Ratio: Roadset/RS67 (LH2) & $0.754 \pm 0.031$ & $0.729 \pm 0.047$ & $0.994 \pm 0.028$ & $1.035 \pm 0.040$ & $0.956 \pm 0.020$ \\ \hline
        Cross-Section Ratio: Roadset/RS67 (LD2) & $0.731 \pm 0.020$ & $0.677 \pm 0.033$ & $0.973 \pm 0.021$ & $0.986 \pm 0.029$ & $0.946 \pm 0.016$ \\ \hline
    \end{tabular}
    }
\end{table}
\vspace{0.2cm}
\footnotesize
\textbf{Note:} These cross-section ratios are consistent with Harsha's numbers mentioned in \href{https://seaquest-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=11539&filename=dy_csr_runs_2_3.pdf&version=2}{DocDB: 11539}
\normalsize
"""

    next_steps_latex = r"""
\begin{itemize}
    \setlength{\itemsep}{1em}
    \item Need to work on unfolding (work on progress).
    \item Planning to present preliminary plots during upcoming APS DPF meeting.
    \item Analysis note will be written for this Run 2-3 release, as a preliminary result.
\end{itemize}
"""
    
    reco_efficiency_explanation_latex = r"Reco efficiency calculated using D1 occupancy and global efficiency curve."
    hodo_efficiency_explanation_latex = r"Hodoscope efficiencies propagated from latest paddle efficiency tables (DocDB 11467)."
    final_efficiency_explanation_latex = r"Final efficiency determined as the product of reco and hodo components."
    epsilon_signal_formula_latex = r"Calculation of signal efficiency for target background subtraction."
    signal_candidates_latex = r"Determination of corrected yield for LH2 and LD2 signal candidates."
    epsilon_signal_latex = r"\begin{equation*} \epsilon_{\text{signal}} = \frac{\langle \epsilon \rangle_{\text{total}} Y_{\text{total}} - \langle \epsilon \rangle_{\text{mix}} Y_{\text{mix}}}{Y_{\text{total}} - Y_{\text{mix}}} \end{equation*}"
    yield_distributions_latex = r"Yield calculations for LH2, LD2, and Empty Flask targets."

    # Plot sequence mapping
    plot_names = [
        "Y_total_LH2", "Y_total_LD2", "Y_total_Flask",
        "Y_mix_LH2", "Y_mix_LD2", "Y_mix_Flask",
        "E_total_reco_LH2", "E_total_reco_LD2", "E_total_reco_Flask",
        "E_mix_reco_LH2", "E_mix_reco_LD2", "E_mix_reco_Flask",
        "E_total_hodo_LH2", "E_total_hodo_LD2", "E_total_hodo_Flask",
        "E_mix_hodo_LH2", "E_mix_hodo_LD2", "E_mix_hodo_Flask",
        "E_total_final_LH2", "E_total_final_LD2", "E_total_final_Flask",
        "E_mix_final_LH2", "E_mix_final_LD2", "E_mix_final_Flask",
        "E_final_signal_LH2", "E_final_signal_LD2", "E_final_signal_Flask",
        "Y_corrected_LH2", "Y_corrected_LD2", "Y_corrected_Flask",
        "Y_corrected_Subtracted_LH2", "Y_corrected_Subtracted_LD2",
        "CrossSection_LH2_geom_vs_pT_with_logo",
        "CrossSection_LD2_geom_vs_pT_with_logo",
        "CrossSection_Ratio_pd_2pp_vs_pT_geom_with_logo",
        "Combined_XSec_Ratio_vs_pT_geom",
        "Combined_XSec_Ratio_vs_pT_true_pt"
    ]
    ordered_files = [f"{name}.pdf" for name in plot_names]

    pt2_plots = [
        "CrossSection_LH2_geom_pT2_with_logo.pdf",
        "CrossSection_LD2_geom_pT2_with_logo.pdf",
        "CrossSection_Ratio_pd_2pp_geom_pT2_with_logo.pdf",
        "Combined_XSec_Ratio_geom_pT2.pdf"
    ]

    ratio_pairs = [
        ("Ratio_RS57_vs_RS67_LH2.pdf", "Ratio_RS57_vs_RS67_LD2.pdf"),
        ("Ratio_RS59_vs_RS67_LH2.pdf", "Ratio_RS59_vs_RS67_LD2.pdf"),
        ("Ratio_RS62_vs_RS67_LH2.pdf", "Ratio_RS62_vs_RS67_LD2.pdf"),
        ("Ratio_RS70_vs_RS67_LH2.pdf", "Ratio_RS70_vs_RS67_LD2.pdf"),
        ("Ratio_RS57-70_vs_RS67_LH2.pdf", "Ratio_RS57-70_vs_RS67_LD2.pdf")
    ]

    def tex_escape(text):
        return text.replace("_", r"\_")
    
    def clean_title(text):
        return text.replace(".pdf", "").replace("_", " ")

    beamer_file = "slides.tex"
    with open(beamer_file, "w") as f:
        # Preamble and Title Slide
        f.write(r"""\documentclass[aspectratio=169]{beamer}
\usetheme{Madrid}\usepackage{graphicx}\usepackage{multicol}\usepackage{hyperref}
\setbeamertemplate{footline}{\leavevmode\hbox{\begin{beamercolorbox}[wd=.333333\paperwidth,ht=2.25ex,dp=1ex,center]{author in head/foot}%s\end{beamercolorbox}\begin{beamercolorbox}[wd=.333333\paperwidth,ht=2.25ex,dp=1ex,center]{title in head/foot}%s\end{beamercolorbox}\begin{beamercolorbox}[wd=.333333\paperwidth,ht=2.25ex,dp=1ex,left]{date in head/foot}\hspace*{\fill}%s\hspace*{\fill}\rlap{\insertframenumber/\inserttotalframenumber}\hspace*{2ex}\end{beamercolorbox}}\vskip0pt}
\graphicspath{{%s/}{%s/}{%s/}{%s/}{%s/}{%s/}}

\title{%s}
\author{%s}
\institute{New Mexico State University\\ \vspace{0.1cm} SeaQuest Experiment (E906)}
\date{\today}

\begin{document}
\begin{frame}\titlepage\end{frame}
\begin{frame}{Overview}%s\end{frame}
\begin{frame}{Events from runs 2 and 3 (Input Files Used)}%s\end{frame}
\begin{frame}{Kinematic Phase Space}%s\end{frame}
\begin{frame}{Event Selection Criteria (Chuck Cuts)}%s\end{frame}
\begin{frame}{POT Values Used}%s\end{frame}
\begin{frame}{Comments Received for DocDB: 11535}%s\end{frame}
""" % (author, footer_center_text, today_str, base_dir, backup_dir, acceptance_dir, pt2_dir, xf_dir, ratio_comp_dir, 
       title, author, overview_latex, inputs_latex, kinematic_bins_latex, event_selection_latex, pot_latex, docdb_comments_latex))

        # Main Plot Sequence
        for idx, file in enumerate(ordered_files):
            if idx >= 32 and idx <= 35:
                frame_title = clean_title(file)
                f.write(r"\begin{frame}{%s}\begin{center}\includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{%s}\end{center}\end{frame}" % (tex_escape(frame_title), file))
                
                if idx == 35:
                    for pfile in pt2_plots:
                        f.write(r"\begin{frame}{%s}\begin{center}\includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{%s}\end{center}\end{frame}" % (clean_title(pfile), pfile))
                    
                    for i in [0] + list(range(1, 16)):
                        low_bound = i * 0.05
                        high_bound = (i + 1) * 0.05
                        xf_title = f"Single Differential Cross-Sections ${low_bound:0.2f} < xF < {high_bound:0.2f}$"
                        xfile = f"Combined_XSec_Ratio_xF_{i}_vs_pT_geom.pdf"
                        f.write(r"\begin{frame}{%s}\begin{center}\includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{%s}\end{center}\end{frame}" % (xf_title, xfile))
                    
                    f.write(r"\begin{frame}{Cross-Section ratio Vs xF}\begin{center}\includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{CrossSection_Ratio_pd_2pp_vs_xF_geom_with_logo.pdf}\end{center}\end{frame}")
                    
                    for lh2, ld2 in ratio_pairs:
                        match = re.search(r"Ratio_(RS[0-9-]+)_vs", lh2)
                        rs_num = match.group(1) if match else "Unknown"
                        rs_title = f"Cross-Section Ratio {rs_num}/RS67"
                        f.write(r"\begin{frame}{%s}\begin{columns}\column{0.5\textwidth}\centering \textbf{LH2 Target}\\\includegraphics[width=\textwidth]{%s}\column{0.5\textwidth}\centering \textbf{LD2 Target}\\\includegraphics[width=\textwidth]{%s}\end{columns}\end{frame}" % (rs_title, lh2, ld2))
                    
                    f.write(r"\begin{frame}{Summary}%s\end{frame}" % summary_table_latex)

        f.write(r"\begin{frame}{Next Steps}%s\end{frame}" % next_steps_latex)

        f.write(r"\appendix \begin{frame}\centering\Huge Backup Slides\end{frame}")
        f.write(r"\begin{frame}{Epsilon Signal}%s\end{frame}\begin{frame}{Yield Distributions}%s\end{frame}" % (epsilon_signal_latex, yield_distributions_latex))

        for idx, file in enumerate(ordered_files):
            if idx < 32:
                if os.path.exists(os.path.join(base_dir, file)):
                    frame_title = clean_title(file)
                    f.write(r"\begin{frame}{%s}\begin{center}\includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{%s}\end{center}\end{frame}" % (tex_escape(frame_title), file))
                if idx == 5: f.write(r"\begin{frame}{Reco Explanation}%s\end{frame}" % reco_efficiency_explanation_latex)
                if idx == 11: f.write(r"\begin{frame}{Hodo Explanation}%s\end{frame}" % hodo_efficiency_explanation_latex)

        file_36 = ordered_files[36]
        f.write(r"\begin{frame}{%s}\begin{center}\includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{%s}\end{center}\end{frame}" % (tex_escape(clean_title(file_36)), file_36))
        f.write(r"\end{document}")

    print("Successfully generated slides.tex")

if __name__ == "__main__":
    generate_latex()