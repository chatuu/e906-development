import os

# ==========================================
# Configuration
# ==========================================
old_dir = "../StepByStepPlots"
new_dir = "."

targets = ["LH2", "LD2", "Flask"]
yield_types = ["total", "mix"]
eff_types = [
    ("hodo", "Hodoscope Efficiency Corrections"),
    ("reco", "Reconstruction Efficiency Corrections"),
    ("final", "Final Efficiency Corrections")
]

output_report_file = "efficiency_corrections_comparison.tex"
output_slides_file = "efficiency_corrections_slides.tex"

# ==========================================
# Mathematical Explanations (LaTeX Blocks)
# ==========================================
math_section_report = r"""
\section{Methodology: Efficiency and Uncertainty Calculations}

This section details the updated statistical methodologies used to evaluate the average efficiencies and propagate their associated uncertainties within each kinematic mass and $x_F$ bin.

\subsection{Reconstruction Efficiency and Correlated Uncertainties}
The reconstruction efficiency $E_{{\rm reco}, i}$ for a single event $i$ is interpolated from a global tracking efficiency curve evaluated at the chamber occupancy parameter, $D_1$. Because multiple events within the same kinematic bin frequently fall between the same $D_1$ occupancy, their interpolated uncertainties $\sigma_i$ are mathematically correlated.

To account for this, we construct an $N \times N$ covariance matrix for the $N$ events within a given kinematic bin:
\begin{equation}
    Cov_{i,j} = \rho_{i,j} \sigma_i \sigma_j
\end{equation}
where the correlation coefficient $\rho_{i,j}$ is determined by the shared occupancy:
\begin{itemize}
    \item $\rho_{i,j} = 1$ if event $i$ and event $j$ share both bounding occupancy (same $D_1$ bin).
    \item $\rho_{i,j} = 1$ if event $i$ and event $j$ share exactly one bounding interpolation node (adjacent $D_1$ bins).
    \item $\rho_{i,j} = 0$ otherwise (assumed strictly independent).
\end{itemize}

The mean reconstruction efficiency for the kinematic bin is evaluated as the arithmetic mean of the event-by-event efficiencies:
\begin{equation}
    \bar{E}_{\rm reco} = \frac{1}{N} \sum_{i=1}^{N} E_{{\rm reco}, i}
\end{equation}

By accounting for the off-diagonal covariance elements, the correctly propagated, correlated uncertainty on the mean efficiency is:
\begin{equation}
    \sigma_{\bar{E}_{\rm reco}} = \frac{1}{N} \sqrt{ \sum_{i=1}^{N} \sum_{j=1}^{N} Cov_{i,j} }
\end{equation}

\subsection{Final Efficiency Evaluation}
The mean hodoscope efficiency $\bar{E}_{\rm hodo}$ is derived similarly. However, because its event-by-event evaluations are independent of the global tracking $D_1$ curve, its uncertainty is evaluated assuming strictly independent events (a purely diagonal covariance matrix):
\begin{equation}
    \sigma_{\bar{E}_{\rm hodo}} = \frac{1}{N} \sqrt{ \sum_{i=1}^{N} \sigma_{{\rm hodo}, i}^2 }
\end{equation}

Previously, the final efficiency was calculated by multiplying $E_{\rm reco} \times E_{\rm hodo}$ on an event-by-event basis prior to averaging. To properly isolate and utilize the independent error characteristics of the tracking and hodoscope responses, the current methodology calculates the final efficiency as the product of the independently evaluated bin averages:
\begin{equation}
    \bar{E}_{\rm final} = \bar{E}_{\rm reco} \times \bar{E}_{\rm hodo}
\end{equation}

The final propagated uncertainty on this product is evaluated using standard error propagation for independent variables:
\begin{equation}
    \sigma_{\bar{E}_{\rm final}} = \sqrt{ \left(\bar{E}_{\rm hodo} \cdot \sigma_{\bar{E}_{\rm reco}}\right)^2 + \left(\bar{E}_{\rm reco} \cdot \sigma_{\bar{E}_{\rm hodo}}\right)^2 }
\end{equation}
\clearpage
"""

math_slides = r"""
\section{Methodology}

\begin{frame}{Methodology}{Reconstruction Efficiency and Covariance}
    \textbf{1. Mean Reconstruction Efficiency:}
    Calculated as the arithmetic mean of event-by-event interpolated efficiencies for $N$ events in a bin:
    \begin{equation*}
        \bar{E}_{\rm reco} = \frac{1}{N} \sum_{i=1}^{N} E_{{\rm reco}, i}
    \end{equation*}
    
    \vspace{0.3cm}
    \textbf{2. Correlated Error Propagation:}
    Because efficiencies are interpolated from the same $D_1$ nodes, event uncertainties ($\sigma_i$) are correlated. We construct a covariance matrix:
    \begin{equation*}
        Cov_{i,j} = \rho_{i,j} \sigma_i \sigma_j
    \end{equation*}
    \begin{itemize}
        \item $\rho = 1$ if events share identical or adjacent $D_1$ occupancy.
        \item $\rho = 0$ if events use entirely different nodes.
    \end{itemize}
    \begin{equation*}
        \sigma_{\bar{E}_{\rm reco}} = \frac{1}{N} \sqrt{ \sum_{i=1}^{N} \sum_{j=1}^{N} Cov_{i,j} }
    \end{equation*}
\end{frame}

\begin{frame}{Methodology}{Final Efficiency Corrections}
    \textbf{1. Hodoscope Mean Efficiency:}
    Evaluated with purely independent, standard weighted uncertainty:
    \begin{equation*}
        \sigma_{\bar{E}_{\rm hodo}} = \frac{1}{N} \sqrt{ \sum_{i=1}^{N} \sigma_{{\rm hodo}, i}^2 }
    \end{equation*}
    
    \vspace{0.3cm}
    \textbf{2. Final Bin Efficiency:}
    Computed as the product of the independent bin averages:
    \begin{equation*}
        \bar{E}_{\rm final} = \bar{E}_{\rm reco} \times \bar{E}_{\rm hodo}
    \end{equation*}
    
    \vspace{0.3cm}
    \textbf{3. Final Propagated Uncertainty:}
    \begin{equation*}
        \sigma_{\bar{E}_{\rm final}} = \sqrt{ \left(\bar{E}_{\rm hodo} \cdot \sigma_{\bar{E}_{\rm reco}}\right)^2 + \left(\bar{E}_{\rm reco} \cdot \sigma_{\bar{E}_{\rm hodo}}\right)^2 }
    \end{equation*}
\end{frame}
"""

# ==========================================
# 1. LaTeX Report Generation (Article)
# ==========================================
def generate_report():
    tex_content = r"""\documentclass[12pt, a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{geometry}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{float}
\usepackage{hyperref}

\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,      
    urlcolor=cyan,
    pdftitle={Efficiency Corrections Comparison},
}

\geometry{top=1in, bottom=1in, left=0.5in, right=0.5in}

\title{Efficiency Corrections Comparison: Current vs. Previous Results}
\author{Chatura Kuruppu}
\date{\today}

\begin{document}
\maketitle

\tableofcontents
\clearpage

\listoffigures
\clearpage
"""
    # Inject Math Methodology Section
    tex_content += math_section_report

    for target in targets:
        tex_content += f"\\section{{{target} Efficiency Corrections}}\n\n"
        
        for y_type in yield_types:
            tex_content += f"\\subsection{{{y_type.capitalize()} Yield}}\n\n"
            
            for e_type_prefix, e_type_title in eff_types:
                old_file = f"{old_dir}/E_{y_type}_{e_type_prefix}_{target}.pdf"
                new_file = f"{new_dir}/E_{y_type}_{e_type_prefix}_{target}.pdf"
                
                tex_content += r"""
\begin{figure}[H]
    \centering
    \begin{subfigure}{0.49\textwidth}
        \centering
        \includegraphics[width=\linewidth]{""" + old_file + r"""}
        \caption{Previous}
    \end{subfigure}\hfill
    \begin{subfigure}{0.49\textwidth}
        \centering
        \includegraphics[width=\linewidth]{""" + new_file + r"""}
        \caption{Current}
    \end{subfigure}
    \caption{""" + f"{target} {y_type.capitalize()} - {e_type_title}" + r"""}
\end{figure}
"""
            tex_content += "\\clearpage\n"

        tex_content += f"\\subsection{{Signal Efficiency Corrections}}\n\n"
        
        old_sig_file = f"{old_dir}/E_final_signal_{target}.pdf"
        new_sig_file = f"{new_dir}/E_final_signal_{target}.pdf"
        
        tex_content += r"""
\begin{figure}[H]
    \centering
    \begin{subfigure}{0.49\textwidth}
        \centering
        \includegraphics[width=\linewidth]{""" + old_sig_file + r"""}
        \caption{Previous}
    \end{subfigure}\hfill
    \begin{subfigure}{0.49\textwidth}
        \centering
        \includegraphics[width=\linewidth]{""" + new_sig_file + r"""}
        \caption{Current}
    \end{subfigure}
    \caption{""" + f"{target} - Signal Efficiency Corrections" + r"""}
\end{figure}
\clearpage
"""

    tex_content += "\\end{document}\n"

    with open(output_report_file, "w") as f:
        f.write(tex_content)
    print(f"LaTeX report '{output_report_file}' generated successfully.")


# ==========================================
# 2. LaTeX Slides Generation (Beamer)
# ==========================================
def generate_slides():
    tex_content = r"""\documentclass[aspectratio=169]{beamer}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{subcaption}

% NMSU Brand Colors
\definecolor{AggieCrimson}{HTML}{8C0B42}
\definecolor{DesertMonsoon}{HTML}{6D6E71}

% Clean, professional beamer theme customized for NMSU
\usetheme{Madrid}
\usecolortheme[named=AggieCrimson]{structure}
\setbeamercolor{palette primary}{bg=AggieCrimson,fg=white}
\setbeamercolor{palette secondary}{bg=DesertMonsoon,fg=white}
\setbeamercolor{palette tertiary}{bg=AggieCrimson,fg=white}
\setbeamercolor{palette quaternary}{bg=AggieCrimson,fg=white}
\setbeamercolor{titlelike}{parent=palette primary,bg=AggieCrimson,fg=white}
\setbeamercolor{item}{fg=AggieCrimson}
\setbeamercolor{block title}{bg=AggieCrimson,fg=white}
\setbeamercolor{block body}{bg=gray!10,fg=black}

\setbeamertemplate{navigation symbols}{}

\title{Efficiency Corrections Comparison}
\subtitle{Current vs. Previous Results}
\author{Chatura Kuruppu}
\institute{Department of Physics \\ New Mexico State University}
\date{\today}

\begin{document}

\begin{frame}
    \titlepage
\end{frame}
"""
    
    # Inject Math Methodology Slides
    tex_content += math_slides

    for target in targets:
        # Create a section frame for each target
        tex_content += f"\\section{{{target} Target}}\n"
        tex_content += r"""
\begin{frame}
    \vfill
    \centering
    \begin{beamercolorbox}[sep=8pt,center,shadow=true,rounded=true]{title}
        \usebeamerfont{title}""" + f"{target} Efficiency Corrections" + r"""\par%
    \end{beamercolorbox}
    \vfill
\end{frame}
"""
        
        for y_type in yield_types:
            for e_type_prefix, e_type_title in eff_types:
                old_file = f"{old_dir}/E_{y_type}_{e_type_prefix}_{target}.pdf"
                new_file = f"{new_dir}/E_{y_type}_{e_type_prefix}_{target}.pdf"
                
                tex_content += r"""
\begin{frame}{""" + f"{target} Target: {y_type.capitalize()} Yield" + r"""}{""" + e_type_title + r"""}
    \begin{figure}
        \centering
        \begin{subfigure}{0.49\textwidth}
            \centering
            % Constrain height so it fits on a 16:9 slide
            \includegraphics[width=\linewidth,height=0.65\textheight,keepaspectratio]{""" + old_file + r"""}
            \caption{Previous}
        \end{subfigure}\hfill
        \begin{subfigure}{0.49\textwidth}
            \centering
            \includegraphics[width=\linewidth,height=0.65\textheight,keepaspectratio]{""" + new_file + r"""}
            \caption{Current}
        \end{subfigure}
    \end{figure}
\end{frame}
"""

        # Signal Efficiency Slide
        old_sig_file = f"{old_dir}/E_final_signal_{target}.pdf"
        new_sig_file = f"{new_dir}/E_final_signal_{target}.pdf"
        
        tex_content += r"""
\begin{frame}{""" + f"{target} Target: Signal Efficiency" + r"""}{Final Corrections}
    \begin{figure}
        \centering
        \begin{subfigure}{0.49\textwidth}
            \centering
            \includegraphics[width=\linewidth,height=0.65\textheight,keepaspectratio]{""" + old_sig_file + r"""}
            \caption{Previous}
        \end{subfigure}\hfill
        \begin{subfigure}{0.49\textwidth}
            \centering
            \includegraphics[width=\linewidth,height=0.65\textheight,keepaspectratio]{""" + new_sig_file + r"""}
            \caption{Current}
        \end{subfigure}
    \end{figure}
\end{frame}
"""

    tex_content += "\\end{document}\n"

    with open(output_slides_file, "w") as f:
        f.write(tex_content)
    print(f"LaTeX slides '{output_slides_file}' generated successfully.")


if __name__ == "__main__":
    # Generate both documents
    generate_report()
    generate_slides()