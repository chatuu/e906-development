import os

def main():
    tex_filename = "hodoscope_slides.tex"
    planes = ["H1B", "H1T", "H2B", "H2T", "H3B", "H3T", "H4B", "H4T"]
    
    # Data for the Capped Error Bar Table
    capped_data = [
        ("H1B", 5, "0.1592", "0.0981"),
        ("H1B", 19, "0.0707", "0.0690"),
        ("H1T", 5, "0.0704", "0.0479"),
        ("H1T", 7, "0.0378", "0.0378"),
        ("H2B", 2, "0.5164", "0.4000"),
        ("H2B", 8, "0.0006", "0.0005"),
        ("H2T", 8, "0.0004", "0.0004"),
        ("H2T", 9, "0.0004", "0.0004"),
        ("H3B", 4, "0.0086", "0.0075"),
        ("H3B", 16, "0.0318", "0.0182"),
        ("H3T", 8, "0.0011", "0.0004"),
    ]

    # Start LaTeX document
    latex_code = r"""\documentclass[aspectratio=169]{beamer}
\usepackage{graphicx}
\usepackage{geometry}
\usepackage{booktabs}
\usepackage{amsmath}

% Theme settings
\usetheme{Madrid}
\usecolortheme{default}

\title{Hodoscope Paddle Efficiencies}
\subtitle{Comparison: RS67 vs RS57-70 Final}
\author{Chatura Kuruppu}
\date{\today}

\begin{document}

% Title Slide
\begin{frame}
    \titlepage
\end{frame}

% Slide 01: Input Files
\begin{frame}
    \frametitle{Input Data}
    \textbf{Input files used:}
    \begin{itemize}
        \item \texttt{effies\_rs57.xlsx}
        \item \texttt{effies\_rs59.xlsx}
        \item \texttt{effies\_rs62.xlsx}
        \item \texttt{effies\_rs67.xlsx}
        \item \texttt{effies\_rs70.xlsx}
    \end{itemize}
    
    \vspace{0.5cm}
    \begin{block}{Note}
    The current study is only focused on runs 2 and 3 only. We will do the same study for runs 5 and 6 later as a separate study.
    \end{block}
\end{frame}

% Slide 02: Procedure
\begin{frame}
    \frametitle{Analysis Procedure}
    \begin{itemize}
        \item \textbf{Average Efficiency:} Calculated the arithmetic mean of efficiency across the runs.
        \item \textbf{Propagated Error:} Errors were combined in quadrature:
        \[ \sigma_{prop} = \frac{1}{N} \sqrt{\sum_{i=1}^{N} \sigma_i^2} \]
        \item \textbf{Standard Deviation ($\sigma_{sample}$):} Calculated the sample standard deviation of the efficiency distribution to determine the spread:
        \[ \sigma_{sample} = \sqrt{\frac{1}{N-1} \sum_{i=1}^{N} (\epsilon_i - \bar{\epsilon})^2} \]
        \item \textbf{Error Bar Substitution:} If the standard deviation is larger than the maximum propagated error, both lower and upper error bars are replaced by it:
        \[ \text{If } \sigma_{sample} > \max(\sigma_{prop}^{up}, \sigma_{prop}^{down}) \implies \sigma^{up,down} = \sigma_{sample} \]
        \item \textbf{Physical Limit Capping:} To respect physical boundaries, if the average efficiency plus the upper error bar exceeds 1.0 (100\%), the upper error bar is capped:
        \[ \text{If } \bar{\epsilon} + \sigma^{up} > 1.0 \implies \sigma^{up} = 1.0 - \bar{\epsilon} \]
    \end{itemize}
\end{frame}

% Slide 03: Capped Error Table
\begin{frame}
    \frametitle{Paddles with Upper Error Bar Capped at 1.0 Bound}
    \begin{table}
        \centering
        \begin{tabular}{llcc}
            \toprule
            \textbf{Plane} & \textbf{Paddle ID} & \textbf{Original $\sigma_{up}$} & \textbf{Capped $\sigma_{up}$} \\
            \midrule
"""
    
    # Inject Table Data
    for row in capped_data:
        latex_code += f"            {row[0]} & {row[1]} & {row[2]} & {row[3]} \\\\\n"
        
    latex_code += r"""            \bottomrule
        \end{tabular}
    \end{table}
\end{frame}
"""

    # Slides 04+: Generate a slide for each plane's hodoscope efficiency plot
    for plane in planes:
        plot_file = f"eff_plot_{plane}_root.pdf"
        latex_code += f"""
% Efficiency Plot Slide for {plane}
\\begin{{frame}}
    \\frametitle{{Hodoscope Efficiency: {plane}}}
    \\begin{{center}}
        \\IfFileExists{{{plot_file}}}{{
            \\includegraphics[width=0.9\\textwidth,height=0.8\\textheight,keepaspectratio]{{{plot_file}}}
        }}{{
            \\textbf{{Plot missing: {plot_file}}}
        }}
    \\end{{center}}
\\end{{frame}}
"""

    # Appendix: Efficiency Distributions
    latex_code += r"""
% ---------------------------------------------------------
% APPENDIX
% ---------------------------------------------------------
\appendix

\begin{frame}
    \centering
    \Huge \textbf{Appendix} \\
    \vspace{0.5cm}
    \Large Hodoscope Efficiency Distributions with Mean and RMS
\end{frame}
"""

    # Generate appendix slides looping over planes and paddle counts
    for plane in planes:
        # H1 has 23 paddles, H2-H4 have 16 paddles
        max_paddle = 23 if "H1" in plane else 16
        
        for paddle in range(1, max_paddle + 1):
            dist_plot = f"paddle_distributions/eff_dist_{plane}_paddle_{paddle}.pdf"
            latex_code += f"""
\\begin{{frame}}
    \\frametitle{{Efficiency Distribution: {plane} - Paddle {paddle}}}
    \\begin{{center}}
        \\IfFileExists{{{dist_plot}}}{{
            \\includegraphics[width=0.8\\textwidth,height=0.8\\textheight,keepaspectratio]{{{dist_plot}}}
        }}{{
            \\textbf{{Plot missing: {dist_plot}}}
        }}
    \\end{{center}}
\\end{{frame}}
"""

    # Close document
    latex_code += r"\end{document}"

    # Write to file
    with open(tex_filename, "w") as f:
        f.write(latex_code)
        
    print(f"Success! Generated LaTeX file: {tex_filename}")

if __name__ == "__main__":
    main()