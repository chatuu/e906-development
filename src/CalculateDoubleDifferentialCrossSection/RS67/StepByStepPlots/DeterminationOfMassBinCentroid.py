import numpy as np
import os
import csv

# Define the same bins used in your analysis
xf_bins_np = np.round(np.arange(0.0, 0.85, 0.05), 2)

latex_filename = "Appendix_MassCentroids.tex"

# Start writing the LaTeX file
with open(latex_filename, "w") as tex_file:
    
    # 1. Write the Introductory Text and Equations
    intro_text = r"""\clearpage
\section{Appendix: Data-Driven Mass Bin Centroids}
\label{sec:appendix_mass_centroids}

In order to accurately plot the double-differential cross-sections, the horizontal placement of the data points must reflect the true physical average of the invariant mass within each finite mass bin, rather than the simple geometric center. Because the dimuon mass distribution falls exponentially, the true centroid is typically skewed toward the lower edge of the bin.

To account for differing backgrounds, tracking efficiencies, and target flask contributions, the centroid calculation must be weighted using the fully corrected equations for the yields. 

For the $\text{LH}_2$ target, the corrected yield is given by:
\begin{equation}
    \label{eq:appendix_yield_lh2}
    Y_{\rm corrected}^{\rm LH2}= \frac{Y^{\text{LH2}}_{\text{total}} - Y^{\text{LH2}}_{\text{mixed}}}{\langle\epsilon^{\text{LH2}}_{\text{signal}}\rangle} 
        - \frac{I_{\text{LH2}}}{I_{\text{flask}}} 
        \left( \frac{Y^{\text{flask}}_{\text{total}} - Y^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right)
\end{equation}

For the $\text{LD}_2$ target, the corrected yield accounts for both the empty flask and the hydrogen contribution:
\begin{eqnarray}
    Y_{\rm corrected}^{\rm LD2}&=& 
        \frac{Y^{\text{LD2}}_{\text{total}} - Y^{\text{LD2}}_{\text{mixed}}}{\langle\epsilon^{\text{LD2}}_{\text{signal}}\rangle} 
        - \frac{I_{\text{LD2}}}{I_{\text{flask}}} 
        \left( \frac{Y^{\text{flask}}_{\text{total}} - Y^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right) \nonumber \\
    & & - \frac{T_{HD}}{T_{HH}} \frac{I_{\rm LD2}}{I_{\rm LH2}}\left[ 
        \frac{Y^{\text{LH2}}_{\text{total}} - Y^{\text{LH2}}_{\text{mixed}}}{\langle\epsilon^{\text{LH2}}_{\text{signal}}\rangle} 
        - \frac{I_{\text{LH2}}}{I_{\text{flask}}} 
        \left( \frac{Y^{\text{flask}}_{\text{total}} - Y^{\text{flask}}_{\text{mixed}}}{\langle\epsilon^{\text{flask}}_{\text{signal}}\rangle} \right)
    \right]
    \label{eq:appendix_yield_ld2}
\end{eqnarray}

The corrected total mass for each target ($\sum M_{\rm corrected}$) is calculated by applying identical background subtractions and efficiency corrections to the sum of the invariant masses of the accepted dimuon pairs. The physical mass bin centroid is then simply the ratio of the corrected mass sum to the corrected yield:
\begin{equation}
    \langle M \rangle_{\rm true} = \frac{\sum M_{\rm corrected}}{Y_{\rm corrected}}
\end{equation}

The following subsections document the 1D mass distributions for the signal, mixed background, and flask components used in these calculations. Alongside the distributions are the precise kinematic values, raw counts, and resulting data-driven centroids for each $x_F$ bin.

"""
    tex_file.write(intro_text)

    # 2. Loop through each xF bin to generate plots and tables
    for i_x in range(len(xf_bins_np) - 1):
        x_low = xf_bins_np[i_x]
        x_high = xf_bins_np[i_x+1]
        
        csv_lh2 = f"Table_Kinematics_LH2_xF_{x_low:.2f}_{x_high:.2f}.csv"
        csv_ld2 = f"Table_Kinematics_LD2_xF_{x_low:.2f}_{x_high:.2f}.csv"
        plot_lh2 = f"MassDist_LH2_xF_{x_low:.2f}_{x_high:.2f}.pdf"
        plot_ld2 = f"MassDist_LD2_xF_{x_low:.2f}_{x_high:.2f}.pdf"
        
        # Subsection header
        tex_file.write(f"\n\\clearpage\n\\subsection{{Kinematics and Distributions for $x_F \\in [{x_low:.2f}, {x_high:.2f})$}}\n")
        
        # Add Figures
        figure_block = f"""
\\begin{{figure}}[htbp]
    \\centering
    \\includegraphics[width=0.49\\textwidth]{{{plot_lh2}}}
    \\hfill
    \\includegraphics[width=0.49\\textwidth]{{{plot_ld2}}}
    \\caption{{1D Mass distributions for $\\text{{LH}}_2$ (left) and $\\text{{LD}}_2$ (right) in the range $x_F \\in [{x_low:.2f}, {x_high:.2f})$. Components shown include target total (red), target mix (blue), flask total (green), and flask mix (orange). The values listed above each bin represent the efficiency-weighted corrected mass $\\sum M / \\epsilon$ and the raw counts $(N)$.}}
    \\label{{fig:mass_dist_xf_{i_x}}}
\\end{{figure}}
"""
        tex_file.write(figure_block)

        # Formatting functions
        def fnum(val_str, fmt=".3f"):
            try:
                v = float(val_str)
                return f"{v:{fmt}}"
            except ValueError:
                return "-"
                
        def fint(val_str):
            try:
                v = float(val_str)
                return f"{int(v)}"
            except ValueError:
                return "-"

        # -------------------------------------------------------------
        # Read LH2 CSV and create LH2 LaTeX Table
        # -------------------------------------------------------------
        if os.path.exists(csv_lh2):
            tex_file.write("\\begin{table}[htbp]\n")
            tex_file.write("    \\centering\n")
            tex_file.write("    \\caption{$\\text{LH}_2$ raw counts, corrected masses, yields, and final data-driven centroids for $x_F \\in [" + f"{x_low:.2f}, {x_high:.2f}" + ").}\n")
            tex_file.write(f"    \\label{{tab:kinematics_lh2_xf_{i_x}}}\n")
            tex_file.write("    \\resizebox{\\textwidth}{!}{\n")
            tex_file.write("    \\begin{tabular}{cc|c|cccc|cc}\n")
            tex_file.write("        \\hline\\hline\n")
            tex_file.write("        Mass & Geo. & True $\\langle M \\rangle$ & $N_{\\rm tot}^{\\rm LH2}$ & $N_{\\rm mix}^{\\rm LH2}$ & $N_{\\rm tot}^{\\rm fl}$ & $N_{\\rm mix}^{\\rm fl}$ & $\\sum M_{\\rm corr}^{\\rm LH2}$ & $Y_{\\rm corr}^{\\rm LH2}$ \\\\\n")
            tex_file.write("        Bin (GeV) & $\\langle M \\rangle$ & (LH2) & & & & & (Num.) & (Denom.) \\\\\n")
            tex_file.write("        \\hline\n")
            
            with open(csv_lh2, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    m_bin = row["Mass Bin"]
                    m_geo = fnum(row["Mass Center"])
                    m_lh2 = fnum(row["LH2 Mass Bin Average"])
                    n_lh2_t = fint(row["N_LH2_total"])
                    n_lh2_m = fint(row["N_LH2_mixed"])
                    n_fl_t  = fint(row["N_flask_total"])
                    n_fl_m  = fint(row["N_flask_mixed"])
                    num_lh2 = fnum(row["Corrected Total Mass LH2 (num)"], ".1f")
                    den_lh2 = fnum(row["Corrected Yield LH2 (denom)"], ".1f")
                    
                    line = f"        {m_bin} & {m_geo} & {m_lh2} & {n_lh2_t} & {n_lh2_m} & {n_fl_t} & {n_fl_m} & {num_lh2} & {den_lh2} \\\\\n"
                    tex_file.write(line)
            
            tex_file.write("        \\hline\\hline\n")
            tex_file.write("    \\end{tabular}\n")
            tex_file.write("    }\n")
            tex_file.write("\\end{table}\n\n")
        else:
            tex_file.write(f"\n\\textit{{Table data missing: {csv_lh2} not found.}}\n\n")

        # -------------------------------------------------------------
        # Read LD2 CSV and create LD2 LaTeX Table
        # -------------------------------------------------------------
        if os.path.exists(csv_ld2):
            tex_file.write("\\begin{table}[htbp]\n")
            tex_file.write("    \\centering\n")
            tex_file.write("    \\caption{$\\text{LD}_2$ raw counts, corrected masses, yields, and final data-driven centroids for $x_F \\in [" + f"{x_low:.2f}, {x_high:.2f}" + ").}\n")
            tex_file.write(f"    \\label{{tab:kinematics_ld2_xf_{i_x}}}\n")
            tex_file.write("    \\resizebox{\\textwidth}{!}{\n")
            tex_file.write("    \\begin{tabular}{cc|c|cccccc|cc}\n")
            tex_file.write("        \\hline\\hline\n")
            tex_file.write("        Mass & Geo. & True $\\langle M \\rangle$ & $N_{\\rm tot}^{\\rm LD2}$ & $N_{\\rm mix}^{\\rm LD2}$ & $N_{\\rm tot}^{\\rm LH2}$ & $N_{\\rm mix}^{\\rm LH2}$ & $N_{\\rm tot}^{\\rm fl}$ & $N_{\\rm mix}^{\\rm fl}$ & $\\sum M_{\\rm corr}^{\\rm LD2}$ & $Y_{\\rm corr}^{\\rm LD2}$ \\\\\n")
            tex_file.write("        Bin (GeV) & $\\langle M \\rangle$ & (LD2) & & & & & & & (Num.) & (Denom.) \\\\\n")
            tex_file.write("        \\hline\n")
            
            with open(csv_ld2, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    m_bin = row["Mass Bin"]
                    m_geo = fnum(row["Mass Center"])
                    m_ld2 = fnum(row["LD2 Mass Bin Average"])
                    n_ld2_t = fint(row["N_LD2_total"])
                    n_ld2_m = fint(row["N_LD2_mixed"])
                    n_lh2_t = fint(row["N_LH2_total"])
                    n_lh2_m = fint(row["N_LH2_mixed"])
                    n_fl_t  = fint(row["N_flask_total"])
                    n_fl_m  = fint(row["N_flask_mixed"])
                    num_ld2 = fnum(row["Corrected Total Mass LD2 (num)"], ".1f")
                    den_ld2 = fnum(row["Corrected Yield LD2 (denom)"], ".1f")
                    
                    line = f"        {m_bin} & {m_geo} & {m_ld2} & {n_ld2_t} & {n_ld2_m} & {n_lh2_t} & {n_lh2_m} & {n_fl_t} & {n_fl_m} & {num_ld2} & {den_ld2} \\\\\n"
                    tex_file.write(line)
            
            tex_file.write("        \\hline\\hline\n")
            tex_file.write("    \\end{tabular}\n")
            tex_file.write("    }\n")
            tex_file.write("\\end{table}\n")
        else:
            tex_file.write(f"\n\\textit{{Table data missing: {csv_ld2} not found.}}\n")

print(f"LaTeX appendix successfully generated: {latex_filename}")
print("You can now compile this file or use \\input{Appendix_MassCentroids.tex} in your main document.")