import os

# --- Configuration ---
path_rs67 = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS67/Final/"
path_rs57_70 = "/root/github/e906-development/src/CalculateDoubleDifferentialCrossSection/RS70/"

# --- Target Order Definition ---
# Ensuring LH2 (0) -> LD2 (1) -> Flask/Empty (2)
target_order = ["lh2", "ld2", "flask", "empty"]

def get_target_priority(filename):
    f = filename.lower()
    if "lh2" in f: return 0
    if "ld2" in f: return 1
    if "flask" in f or "empty" in f: return 2
    return 3

# --- Slide 02 File List (Sorted) ---
raw_files = [
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs59/merged_RS59_LD2_3_466.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs59/merged_RS59_LH2_1_465.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs57/merged_RS57_LD2_3_1138.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs57/merged_RS57_LH2_1_1138.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs62/merged_RS62_LH2_1_1234.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs62/merged_RS62_LD2_3_1237.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs70/merged_RS70_LH2_1_264.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs70/merged_RS70_LD2_3_266.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs59/merged_RS59_Empty_2_466.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs57/merged_RS57_Empty_2_1138.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs62/merged_RS62_Empty_2_1234.root",
    #"/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs67/merged_RS67_Empty_2_14.root",
    "/seaquest/users/harshaka/e906_project/e906-root-ana/work_gpvm/scripts/results_runFinalTree/merged_results_e906_mixing/rs70/merged_RS70_Empty_2_267.root",
    "/seaquest/users/apun/e906_projects/rs67_merged_files/merged_RS67_3089LD2.root",
    "/seaquest/users/apun/e906_projects/rs67_merged_files/merged_RS67_3089LH2.root",
    "/seaquest/users/apun/e906_projects/rs67_merged_files/merged_RS67_3089flask.root"
]

rs_order = ["57", "59", "62", "67", "70"]

def file_sort_key(filepath):
    fname = os.path.basename(filepath).lower()
    rs_val = next((i for i, rs in enumerate(rs_order) if f"rs{rs}" in fname), 99)
    tgt_val = get_target_priority(fname)
    return (rs_val, tgt_val)

sorted_filenames = [os.path.basename(f) for f in sorted(raw_files, key=file_sort_key)]

# --- Comprehensive Plot List Construction ---
xf_ranges = [
    "0.00_0.05", "0.05_0.10", "0.10_0.15", "0.15_0.20", "0.20_0.25", 
    "0.25_0.30", "0.30_0.35", "0.35_0.40", "0.40_0.45", "0.45_0.50", 
    "0.50_0.55", "0.55_0.60", "0.60_0.65", "0.65_0.70", "0.70_0.75", "0.75_0.80"
]

raw_plots = [
    "E_final_signal_Flask.pdf", "E_final_signal_LD2.pdf", "E_final_signal_LH2.pdf",
    "E_mix_final_Flask.pdf", "E_mix_final_LD2.pdf", "E_mix_final_LH2.pdf",
    "E_mix_hodo_Flask.pdf", "E_mix_hodo_LD2.pdf", "E_mix_hodo_LH2.pdf",
    "E_mix_reco_Flask.pdf", "E_mix_reco_LD2.pdf", "E_mix_reco_LH2.pdf",
    "E_total_final_Flask.pdf", "E_total_final_LD2.pdf", "E_total_final_LH2.pdf",
    "E_total_hodo_Flask.pdf", "E_total_hodo_LD2.pdf", "E_total_hodo_LH2.pdf",
    "E_total_reco_Flask.pdf", "E_total_reco_LD2.pdf", "E_total_reco_LH2.pdf",
    "Y_corrected_Flask.pdf", "Y_corrected_LD2.pdf", "Y_corrected_LH2.pdf",
    "Y_corrected_Subtracted_LD2.pdf", "Y_corrected_Subtracted_LH2.pdf",
    "Y_mix_Flask.pdf", "Y_mix_LD2.pdf", "Y_mix_LH2.pdf",
    "Y_total_Flask.pdf", "Y_total_LD2.pdf", "Y_total_LH2.pdf"
]

# Append CrossSection plots (both LH2 and LD2, both variants)
for target in ["LH2", "LD2"]:
    for xf in xf_ranges:
        raw_plots.append(f"CrossSection_{target}_xF_{xf}_Centroid_with_logo.pdf")
        raw_plots.append(f"CrossSection_{target}_xF_{xf}_GeoCenter_with_logo.pdf")

def get_plot_sort_tuple(p):
    cat = 999
    if p.startswith("Y_total"): cat = 10
    elif p.startswith("Y_mix"): cat = 11
    elif "reco" in p.lower(): cat = 20
    elif "hodo" in p.lower(): cat = 30
    elif "E_total_final" in p or "E_mix_final" in p: cat = 40
    elif "E_final_signal" in p: cat = 50
    elif "Y_corrected" in p and "Subtracted" not in p: cat = 60
    elif "Y_corrected_Subtracted" in p: cat = 70
    elif "CrossSection_LH2" in p: cat = 80
    elif "CrossSection_LD2" in p: cat = 90
    
    return (cat, get_target_priority(p), p)

ordered_plots = sorted(raw_plots, key=get_plot_sort_tuple)
filtered_plots = [p for p in ordered_plots if ("CrossSection" not in p) or ("with_logo" in p)]

# --- LaTeX Presentation Construction ---
latex = [
    r"\documentclass[aspectratio=169]{beamer}",
    r"\usetheme{Madrid}",
    r"\usecolortheme{whale}",
    r"\usepackage{graphicx}",
    r"\usepackage{booktabs}",
    r"\title{Comparison Between RS67 Vs RS70}",
    r"\author{Chatura Kuruppu}",
    r"\date{\today}",
    r"\begin{document}",
    r"\begin{frame}",
    r"\titlepage",
    r"\end{frame}"
]

# Slide 02: Filenames
latex.append(r"\begin{frame}{Slide 02: Input Files (Roadset \& Target Order)}")
latex.append(r"\begin{columns}")
latex.append(r"\column{0.5\textwidth}\tiny\begin{itemize}")
for i, f in enumerate(sorted_filenames):
    if i == 8: latex.append(r"\end{itemize}\column{0.5\textwidth}\tiny\begin{itemize}")
    latex.append(rf"\item {f.replace('_', r'\_')}")
latex.append(r"\end{itemize}\end{columns}")
latex.append(r"\vspace{0.3cm}\centering\small\textbf{Remark:} Currently we don't process RS5a, RS5b files.")
latex.append(r"\end{frame}")

# Slide 03: Physics Constants Table
latex.append(r"\begin{frame}{Slide 03: Physics Constants \& Normalizations}")
latex.append(r"\centering\footnotesize")
latex.append(r"\begin{tabular}{ll}")
latex.append(r"\toprule \textbf{Parameter} & \textbf{Value} \\ \midrule")
latex.append(r"PoT LH$_2$ & $2.7639 \times 10^{17}$ \\")
latex.append(r"PoT LD$_2$ & $1.3197 \times 10^{17}$ \\")
latex.append(r"PoT Flask & $5.5905 \times 10^{16}$ \\")
latex.append(r"LH$_2$ Target Density & $3.5966$ mol/cm$^2$ \\")
latex.append(r"LD$_2$ Target Density & $8.0431$ mol/cm$^2$ \\")
latex.append(r"LH$_2$ / LD$_2$ Target Length & $50.8$ cm \\")
latex.append(r"Avogadro Constant & $6.022 \times 10^{23}$ \\")
latex.append(r"xF Bin Width & $0.05$ \\")
latex.append(r"Nuclear Int. Length (LH$_2$/LD$_2$) & $52.0$ / $54.7$ g/cm$^2$ \\")
latex.append(r"\bottomrule \end{tabular}")
latex.append(r"\end{frame}")

# Slide 04: Cuts Summary
latex.append(r"\begin{frame}{Slide 04: Event Selection Summary}")
latex.append(r"\small")
latex.append(r"\begin{description}")
latex.append(r"\item[Dimuon Level] $|dx| < 0.25$, mass $\in [4.2, 8.8]$, $xF \in [-0.1, 0.95]$, $xT \in [0.05, 0.58]$, $\chi^2 < 18$.")
latex.append(r"\item[Track Level] $\chi^2_{target} < 15$, $pz_{st1} \in [9, 75]$, Hits $> 13$, Station 1/3 $y$-ratio $< 1$.")
latex.append(r"\item[Occupancy] $D1, D2, D3 < 400$ hits each; Total $D1+D2+D3 < 1000$.")
latex.append(r"\item[Physics] $Z_{vertex} \in [-320, -5]$, Track Separation $< 270$, $|\cos\theta| < 0.5$.")
latex.append(r"\end{description}")
latex.append(r"\end{frame}")

# Comparison Slides Generation Loop
for plot in filtered_plots:
    # 1. Strip the .pdf extension
    clean_title = plot.replace(".pdf", "")
    
    # 2. Reformat CrossSection titles based on exact request
    if clean_title.startswith("CrossSection"):
        parts = clean_title.split('_')
        try:
            target = parts[1] # LH2 or LD2
            xf_min = parts[3] # e.g. 0.10
            xf_max = parts[4] # e.g. 0.15
            variant = parts[5] # Centroid or GeoCenter
            # Format: CrossSection <Target Type: LH2/LD2> <xF min> <=xF < <xF max>
            # (Added the variant at the end so the slide titles are unique)
            display_title = f"CrossSection {target} {xf_min} $\\le$ xF $<$ {xf_max} ({variant})"
        except:
            display_title = clean_title.replace("_", r"\_")
    else:
        # For non-cross-section plots, just replace underscores so LaTeX doesn't crash
        display_title = clean_title.replace("_", r"\_")

    latex.append(rf"\begin{{frame}}{{{display_title}}}")
    latex.append(r"\begin{columns}")
    latex.append(r"\column{0.5\textwidth}\centering \textbf{RS67} \\ \includegraphics[width=\textwidth,height=0.75\textheight,keepaspectratio]{" + path_rs67 + plot + "}")
    latex.append(r"\column{0.5\textwidth}\centering \textbf{RS70} \\ \includegraphics[width=\textwidth,height=0.75\textheight,keepaspectratio]{" + path_rs57_70 + plot + "}")
    latex.append(r"\end{columns}\end{frame}")

latex.append(r"\end{document}")

# Write to file
with open("presentation.tex", "w") as f:
    f.write("\n".join(latex))