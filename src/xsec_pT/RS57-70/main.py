"""
main.py
Execution script orchestrating Drell-Yan Analysis.
"""

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn
import config
from analyzer import DYCrossSectionAnalyzer

def main():
    console = Console()
    
    # --- Input Files for RS57 - RS70 ---
    lh2_files = [
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LH2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LH2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/merged_RS62_LH2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LH2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LH2_recoeff_hodoeff.root"
    ]

    ld2_files = [
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LD2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LD2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/merged_RS62_LD2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LD2_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LD2_recoeff_hodoeff.root"
    ]

    flask_files = [
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_Flask_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/merged_RS62_Flask_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_Flask_recoeff_hodoeff.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_Flask_recoeff_hodoeff.root"
    ]

    # 1. Print Config
    config.print_physics_constants()
    console.print("\n[bold blue]Initializing DY Cross-Section Analyzer for RS57-70...[/bold blue]")

    # 2. Setup Progress Bar for multi-stage processing
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console
    ) as progress:
        
        # Define Tasks
        task_kin = progress.add_task("[cyan]Extracting Kinematics & Plotting Mass...", total=100)
        task_sub = progress.add_task("[magenta]Subtracting Yields & Flask Backgrounds...", total=100)
        task_xsec = progress.add_task("[green]Calculating Absolute Cross-Sections...", total=100)
        task_latex = progress.add_task("[yellow]Generating LaTeX Appendix...", total=100)

        # Stage 1: Instantiate OOP Analyzer with the lists of files
        analyzer = DYCrossSectionAnalyzer(
            lh2_files=lh2_files,
            ld2_files=ld2_files,
            flask_files=flask_files,
            out_filename="All_XSec_Objects.root"
        )

        # Stage 2: Kinematics (Files are read, filtered, and concatenated here)
        analyzer.process_kinematics()
        progress.update(task_kin, completed=100)

        # Stage 3: Subtractions (Implicitly handled inside calculate_cross_sections in updated analyzer)
        progress.update(task_sub, completed=100)

        # Stage 4: Cross Sections & Ratios (Calculates for both geometric & true pT centroids)
        analyzer.calculate_cross_sections()
        progress.update(task_xsec, completed=100)

        # Stage 5: LaTeX Generation
        analyzer.generate_latex_appendix()
        progress.update(task_latex, completed=100)

    # 3. Finalize
    analyzer.finalize()
    console.print("\n[bold green]✔ All histograms, tables, 6 single differential cross-section plots (geom & true pT), ratio plots, and overlays generated successfully.[/bold green]")

if __name__ == "__main__":
    main()