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

    config.print_physics_constants()
    console.print("\n[bold blue]Initializing DY Analyzer (pT distributions per xF bin)...[/bold blue]")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console
    ) as progress:
        
        task_kin = progress.add_task("[cyan]Extracting Kinematics (Binned in xF & pT)...", total=100)
        task_xsec = progress.add_task("[green]Calculating Cross-Sections & PDF Plots...", total=100)

        analyzer = DYCrossSectionAnalyzer(
            lh2_files=lh2_files,
            ld2_files=ld2_files,
            flask_files=flask_files,
            out_filename="All_XSec_Objects_Binned_xF.root"
        )

        # Stage 1: Extraction & Histogram Filling
        analyzer.process_kinematics()
        progress.update(task_kin, completed=100)

        # Stage 2: Background subtraction, Math, and ROOT Canvas Plotting
        analyzer.calculate_cross_sections()
        progress.update(task_xsec, completed=100)

    analyzer.finalize()
    console.print("\n[bold green]✔ Complete! All individual cross-sections, ratio plots, and combined overlays have been generated as PDFs.[/bold green]")

if __name__ == "__main__":
    main()