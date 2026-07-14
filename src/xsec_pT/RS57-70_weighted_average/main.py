"""
main.py
Execution script orchestrating Drell-Yan Analysis per roadset.
"""

import argparse
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn
import config
from analyzer import DYCrossSectionAnalyzer

# --- File Paths mapped by Roadset ---
FILE_MAP = {
    "RS57": {
        "lh2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LH2_recoeff_hodoeff.root",
        "ld2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LD2_recoeff_hodoeff.root",
        "flask": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_Flask_recoeff_hodoeff.root"
    },
    "RS59": {
        "lh2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LH2_recoeff_hodoeff.root",
        "ld2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LD2_recoeff_hodoeff.root",
        "flask": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_Flask_recoeff_hodoeff.root"
    },
    "RS62": {
        "lh2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/trimmed_RS62_LH2_recoeff_hodoeff.root",
        "ld2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/trimmed_RS62_LD2_recoeff_hodoeff.root",
        "flask": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/trimmed_RS62_Flask_recoeff_hodoeff.root"
    },
    "RS67": {
        "lh2": "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LH2_recoeff_hodoeff.root",
        "ld2": "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LD2_recoeff_hodoeff.root",
        "flask": "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_Flask_recoeff_hodoeff.root"
    },
    "RS70": {
        "lh2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LH2_recoeff_hodoeff.root",
        "ld2": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LD2_recoeff_hodoeff.root",
        "flask": "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_Flask_recoeff_hodoeff.root"
    }
}

def main():
    # Setup Argument Parser
    parser = argparse.ArgumentParser(description="Run DY Analysis for a specific Roadset.")
    parser.add_argument("--roadset", type=str, required=True, choices=FILE_MAP.keys(),
                        help="Specify the roadset to analyze (e.g., RS57)")
    args = parser.parse_args()
    rs = args.roadset

    # OVERRIDE config constants for the selected roadset
    config.set_roadset(rs)

    console = Console()
    
    # Extract specific file lists for the requested roadset
    lh2_files = [FILE_MAP[rs]["lh2"]]
    ld2_files = [FILE_MAP[rs]["ld2"]]
    flask_files = [FILE_MAP[rs]["flask"]]

    out_root = f"XSec_{rs}_Objects.root"

    # 1. Print Config
    config.print_physics_constants()
    console.print(f"\n[bold blue]Initializing DY Cross-Section Analyzer for {rs}...[/bold blue]")

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

        # Stage 1: Instantiate OOP Analyzer with specific files
        analyzer = DYCrossSectionAnalyzer(
            lh2_files=lh2_files,
            ld2_files=ld2_files,
            flask_files=flask_files,
            out_filename=out_root
        )

        # Stage 2: Kinematics
        analyzer.process_kinematics()
        progress.update(task_kin, completed=100)

        # Stage 3: Subtractions
        progress.update(task_sub, completed=100)

        # Stage 4: Cross Sections & Ratios
        analyzer.calculate_cross_sections()
        progress.update(task_xsec, completed=100)

        # Stage 5: LaTeX Generation
        analyzer.generate_latex_appendix()
        progress.update(task_latex, completed=100)

    # 3. Finalize
    analyzer.finalize()
    console.print(f"\n[bold green]✔ All data for {rs} generated successfully. Saved to {out_root}.[/bold green]")

if __name__ == "__main__":
    main()