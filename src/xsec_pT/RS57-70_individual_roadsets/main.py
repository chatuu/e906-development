"""
main.py
Execution script orchestrating Drell-Yan Analysis per roadset or full Combined analysis.
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

def run_pipeline(rs_name, lh2_files, ld2_files, flask_files, console):
    """Executes the complete DY cross-section extraction pipeline for given inputs."""
    # OVERRIDE config constants for the selected dataset
    config.set_roadset(rs_name)
    out_root = f"XSec_{rs_name}_Objects.root"

    # 1. Print Config
    config.print_physics_constants()
    console.print(f"\n[bold blue]Initializing DY Cross-Section Analyzer for {rs_name}...[/bold blue]")

    # 2. Setup Progress Bar for multi-stage processing
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console
    ) as progress:
        
        task_kin = progress.add_task("[cyan]Extracting Kinematics & Plotting Mass...", total=100)
        task_sub = progress.add_task("[magenta]Subtracting Yields & Flask Backgrounds...", total=100)
        task_xsec = progress.add_task("[green]Calculating Absolute Cross-Sections...", total=100)
        task_latex = progress.add_task("[yellow]Generating LaTeX Appendix...", total=100)

        # Stage 1: Instantiate OOP Analyzer
        analyzer = DYCrossSectionAnalyzer(
            lh2_files=lh2_files,
            ld2_files=ld2_files,
            flask_files=flask_files,
            out_filename=out_root
        )

        # Stage 2: Kinematics
        analyzer.process_kinematics()
        progress.update(task_kin, completed=100)

        # Stage 3: Subtractions (Implicit in cross-section call, updated visually here)
        progress.update(task_sub, completed=100)

        # Stage 4: Cross Sections & Ratios
        analyzer.calculate_cross_sections()
        progress.update(task_xsec, completed=100)

        # Stage 5: LaTeX Generation
        analyzer.generate_latex_appendix()
        progress.update(task_latex, completed=100)

    # 3. Finalize
    analyzer.finalize()
    console.print(f"[bold green]✔ All data for {rs_name} generated successfully. Saved to {out_root}.[/bold green]\n")

def main():
    # Setup Argument Parser
    parser = argparse.ArgumentParser(description="Run DY Analysis for a specific Roadset or combine all.")
    choices = list(FILE_MAP.keys()) + ["All"]
    parser.add_argument("--roadset", type=str, required=True, choices=choices,
                        help="Specify the roadset to analyze (e.g., RS57), or 'All' to run and plot the combined weighted average overlay.")
    args = parser.parse_args()
    
    console = Console()

    if args.roadset == "All":
        roadsets = list(FILE_MAP.keys())
        
        # 1. Process Individual Roadsets
        for rs in roadsets:
            run_pipeline(
                rs_name=rs,
                lh2_files=[FILE_MAP[rs]["lh2"]],
                ld2_files=[FILE_MAP[rs]["ld2"]],
                flask_files=[FILE_MAP[rs]["flask"]],
                console=console
            )

        # 2. Process Combined (Weighted Average statistically calculated via yielding over total sum PoT)
        console.print("[bold magenta]===== Running Combined Analysis (Weighted Average) =====[/bold magenta]")
        all_lh2 = [FILE_MAP[rs]["lh2"] for rs in roadsets]
        all_ld2 = [FILE_MAP[rs]["ld2"] for rs in roadsets]
        all_flask = [FILE_MAP[rs]["flask"] for rs in roadsets]
        
        run_pipeline(
            rs_name="Combined",
            lh2_files=all_lh2,
            ld2_files=all_ld2,
            flask_files=all_flask,
            console=console
        )

        # 3. Generate Overlays
        console.print("[bold cyan]===== Generating Roadset Overlays vs Weighted Average =====[/bold cyan]")
        DYCrossSectionAnalyzer.plot_roadset_overlays(roadsets=roadsets, combined_label="Combined")
        console.print("[bold green]✔ All overlays generated successfully.[/bold green]")

    else:
        # Standard Single Run
        rs = args.roadset
        run_pipeline(
            rs_name=rs,
            lh2_files=[FILE_MAP[rs]["lh2"]],
            ld2_files=[FILE_MAP[rs]["ld2"]],
            flask_files=[FILE_MAP[rs]["flask"]],
            console=console
        )

if __name__ == "__main__":
    main()