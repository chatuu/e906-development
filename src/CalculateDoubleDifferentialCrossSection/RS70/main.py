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
    
    # 1. Print Config
    config.print_physics_constants()
    console.print("\n[bold blue]Initializing DY Cross-Section Analyzer...[/bold blue]")

    # 2. Instantiate OOP Analyzer
    analyzer = DYCrossSectionAnalyzer(
        lh2_file="../../HodoEfficiency/RS57-70/RS70/merged_RS70_LH2_recoeff_hodoeff.root",
        ld2_file="../../HodoEfficiency/RS57-70/RS70/merged_RS70_LD2_recoeff_hodoeff.root",
        flask_file="../../HodoEfficiency/RS57-70/RS70/merged_RS70_Flask_recoeff_hodoeff.root",
        out_filename="All_XSec_Objects.root"
    )

    # 3. Setup Progress Bar for multi-stage processing
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

        # Stage 1: Kinematics (You can map progress.advance inside the loops of your class)
        analyzer.process_kinematics()
        progress.update(task_kin, completed=100)

        # Stage 2: Subtractions
        # analyzer.generate_subtracted_plots()
        progress.update(task_sub, completed=100)

        # Stage 3: Cross Sections
        analyzer.calculate_cross_sections()
        progress.update(task_xsec, completed=100)

        # Stage 4: LaTeX Generation
        analyzer.generate_latex_appendix()
        progress.update(task_latex, completed=100)

    # 4. Finalize
    analyzer.finalize()
    console.print("\n[bold green]✔ All histograms, tables, cross-section plots, and overlays generated successfully.[/bold green]")

if __name__ == "__main__":
    main()
