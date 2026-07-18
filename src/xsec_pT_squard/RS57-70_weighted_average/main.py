"""
main.py
Execution script orchestrating Drell-Yan Analysis per roadset sequentially.
"""

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
    console = Console()
    
    # Iterate over all roadsets automatically
    for rs, paths in FILE_MAP.items():
        # OVERRIDE config constants for the selected roadset
        config.set_roadset(rs)

        lh2_files = [paths["lh2"]]
        ld2_files = [paths["ld2"]]
        flask_files = [paths["flask"]]

        out_root = f"XSec_{rs}_Objects.root"

        console.print(f"\n[bold blue]========================================================[/bold blue]")
        console.print(f"[bold blue]Initializing DY Cross-Section Analyzer for {rs}...[/bold blue]")
        config.print_physics_constants()

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            
            task_kin = progress.add_task(f"[cyan]({rs}) Extracting Kinematics...", total=100)
            task_sub = progress.add_task(f"[magenta]({rs}) Subtracting Yields...", total=100)
            task_xsec = progress.add_task(f"[green]({rs}) Calculating Absolute XSec...", total=100)
            task_latex = progress.add_task(f"[yellow]({rs}) Generating LaTeX Appendix...", total=100)

            analyzer = DYCrossSectionAnalyzer(
                lh2_files=lh2_files,
                ld2_files=ld2_files,
                flask_files=flask_files,
                out_filename=out_root
            )

            analyzer.process_kinematics()
            progress.update(task_kin, completed=100)

            progress.update(task_sub, completed=100)

            analyzer.calculate_cross_sections()
            progress.update(task_xsec, completed=100)

            analyzer.generate_latex_appendix()
            progress.update(task_latex, completed=100)

        analyzer.finalize()
        console.print(f"[bold green]✔ Data for {rs} generated successfully. Saved to {out_root}.[/bold green]")

    console.print("\n[bold green]✔ Execution complete for all Roadsets.[/bold green]")

if __name__ == "__main__":
    main()