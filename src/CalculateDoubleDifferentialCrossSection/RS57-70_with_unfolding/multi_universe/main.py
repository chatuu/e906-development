"""
main.py
Execution script orchestrating Drell-Yan Analysis.
"""
import sys
import traceback
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn
import config
from analyzer import DYCrossSectionAnalyzer

def main():
    console = Console()
    
    # --- Input Files for RS57 - RS70 ---
    lh2_files = [
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LH2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LH2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/merged_RS62_LH2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LH2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LH2_recoeff_hodoeff_unfolding.root"
    ]

    ld2_files = [
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LD2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LD2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/merged_RS62_LD2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LD2_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LD2_recoeff_hodoeff_unfolding.root"
    ]

    flask_files = [
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_Empty_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/merged_RS62_Empty_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_Flask_recoeff_hodoeff_unfolding.root",
        "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_Empty_recoeff_hodoeff_unfolding.root"
    ]

    mc_messy_files = {
        "LH2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        "LD2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LD2_M027_S001_messy_occ_pTxFweight_v2.root"
    }
    
    mc_clean_files = {
        "LH2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_clean_occ_pTxFweight_v2.root",
        "LD2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LD2_M027_S001_clean_occ_pTxFweight_v2.root"
    }

    config.print_physics_constants()

    try:
        analyzer = DYCrossSectionAnalyzer(
            lh2_files=lh2_files,
            ld2_files=ld2_files,
            flask_files=flask_files,
            mc_messy_files=mc_messy_files,
            mc_clean_files=mc_clean_files,
            out_filename="All_XSec_Objects.root"
        )

        console.print("\n[bold blue]Initializing DY Cross-Section Analyzer for RS57-70...[/bold blue]")

        with Progress(
            SpinnerColumn(),
            TextColumn("{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            
            task_kin = progress.add_task("[cyan]Extracting kinematics...", total=100)
            task_mat = progress.add_task("[magenta]building matrices...", total=100)
            task_unf = progress.add_task("[green]unfolding and bootstrapping...", total=100)

            # 1. Extracting kinematics
            analyzer.process_kinematics()
            progress.update(task_kin, completed=100)

            # 2. Building matrices
            analyzer.build_response_matrix()
            progress.update(task_mat, completed=100)

            # 3. Unfolding and Bootstrapping
            analyzer.calculate_cross_sections()
            analyzer.generate_latex_appendix()
            analyzer.run_unfolding_bootstrap(n_toys=100, fraction=0.7, progress=progress, task_id=task_unf)

        analyzer.finalize()
        console.print("\n[bold green]✔ All histograms, tables, cross-section plots, and overlays[/bold green]")
        console.print("[bold green]generated successfully.[/bold green]")

    except BaseException as e:
        console.print(f"\n[bold red]FATAL CRASH: {type(e).__name__} - {e}[/bold red]")
        sys.exit(1)

if __name__ == "__main__":
    main()