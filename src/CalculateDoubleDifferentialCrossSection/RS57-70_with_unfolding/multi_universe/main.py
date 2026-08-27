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

    # --- Monte Carlo Files for RooUnfold Response ---
    mc_messy_files = {
        "LH2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        "LD2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LD2_M027_S001_messy_occ_pTxFweight_v2.root"
    }
    
    mc_clean_files = {
        "LH2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_clean_occ_pTxFweight_v2.root",
        "LD2": "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LD2_M027_S001_clean_occ_pTxFweight_v2.root"
    }

    config.print_physics_constants()
    console.print("\n[bold blue]Initializing DY Cross-Section Analyzer for RS57-70 (with RooUnfold)...[/bold blue]")

    try:
        # Move instantiation OUTSIDE the progress bar to guarantee init errors print cleanly
        analyzer = DYCrossSectionAnalyzer(
            lh2_files=lh2_files,
            ld2_files=ld2_files,
            flask_files=flask_files,
            mc_messy_files=mc_messy_files,
            mc_clean_files=mc_clean_files,
            out_filename="All_XSec_Objects.root",
            console=console
        )

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        ) as progress:
            
            task_kin = progress.add_task("[cyan]Extracting Kinematics & Plotting Mass...", total=100)
            task_rm  = progress.add_task("[magenta]Building 1D Target-Specific Response Matrices (Messy & Clean)...", total=100)
            task_xsec = progress.add_task("[green]Unfolding & Calculating Cross-Sections...", total=100)
            task_latex = progress.add_task("[yellow]Generating LaTeX Appendix...", total=100)
            task_boot = progress.add_task("[blue]Running 100 Unfolding Toys (Systematics)...", total=100)

            # Stage 2: Kinematics 
            analyzer.process_kinematics()
            progress.update(task_kin, completed=100)

            # Stage 3: Build Separate Response Matrices from MC
            analyzer.build_response_matrix()
            progress.update(task_rm, completed=100)

            # Stage 4: Subtractions, Unfolding, & Cross Sections
            analyzer.calculate_cross_sections()
            progress.update(task_xsec, completed=100)

            # Stage 5: LaTeX Generation
            analyzer.generate_latex_appendix()
            progress.update(task_latex, completed=100)
            
            # Stage 6: Bootstrapped Unfolding Systematics
            analyzer.run_unfolding_bootstrap(n_toys=100, fraction=0.7, progress=progress, task_id=task_boot)
            progress.update(task_boot, completed=100)

        # 3. Finalize
        analyzer.finalize()
        console.print("\n[bold green]✔ All histograms, tables, unfolded cross-section plots, and overlays generated successfully.[/bold green]")
        console.print("[bold cyan]✔ Response Matrices (Messy & Clean) separately saved to 'DY_ResponseMatrices.root'[/bold cyan]")
        console.print("[bold magenta]✔ Toy outputs for systematics saved to 'Unfolding_Toys' inside main output ROOT file.[/bold magenta]")

    except BaseException as e:
        console.print(f"\n[bold red]FATAL CRASH:[/bold red] {type(e).__name__} - {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()