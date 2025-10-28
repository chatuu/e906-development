# Mass-Binned D2 Efficiency Curves

This project provides a set of scripts to calculate D2 efficiency (from "messy" and "clean" Monte Carlo files) binned by dimuon mass. The calculation is parallelized by mass bin to significantly speed up processing, and a final script merges all results into a single, convenient ROOT file.

## Requirements

* A C++ compiler (e.g., `g++`)
* [ROOT (CERN)](https://root.cern/) installed and configured (i.e., `root-config` is in your `PATH`).
* The `chuckcuts.h` header file (this is a local dependency, assumed to be in the same directory).

## How to Use

The entire process is automated by the `run_parallel.sh` bash script.

1.  Make the script executable:
    ```bash
    chmod +x run_parallel.sh
    ```

2.  Run the script:
    ```bash
    ./run_parallel.sh
    ```

### What the Script Does

1.  **Compiles:** It first compiles `Efficiency_Mass_Parallel.cpp` into an executable named `runEfficiency` using `g++` and your `root-config` flags.
2.  **Creates Directories:** It creates the output directories `D2_occ/` for plots/files and `D2_occ/logs/` for log files.
3.  **Executes in Parallel:** It loops through all 11 mass bins and launches a separate, backgrounded process for each one (e.g., `./runEfficiency 0 &`, `./runEfficiency 1 &`, ...). This uses a single core for each job.
4.  **Waits:** The script waits for all parallel jobs to complete.
5.  **Merges:** Once all jobs are finished, it executes the `merge_graphs.C` ROOT macro to collect all the results into one file.
6.  **Cleans Up (Optional):** The script includes a commented-out line to remove the intermediate `.root` files from the `D2_occ/` directory after merging.

## File Descriptions

### 1. `run_parallel.sh`

This is the main "entry point" script. It is a bash script responsible for:
* Compiling the C++ program.
* Launching one instance of the compiled program for each mass bin *in parallel*.
* Redirecting the `stdout` and `stderr` of each job to a separate log file in `D2_occ/logs/`.
* Waiting for all parallel jobs to finish.
* Running the `merge_graphs.C` macro to produce the final output file.

### 2. `Efficiency_Mass_Parallel.cpp`

This is the core logic of the analysis, written as a compilable C++ program. It is **not** a ROOT macro.
* It is designed to be run from the command line with **one argument**: the mass bin index (e.g., `0`, `1`, `2`...).
* It opens the "messy" and "clean" ROOT files.
* It applies the base analysis cuts and the specific mass cut corresponding to the index it was given.
* It generates the `TGraphAsymmErrors` (named `eff_original`) for that single bin.
* It saves the graph to an intermediate, bin-specific ROOT file (e.g., `D2_occ/D2_Efficiency_M4.2to4.5.root`) and also saves an `.eps` plot.
* It properly cleans up all ROOT objects from memory to prevent leaks.

### 3. `merge_graphs.C`

This is a simple ROOT macro designed to run *after* all parallel jobs are complete.
* It opens the final output file, `EfficiencyCurves_All_Mass_bins.root`, in `RECREATE` mode.
* It loops from `imass = 0` to `10`.
* In each loop, it opens the corresponding intermediate file (e.g., `D2_occ/D2_Efficiency_M4.2to4.5.root`).
* It retrieves the `TGraphAsymmErrors` object `eff_original`.
* It renames the graph to the standard format: **`mass_bin_<index>`** (e.g., `mass_bin_0`).
* It writes this newly named graph to the final output file.
* Finally, it closes all files.

## Final Output

### `EfficiencyCurves_All_Mass_bins.root`

This is the primary output of the entire process. It is a single ROOT file that contains the 11 `TGraphAsymmErrors` objects, one for each mass bin. The objects are named according to their bin index:

* `mass_bin_0`
* `mass_bin_1`
* `mass_bin_2`
* ...
* `mass_bin_10`

You can open this file in ROOT and easily draw or compare the different efficiency curves.
