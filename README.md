# E906 Development - High Energy Physics Data Analysis for DY Absolute Cross-Section

This repository contains Python scripts and potentially other utilities for processing, analyzing, and managing data related to High Energy Physics experiments, with a likely focus on the E906/SeaQuest experiment. The primary script facilitates the generation of cut flow tables from experimental and Monte Carlo (MC) data stored in ROOT files.

## Features

* **Data Processing from ROOT files**: Reads particle physics data (e.g., dimuon events) using `uproot`.
* **Handles Multiple Datasets**: Processes experimental data, various MC simulations, mixed event samples, and empty flask data.
* **Sequential Cut Application**: Applies a series of predefined selection cuts to filter events.
* **Weighted MC Events**: Calculates weighted sums for MC events using a 'ReWeight' column.
* **Conditional Cuts**: Supports specialized cuts like the "Good Spills Cut" which is applied conditionally to specific datasets based on an external `spillID` list.
* **Beam Offset Correction**: Applies run-independent beam offset corrections (currently fixed).
* **Parallel Processing**: Utilizes `concurrent.futures.ProcessPoolExecutor` for faster processing of multiple files.
* **Cut Flow Table Generation**:
    * Aggregates results into a pandas DataFrame.
    * Calculates derived columns such as 'Total MC', 'Purity (DY MC)', 'Efficiency (DY MC)'.
    * Prints a formatted cut flow table to the console using `tabulate`.
    * Saves the resulting cut flow table to a CSV file with a defined column and row order.
* **Configuration**: Easily configurable through Python dictionaries for cuts, variable lists, and file paths.

## Project Structure (Subjected to change)

While the exact structure can evolve, a typical setup implied by the scripts might be:

### Key Directories

* **`analysis/`**: Contains the primary analysis scripts (e.g., `process_cutflow.py` for generating cut flow tables) and Jupyter notebooks (`notebooks/`) for exploratory data analysis, visualizations, and specific studies.
* **`config/`**: (Optional) If you externalize configurations (e.g., cut parameters, file lists, run settings), they would reside here, perhaps as YAML or JSON files.
* **`data/`**: Intended for small, version-controlled resource files essential for running the analysis, such as the `good_spills/GoodSpills.txt`.
    * **Note on Large Data**: Large experimental data files (e.g., ROOT files) should generally **not** be stored directly in the Git repository. This directory might contain instructions (e.g., a `data/README.md`) on how to obtain or link to the necessary large datasets, or it might hold very small, representative samples for testing if feasible.
* **`docs/`**: For more detailed documentation that goes beyond the `README.md`, such as theoretical notes, specific setup guides, or in-depth explanations of analysis methodologies.
* **`output/`**: A suggested directory for saving generated files like tables (CSVs) and plots. It's common practice to add this directory or its contents to `.gitignore` if the outputs are large or frequently regenerated.
* **`tests/`**: For automated tests (e.g., using `pytest` or `unittest`) to ensure code correctness and prevent regressions.
* **`utils/`**: For helper Python modules or scripts containing reusable functions, such as common physics calculations, data loading utilities, or plotting helpers.
* **Root Files (`.gitignore`, `LICENSE`, `README.md`)**: Standard top-level repository files.

## Prerequisites

* Python 3.x
* The following Python libraries are required:
    * `uproot`: For reading ROOT files.
    * `pandas`: For data manipulation and analysis.
    * `numpy`: For numerical operations.
    * `tabulate`: For pretty-printing tables to the console.

You can install these dependencies using pip:

```bash
pip install uproot pandas numpy tabulate