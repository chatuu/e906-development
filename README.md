# E906 Development - High Energy Physics Data Analysis Utilities

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

## Project Structure (Example)

While the exact structure can evolve, a typical setup implied by the scripts might be:

e906-development/
│
├── scripts/                      # Main analysis scripts
│   └── process_data.py           # Example name for the main script we've developed
│
├── config/                       # Configuration files (if any, or embedded in scripts)
│
├── res/                          # Resource files
│   ├── GoodSpills/
│   │   └── GoodSpills.txt        # Text file containing list of good spillIDs
│   └── ROOTFiles/                # (Likely external or linked) Directory for input ROOT files
│       ├── Hugo/
│       │   └── ... (data and MC .root files)
│       └── MixedEvents/
│           └── ... (.root files for mixed and empty flask data)
│
├── output/                       # Output files like CSV tables
│   └── cut_table_ordered_conditional_goodspills.csv
│
├── notebooks/                    # Jupyter notebooks for exploratory analysis (optional)
│
├── .gitignore
└── README.md

*(Please adjust the above structure to reflect your actual repository layout.)*

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