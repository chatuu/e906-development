import pandas as pd
import numpy as np
import sys

def generate_latex_table(csv_path, output_path):
    """
    Reads efficiency data from a CSV, filters it, and generates a 
    multi-page LaTeX table with continuous vertical and full-width 
    horizontal lines.

    Args:
        csv_path (str): The path to the input CSV file.
        output_path (str): The path for the output .tex file.
    """
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{csv_path}'")
        return

    # --- MODIFICATION: Filter the DataFrame ---
    # Keep only the rows where 'avg_efficiency' is greater than 0 and less than or equal to 1.
    original_row_count = len(df)
    df = df[(df['avg_efficiency'] > 0) & (df['avg_efficiency'] <= 1)].copy()
    filtered_row_count = len(df)
    print(f"Filtered out {original_row_count - filtered_row_count} rows based on avg_efficiency.")
    # --- END MODIFICATION ---

    col_specs = "| l | l | r | r | r | r | r | r |"

    with open(output_path, 'w') as f:
        # --- longtable setup ---
        f.write(f"\\begin{{longtable}}{{{col_specs}}}\n")
        
        # Caption and label for the first page
        f.write("\\caption{Average Efficiency and Errors for Bins in $x_F$ and Mass}\n")
        f.write("\\label{tab:efficiency}\n")
        
        # Header for the first page
        f.write("\\hline\n")
        f.write(
            "        $x_F$ Bin & Mass Bin (GeV/$c^2$) & $N_{\\text{events}}$ & "
            "$<\\epsilon>$ & $\\delta_{\\text{stat}} <\\epsilon>$ & "
            "$\\delta_{\\text{prop}} <\\epsilon>$ & $1/<\\epsilon>$ & "
            "$\\delta(1/<\\epsilon>)$ \\\\\n"
        )
        f.write("\\hline\n")
        f.write("\\endfirsthead\n\n")

        # Header for all subsequent pages
        f.write("\\caption[]{{(Continued)}}\n")
        f.write("\\hline\n")
        f.write(
            "        $x_F$ Bin & Mass Bin (GeV/$c^2$) & $N_{\\text{events}}$ & "
            "$<\\epsilon>$ & $\\delta_{\\text{stat}} <\\epsilon>$ & "
            "$\\delta_{\\text{prop}} <\\epsilon>$ & $1/<\\epsilon>$ & "
            "$\\delta(1/<\\epsilon>)$ \\\\\n"
        )
        f.write("\\hline\n")
        f.write("\\endhead\n\n")

        # Footer for all pages except the last
        f.write("\\hline\n")
        f.write("\\multicolumn{8}{r}{{Continued on next page}} \\\\\n")
        f.write("\\endfoot\n\n")
        
        # Footer for the very last page
        f.write("\\hline\n")
        f.write("\\endlastfoot\n\n")

        # --- Table Body ---
        for _, row in df.iterrows():
            xf_bin = row['xF_bin'].replace('[', '$[$').replace(')', '$)$')
            mass_bin = row['mass_bin'].replace('[', '$[$').replace(')', '$)$')
            n_events = str(int(row['N_events'])) if pd.notna(row['N_events']) else "--"
            
            row_data = [
                xf_bin, mass_bin, n_events,
                f"{row['avg_efficiency']:.4f}" if pd.notna(row['avg_efficiency']) else "--",
                f"{row['stat_error']:.4f}" if pd.notna(row['stat_error']) else "--",
                f"{row['propagated_error']:.4f}" if pd.notna(row['propagated_error']) else "--",
                f"{row['inv_avg_efficiency']:.3f}" if pd.notna(row['inv_avg_efficiency']) else "--",
                f"{row['inv_propagated_error']:.3f}" if pd.notna(row['inv_propagated_error']) else "--",
            ]
            f.write("        " + " & ".join(row_data) + " \\\\\n")

        # --- End longtable ---
        f.write("\\end{longtable}\n")

    print(f"Successfully generated LaTeX longtable at '{output_path}'")

# --- Main execution block ---
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_latex_table.py <input_csv_path> <output_tex_path>")
        sys.exit(1)
        
    input_csv = sys.argv[1]
    output_tex = sys.argv[2]
    generate_latex_table(input_csv, output_tex)
