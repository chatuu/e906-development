import ROOT

def skim_trees(input_filepath, output_filepath):
    # Enable multi-threading for faster processing (optional but recommended)
    ROOT.ROOT.EnableImplicitMT()
    
    # Process the first TTree: 'result'
    # By default, Snapshot creates a new file (RECREATE mode)
    df_result = ROOT.RDataFrame("result", input_filepath)
    df_result_filtered = df_result.Filter("runID > 11500")
    df_result_filtered.Snapshot("result", output_filepath)
    print(f"Saved filtered 'result' tree to {output_filepath}")

    # Process the second TTree: 'result_mix'
    df_mix = ROOT.RDataFrame("result_mix", input_filepath)
    df_mix_filtered = df_mix.Filter("runID > 11500")
    
    # To save the second tree to the SAME file, we must set the snapshot mode to "UPDATE"
    opts = ROOT.RDF.RSnapshotOptions()
    opts.fMode = "UPDATE"
    df_mix_filtered.Snapshot("result_mix", output_filepath, "", opts)
    print(f"Saved filtered 'result_mix' tree to {output_filepath}")

if __name__ == "__main__":
    input_file = "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/merged_RS62_LH2_recoeff_hodoeff.root"
    output_file = "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/trimmed_RS62_LH2_recoeff_hodoeff.root" # Replace with your desired output file name
    
    skim_trees(input_file, output_file)
    print("Done!")