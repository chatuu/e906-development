// ---------------- Mass bins (Must match your C++ code) ----------------
const int NBINS = 11;
const Double_t edges[NBINS + 1] = {4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7};

/**
 * @brief Generates a string label for a given mass bin.
 * @param iMass The index of the mass bin.
 * @return A TString label (e.g., "M4.2to4.5").
 */
TString get_mass_label(int iMass) {
    return Form("M%.1fto%.1f", edges[iMass], edges[iMass+1]);
}

/**
 * @brief Main macro function to merge graphs.
 */
void merge_graphs() {
    TString finalFileName = "EfficiencyCurves_All_Mass_bins.root";
    TString outputDir = "D2_occ";

    TFile *finalFile = new TFile(finalFileName, "RECREATE");
    if (!finalFile || finalFile->IsZombie()) {
        cout << "Error: Could not create final file: " << finalFileName << endl;
        return;
    }

    cout << "Creating final file: " << finalFileName << endl;

    for (int imass = 0; imass < NBINS; imass++) {
        TString masslabel = get_mass_label(imass);
        TString inFileName = Form("%s/D2_Efficiency_%s.root", outputDir.Data(), masslabel.Data());
        
        TFile *inFile = TFile::Open(inFileName);
        if (!inFile || inFile->IsZombie()) {
            cout << "Warning: Could not open input file: " << inFileName << endl;
            continue;
        }

        // Get the graph from the intermediate file
        TGraphAsymmErrors *gr = (TGraphAsymmErrors*)inFile->Get("eff_original");
        
        if (!gr) {
            cout << "Warning: 'eff_original' not found in " << inFileName << endl;
            inFile->Close();
            delete inFile;
            continue;
        }

        // Set the new name as requested
        TString newName = Form("mass_bin_%d", imass);
        gr->SetName(newName);
        gr->SetTitle(Form("Efficiency for %.2f #leq Mass < %.2f", edges[imass], edges[imass+1])); // Optional: set a useful title

        // Switch to the final file and write the graph
        finalFile->cd();
        gr->Write();

        cout << "  -> Added graph '" << newName << "' from " << inFileName << endl;

        // Clean up
        inFile->Close();
        delete inFile; 
        // Note: 'gr' is now owned by 'finalFile', so we don't delete it
    }

    // Save and close the final file
    finalFile->Close();
    delete finalFile;

    cout << "Merging complete." << endl;
}