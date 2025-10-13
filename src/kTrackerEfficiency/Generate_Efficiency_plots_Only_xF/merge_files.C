#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <iostream>

/**
 * @brief Merges individual TGraph ROOT files into a single file.
 * @param[in] outputFileName The name of the final, combined ROOT file.
 */
void merge_files(const char* outputFileName = "AllGraphs_D2_Efficiency.root") {
    // Define the same xF binning as in the processing script to reconstruct filenames
    const int nXfBins = 17;
    const double xfBinEdges[nXfBins + 1] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85
    };

    // Create the final output file
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }
    std::cout << "Created output file: " << outputFileName << std::endl;

    // Loop over all the bins to find their corresponding files
    for (int ibin = 0; ibin < nXfBins; ++ibin) {
        double xfLow = xfBinEdges[ibin];
        double xfHigh = xfBinEdges[ibin + 1];

        // Reconstruct the individual file name created by the parallel script
        TString inputFileName = TString::Format("D2_occ/D2_Efficiency_xF_%.2f_to_%.2f.root", xfLow, xfHigh);

        TFile* inputFile = TFile::Open(inputFileName, "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Warning: Could not open input file " << inputFileName << ". Skipping." << std::endl;
            continue;
        }

        // Get the graph from the individual file
        TGraphAsymmErrors* graph = (TGraphAsymmErrors*)inputFile->Get("eff_D2_vs_xF");

        if (!graph) {
            std::cerr << "Warning: Could not find graph 'eff_D2_vs_xF' in " << inputFileName << ". Skipping." << std::endl;
            inputFile->Close();
            continue;
        }

        // Set a new, unique name for the graph in the combined file
        TString newGraphName = TString::Format("eff_xF_bin_%d", ibin);
        graph->SetName(newGraphName);

        // Change directory to the output file and write the graph
        outputFile->cd();
        graph->Write();

        std::cout << "  -> Merged graph from " << inputFileName << " as " << newGraphName << std::endl;

        // Clean up
        inputFile->Close();
        delete inputFile;
    }

    // Close the final output file
    outputFile->Close();
    delete outputFile;
    std::cout << "\nMerging complete." << std::endl;
}
