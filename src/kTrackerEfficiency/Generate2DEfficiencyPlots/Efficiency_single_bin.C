/**
 * @file Efficiency_single_bin.C
 * @brief A ROOT script to calculate and plot D2 efficiency for a single (xF, mass) bin.
 * @author Gemini
 * @date 18 September 2025
 *
 * This script is designed to be executed by a shell script that iterates through all bins,
 * allowing for parallel processing at the process level.
 * It takes two integer arguments (ixf, imass) to specify which bin to process.
 * The output is a PDF plot and a ROOT file containing the final efficiency graph.
 */

#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include "chuckcuts.h" // Assumes this header file is in the correct path

using namespace std;

// --- GLOBAL VARIABLE DEFINITIONS ---

// Input files
/// @brief Path to the input 'messy' ROOT file containing reconstructed events.
TString messyFileName = "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_messy_occ_pTxFweight_v2.root";
/// @brief Path to the input 'clean' ROOT file containing generated events.
TString cleanFileName = "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_clean_occ_pTxFweight_v2.root";

// Binning
/// @brief The variable to be plotted on the x-axis of the efficiency curve.
TString binVar = "D2";
/// @brief The lower edge of the histogram range for the binned variable.
double binLow = 0.0;
/// @brief The upper edge of the histogram range for the binned variable.
double binHigh = 400.0;
/// @brief The number of bins to use for the efficiency curve.
int binCount = 16;

// TTree pointers
/// @brief Pointer to the TTree in the 'clean' file. Initialized in the main driver.
TTree *cleanTree = nullptr;
/// @brief Pointer to the TTree in the 'messy' file. Initialized in the main driver.
TTree *messyTree = nullptr;

// Mass bins
/// @brief The total number of mass bins.
const int NBINS = 11;
/// @brief Array defining the edges of the mass bins in GeV/c^2.
Double_t edges[NBINS + 1] = {4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7};

// xF bins
/// @brief The total number of xF bins.
const int nXfBins = 17;
/// @brief Array defining the edges of the xF bins.
const double xF_edges[nXfBins + 1] = {
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85
};


// --- HELPER FUNCTIONS ---

/**
 * @brief Generates a string label for a given mass bin index.
 * @param[in] iMass The index of the mass bin.
 * @return A TString object representing the label (e.g., "M4.2to4.5").
 */
TString get_mass_label(int iMass) {
    return Form("M%.1fto%.1f", edges[iMass], edges[iMass+1]);
}

/**
 * @brief Generates a TCut object for a given mass bin index.
 * @param[in] iMass The index of the mass bin.
 * @return A TCut object for selecting events within the mass bin.
 */
TCut get_mass_cut(int iMass) {
    return Form("mass > %f && mass < %f", edges[iMass], edges[iMass+1]);
}

/**
 * @brief Generates a TCut object for a given xF bin index.
 * @param[in] binNumber The index of the xF bin.
 * @return A TCut object for selecting events within the xF bin.
 */
TCut get_xf_cut(int binNumber){
    if(binNumber >= 0 && binNumber < nXfBins) {
        return TString::Format("xF >= %f && xF < %f", xF_edges[binNumber], xF_edges[binNumber+1]).Data();
    }
    return "";
}

/**
 * @brief Generates a string label for a given xF bin index.
 * @param[in] binNumber The index of the xF bin.
 * @return A TString object representing the label (e.g., "0.05to0.10").
 */
TString get_xf_label(int binNumber){
    if(binNumber >= 0 && binNumber < nXfBins) {
        return TString::Format("%.2fto%.2f", xF_edges[binNumber], xF_edges[binNumber+1]);
    }
    return "";
}


// --- WORKER FUNCTION ---

/**
 * @brief Processes a single (xF, mass) bin to calculate and plot D2 efficiency.
 *
 * This is the core worker function. For a given bin, it performs the following steps:
 * 1. Defines the kinematic cuts.
 * 2. Fills histograms from the 'clean' and 'messy' TTrees.
 * 3. Calculates the efficiency as the ratio of the histograms.
 * 4. Calculates weighted Wilson score errors for each point.
 * 5. Creates a TGraphAsymmErrors object with the results.
 * 6. Formats and styles the plot (titles, colors, ticks, etc.).
 * 7. Saves the resulting plot as a PDF and the graph object to a ROOT file.
 *
 * @param[in] ixf The index of the xF bin to process.
 * @param[in] imass The index of the mass bin to process.
 */
void process_and_plot_bin(int ixf, int imass) {
    // 1. Get bin-specific information
    TCut xfcut = get_xf_cut(ixf);
    TString xflabel = get_xf_label(ixf);
    if (xflabel == "") return;

    TCut masscut = get_mass_cut(imass);
    cout << "Processing xF bin " << ixf << ", mass bin " << imass << endl;

    // 2. Define cuts
    TCut baseCcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && DYCut_MC;
    TCut baseMcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp &&
                    occCuts_2111v42 && DYCut_MC;
    TCut CcutsTemp = baseCcut && xfcut && masscut;
    TCut McutsTemp = baseMcut && xfcut && masscut;

    // 3. Create unique, local histograms
    TString hCname = Form("hClean_xf%d_mass%d", ixf, imass);
    TString hMname = Form("hMessy_xf%d_mass%d", ixf, imass);
    TH1D *hClean = new TH1D(hCname, "", binCount, binLow, binHigh);
    TH1D *hMessy = new TH1D(hMname, "", binCount, binLow, binHigh);

    // 4. Fill histograms
    cleanTree->Draw((binVar + ">>" + hClean->GetName()).Data(), ("ReWeight*(" + TString(CcutsTemp) + ")").Data(), "goff");
    messyTree->Draw((binVar + ">>" + hMessy->GetName()).Data(), ("ReWeight*(" + TString(McutsTemp) + ")").Data(), "goff");

    // 5. Calculate efficiency and errors
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(binCount);

    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerStyle(21);
    gr->SetMinimum(0);
    gr->SetMaximum(1.2);

    TString hWeightTempName = Form("hWeightTemp_xf%d_mass%d", ixf, imass);
    TH1D* hWeightTemp = new TH1D(hWeightTempName, "", 100, 0, 1e10);

    for (int i = 1; i <= binCount; i++) {
        double lowIntensityEdge = hClean->GetBinLowEdge(i);
        double highIntensityEdge = lowIntensityEdge + hClean->GetBinWidth(i);
        TString intensityBin = TString::Format("D2 > %f && D2 < %f", lowIntensityEdge, highIntensityEdge);

        cleanTree->Draw(("ReWeight>>" + hWeightTempName).Data(), CcutsTemp && intensityBin, "goff");
        double weightSumSqClean = hWeightTemp->GetMean() * hWeightTemp->GetEntries();
        weightSumSqClean = weightSumSqClean * weightSumSqClean;

        cleanTree->Draw(("ReWeight*ReWeight>>" + hWeightTempName).Data(), CcutsTemp && intensityBin, "goff");
        double weightSqClean = hWeightTemp->GetMean() * hWeightTemp->GetEntries();

        double x = hMessy->GetBinCenter(i);
        double y = (hClean->GetBinContent(i) != 0) ? hMessy->GetBinContent(i) / hClean->GetBinContent(i) : 0.0;
        double w = (weightSumSqClean != 0) ? weightSqClean / weightSumSqClean : 0.0;

        double midWilson = (y + 0.5 * w) / (1 + w);
        double eyHigh = 0, eyLow = 0;
        if (!(y == 0 && w == 0)) {
            eyHigh = midWilson + TMath::Sqrt(y * (1 - y) * w + 0.25 * w * w) / (1 + w) - y;
            eyLow = y - (midWilson - TMath::Sqrt(y * (1 - y) * w + 0.25 * w * w) / (1 + w));
        }

        gr->SetPoint(i - 1, x, y);
        gr->SetPointError(i - 1, 0, 0, eyLow, eyHigh);
    }

    // 7. Draw and save results
    TString canvasName = Form("canvas_xf%d_mass%d", ixf, imass);
    TCanvas* c1 = new TCanvas(canvasName, "", 800, 600);

    c1->SetTickx(1);
    c1->SetTicky(1);

    TString plotTitle = TString::Format("D2 Efficiency %.2f <= xF < %.2f and %.1f <= mass < %.1f GeV/c^{2}",
                                        xF_edges[ixf], xF_edges[ixf+1], edges[imass], edges[imass+1]);
    gr->SetTitle(plotTitle + ";D2;Efficiency");

    gr->GetXaxis()->CenterTitle(true);
    gr->GetYaxis()->CenterTitle(true);

    gr->Draw("ALP");

    gSystem->Exec("mkdir -p D2_occ");
    TString saveName = TString::Format("D2_Efficiency_xF%d_mass%d", ixf, imass);
    c1->SaveAs("D2_occ/" + saveName + ".pdf");

    TFile* outFile = new TFile("D2_occ/" + saveName + ".root", "RECREATE");
    gr->Write("eff_graph");
    outFile->Close();

    // 8. Clean up memory
    delete hClean;
    delete hMessy;
    delete hWeightTemp;
    delete gr;
    delete c1;
    delete outFile;
}

// --- MAIN DRIVER ---

/**
 * @brief Main driver function and entry point for the script.
 *
 * This function is called by the ROOT CINT/Cling interpreter when the script is executed.
 * It takes the bin indices as arguments, opens the necessary ROOT files,
 * and calls the worker function `process_and_plot_bin` to perform the analysis.
 *
 * @param[in] ixf The index of the xF bin to process.
 * @param[in] imass The index of the mass bin to process.
 */
void Efficiency_single_bin(int ixf, int imass) {
    // Open files and get TTree pointers.
    TFile* cleanFile = TFile::Open(cleanFileName);
    if (!cleanFile || cleanFile->IsZombie()) {
        cerr << "Error: Cannot open clean file: " << cleanFileName << endl;
        return;
    }
    cleanTree = (TTree*)cleanFile->Get("Tree");

    TFile* messyFile = TFile::Open(messyFileName);
    if (!messyFile || messyFile->IsZombie()) {
        cerr << "Error: Cannot open messy file: " << messyFileName << endl;
        delete cleanFile;
        return;
    }
    messyTree = (TTree*)messyFile->Get("Tree");

    if (!cleanTree || !messyTree) {
        cerr << "Error: Could not find TTree named 'Tree' in one of the files." << endl;
        delete cleanFile;
        delete messyFile;
        return;
    }

    // Process the single bin specified by the function arguments.
    process_and_plot_bin(ixf, imass);

    // Clean up TFile objects
    delete cleanFile;
    delete messyFile;
}
