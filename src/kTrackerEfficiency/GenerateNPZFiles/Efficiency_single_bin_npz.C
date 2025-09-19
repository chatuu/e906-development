/**
 * @file Efficiency_single_bin_npz.C
 * @brief A ROOT script to calculate D2 efficiency and save it as an .npz file.
 * @author Gemini
 * @date 18 September 2025
 *
 * This script is designed to be executed for a single (xF, mass) bin. It calculates
 * the efficiency and its errors, then saves the resulting data arrays (x, y, y_error_low, y_error_high)
 * into a compressed NumPy (.npz) file.
 *
 * @note This script requires Python 3, NumPy, and Pandas to be installed in the environment.
 */

#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream> // Required for file output
#include <iostream>
#include <vector> // Required for std::vector
#include "chuckcuts.h" // Assumes this header file is in the correct path

using namespace std;

// --- GLOBAL VARIABLE DEFINITIONS ---

// Input files
/// @brief Path to the input 'messy' ROOT file containing reconstructed events.
TString messyFileName = "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_messy_occ_pTxFweight_v2.root";
/// @brief Path to the input 'clean' ROOT file containing generated events.
TString cleanFileName = "/root/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_clean_occ_pTxFweight_v2.root";

// Binning
/// @brief The variable used for efficiency calculation.
TString binVar = "D2";
/// @brief The lower edge of the histogram range.
double binLow = 0.0;
/// @brief The upper edge of the histogram range.
double binHigh = 400.0;
/// @brief The number of bins to use for the efficiency curve.
int binCount = 16;

// TTree pointers
/// @brief Pointer to the TTree in the 'clean' file.
TTree *cleanTree = nullptr;
/// @brief Pointer to the TTree in the 'messy' file.
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


// --- NPZ SAVING FUNCTION ---

/**
 * @brief Saves the data from a TGraphAsymmErrors object into a compressed NumPy .npz file.
 * @param[in] g The TGraphAsymmErrors object containing the data.
 * @param[in] outName The base name for the output file (without extension).
 */
void saveNPZ(TGraphAsymmErrors* g, TString outName) {
    int n = g->GetN();
    vector<double> vx, vy, veyl, veyh;

    for (int i=0; i<n; i++) {
        double x,y;
        g->GetPoint(i,x,y);
        vx.push_back(x);
        vy.push_back(y);
        veyl.push_back(g->GetErrorYlow(i));
        veyh.push_back(g->GetErrorYhigh(i));
    }

    // Write a temporary CSV file to be read by Python
    TString csvName = outName + ".csv";
    ofstream fout(csvName.Data());
    fout << "x,y,y_error_low,y_error_high\n";
    for (size_t i=0; i<vx.size(); i++) {
        fout << vx[i] << "," << vy[i] << "," << veyl[i] << "," << veyh[i] << "\n";
    }
    fout.close();

    // Call an inline Python script to convert the CSV to NPZ
    TString cmd = Form("python3 -c \"import numpy as np; import pandas as pd; "
                       "df = pd.read_csv('%s'); "
                       "np.savez('%s.npz', x=df['x'].to_numpy(), y=df['y'].to_numpy(), "
                       "y_error_low=df['y_error_low'].to_numpy(), y_error_high=df['y_error_high'].to_numpy())\"",
                       csvName.Data(), outName.Data());
    gSystem->Exec(cmd);

    // Clean up the temporary CSV file
    gSystem->Exec(Form("rm %s", csvName.Data()));
}


// --- WORKER FUNCTION ---

/**
 * @brief Processes a single (xF, mass) bin to calculate efficiency and save to .npz.
 * @param[in] ixf The index of the xF bin to process.
 * @param[in] imass The index of the mass bin to process.
 */
void process_and_save_bin(int ixf, int imass) {
    // 1. Get bin-specific information
    TCut xfcut = get_xf_cut(ixf);
    if (TString(xfcut.GetTitle()) == "") return;

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

    // 7. Save the efficiency data to an NPZ file
    gSystem->Exec("mkdir -p D2_npz");
    TString saveName = TString::Format("D2_npz/D2_Efficiency_xF%d_mass%d", ixf, imass);
    saveNPZ(gr, saveName);

    // 8. Clean up memory
    delete hClean;
    delete hMessy;
    delete hWeightTemp;
    delete gr;
}

// --- MAIN DRIVER ---

/**
 * @brief Main driver function and entry point for the script.
 * @param[in] ixf The index of the xF bin to process.
 * @param[in] imass The index of the mass bin to process.
 */
void Efficiency_single_bin_npz(int ixf, int imass) {
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
    process_and_save_bin(ixf, imass);

    // Clean up TFile objects
    delete cleanFile;
    delete messyFile;
}
