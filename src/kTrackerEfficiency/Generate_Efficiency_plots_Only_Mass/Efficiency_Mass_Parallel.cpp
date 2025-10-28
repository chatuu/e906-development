/**
 * @file Efficiency_Mass_Parallel.cpp
 * @brief Calculates D2 efficiency for a single, specified mass bin.
 * @details This program is intended to be called by a parallel processing script.
 * It takes one command-line argument: the mass bin index (0-10).
 * It reads from "messy" and "clean" ROOT files, applies cuts for
 * the specified mass bin, and generates efficiency plots and a ROOT file.
 */

#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCut.h>      // <--- THIS IS THE FIX
#include <iostream>
#include <stdexcept> // For exception handling
#include "chuckcuts.h" // Assumed to be in the include path

using namespace std;

// ---------------- Input files ----------------
const TString messyFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_messy_occ_pTxFweight_v2.root";
const TString cleanFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_clean_occ_pTxFweight_v2.root";

// ---------------- Binning ----------------
const TString binVar = "D2";
const double binLow = 0.0;
const double binHigh = 400.0;
const int binCount = 16;

// ---------------- Mass bins ----------------
const int NBINS = 11;
const Double_t edges[NBINS + 1] = {4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7};

// ---------------- Forward declarations ----------------
void plotAbsolute(TCut CcutsTemp, TCut McutsTemp, TString outTag, TTree* cleanTree, TTree* messyTree);
void getErrors(TCut cutsTemp, TString outTag, TH1D* hClean, TH1D* hMessy, TTree* cleanTree);

/**
 * @brief Generates a string label for a given mass bin.
 * @param iMass The index of the mass bin.
 * @return A TString label (e.g., "M4.2to4.5").
 */
TString get_mass_label(int iMass) {
    return Form("%.2f #get Mass #lt %.2f", edges[iMass], edges[iMass+1]);
}

/**
 * @brief Generates a TCut object for a given mass bin.
 * @param iMass The index of the mass bin.
 * @return A TCut object defining the mass range.
 */
TCut get_mass_cut(int iMass) {
    return Form("mass >= %f && mass < %f", edges[iMass], edges[iMass+1]);
}

/**
 * @brief Main function: runs the analysis for a single mass bin.
 */
int main(int argc, char* argv[]) {
    // --- Argument Parsing ---
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <mass_bin_index>" << endl;
        cerr << "Example: " << argv[0] << " 0" << endl;
        return 1;
    }

    int imass = -1;
    try {
        imass = std::stoi(argv[1]);
    } catch (const std::exception& e) {
        cerr << "Error: Invalid mass bin index '" << argv[1] << "'. Must be an integer." << endl;
        return 1;
    }

    if (imass < 0 || imass >= NBINS) {
        cerr << "Error: Mass bin index " << imass << " is out of range. Must be 0-" << (NBINS - 1) << "." << endl;
        return 1;
    }

    cout << "Processing mass bin " << imass << "..." << endl;

    // --- Base cuts ---
    TCut baseCcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && DYCut_MC;

    TCut baseMcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp &&
                    occCuts_2111v42 && DYCut_MC;

    // --- Open files ---
    TFile* cleanFile = TFile::Open(cleanFileName, "READ");
    if (!cleanFile || cleanFile->IsZombie()) {
        cerr << "Error opening clean file: " << cleanFileName << endl;
        return 1;
    }
    TTree* cleanTree = (TTree*)cleanFile->Get("Tree");
    if (!cleanTree) {
        cerr << "Error getting 'Tree' from clean file." << endl;
        cleanFile->Close();
        delete cleanFile;
        return 1;
    }

    TFile* messyFile = TFile::Open(messyFileName, "READ");
    if (!messyFile || messyFile->IsZombie()) {
        cerr << "Error opening messy file: " << messyFileName << endl;
        cleanFile->Close();
        delete cleanFile;
        return 1;
    }
    TTree* messyTree = (TTree*)messyFile->Get("Tree");
    if (!messyTree) {
        cerr << "Error getting 'Tree' from messy file." << endl;
        messyFile->Close();
        delete messyFile;
        cleanFile->Close();
        delete cleanFile;
        return 1;
    }

    // --- Process the single mass bin ---
    TCut masscut = get_mass_cut(imass);
    TString masslabel = get_mass_label(imass);

    TCut CcutsTemp = baseCcut && masscut;
    TCut McutsTemp = baseMcut && masscut;

    plotAbsolute(CcutsTemp, McutsTemp, masslabel, cleanTree, messyTree);

    // --- Cleanup ---
    messyFile->Close();
    delete messyFile;
    
    cleanFile->Close();
    delete cleanFile;

    cout << "Finished processing mass bin " << imass << "." << endl;
    return 0;
}

/**
 * @brief Fills histograms and calls the error calculation.
 */
void plotAbsolute(TCut CcutsTemp, TCut McutsTemp, TString outTag, TTree* cleanTree, TTree* messyTree) {
    TString hCname = "hClean_" + outTag;
    TString hMname = "hMessy_" + outTag;
    TH1D *hClean = new TH1D(hCname, hCname, binCount, binLow, binHigh);
    TH1D *hMessy = new TH1D(hMname, hMname, binCount, binLow, binHigh);

    cleanTree->Draw(binVar + ">>" + hClean->GetName(), "ReWeight*" + TString(CcutsTemp));
    messyTree->Draw(binVar + ">>" + hMname, "ReWeight*" + TString(McutsTemp));
    
    getErrors(CcutsTemp, outTag, hClean, hMessy, cleanTree);

    delete hClean;
    delete hMessy;
}

/**
 * @brief Calculates efficiencies, errors, and generates plots/files.
 */
void getErrors(TCut cutsTemp, TString outTag, TH1D* hClean, TH1D* hMessy, TTree* cleanTree) {
    double weightSqClean = 0;
    double weightSumSqClean = 0;

    Double_t lowIntensityEdge, highIntensityEdge;
    TString intensityBin;

    TH1D* hWeightTemp = new TH1D("hWeightTemp", "hWeightTemp", 100, 0, 1e10);

    TGraphAsymmErrors *gr = new TGraphAsymmErrors(binCount);
    gr->SetMinimum(0);
    gr->SetMaximum(1.2);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->SetTitle(" ");

    for (int i = 1; i <= binCount; i++) {
        lowIntensityEdge = hClean->GetBinLowEdge(i);
        highIntensityEdge = lowIntensityEdge + hClean->GetBinWidth(i);

        intensityBin = TString::Format("D2 > %f && D2 < %f", lowIntensityEdge, highIntensityEdge);

        cleanTree->Draw("ReWeight>>hWeightTemp", cutsTemp && intensityBin, "goff");
        weightSumSqClean = hWeightTemp->GetMean() * hWeightTemp->GetEntries();
        weightSumSqClean = weightSumSqClean * weightSumSqClean;

        cleanTree->Draw("ReWeight*ReWeight>>hWeightTemp", cutsTemp && intensityBin, "goff");
        weightSqClean = hWeightTemp->GetMean() * hWeightTemp->GetEntries();

        double x = hMessy->GetBinCenter(i);
        double y = (hClean->GetBinContent(i) != 0) ? hMessy->GetBinContent(i) / hClean->GetBinContent(i) : 0.0;

        double w = (weightSumSqClean != 0) ? weightSqClean / weightSumSqClean : 0.0;

        double midWilson = (y + 0.5 * w) / (1 + w);
        double eyHigh = 0, eyLow = 0;
        if (!(y == 0 && w == 0)) {
            eyHigh = midWilson + TMath::Sqrt(y*(1-y)*w + 0.25*w*w)/(1+w) - y;
            eyLow  = y - (midWilson - TMath::Sqrt(y*(1-y)*w + 0.25*w*w)/(1+w));
        }

        gr->SetPoint(i - 1, x, y);
        gr->SetPointError(i - 1, 0, 0, eyLow, eyHigh);
    }

    // --- Draw results ---
    TCanvas* c1 = new TCanvas();
    TString title = "D2_Efficiency_" + outTag;
    gr->SetTitle(title + ";D2;Efficiency");
    gr->Draw("AP");

    c1->SaveAs("D2_occ/"+title+".eps");

    // --- Save to file ---
    TFile* outFile = new TFile("D2_occ/"+title+".root","recreate");
    if (outFile && !outFile->IsZombie()) {
        gr->Write("eff_original"); 
        outFile->Close();
    } else {
        cerr << "Error creating output file: " << "D2_occ/"+title+".root" << endl;
    }

    // --- Cleanup local objects ---
    delete hWeightTemp; 
    delete outFile;     
    delete c1;          
    delete gr;          
}