/**
 * @file Efficiency_xF.C
 * @brief Processes a SINGLE xF bin to calculate D2 efficiency.
 * @details This script is designed to be called from a parallel-processing
 * bash script. It takes an integer bin index as an argument and processes

 * only that specific bin.
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>

// NOTE: We do NOT enable ImplicitMT here, as the parallelization
// will be handled by the bash script launching multiple processes.
#include "chuckcuts.h"

// --- Function Declarations ---
void getErrorsAndSave(const TH1D* hMessy, const TH1D* hClean, TTree* cleanTree, TCut cutsTemp, const TString& fileTag, const TString& plotTitle);


/**
 * @brief Main function to process a single xF bin.
 * @param[in] ibin The index of the xF bin to process (e.g., 0 to 16).
 */
void Efficiency_xF(int ibin) {
    // --- Configuration ---
    TString messyFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_messy_occ_pTxFweight_v2.root";
    TString cleanFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_clean_occ_pTxFweight_v2.root";
    const TString binVar = "D2";
    const double binLow = 0.0;
    const double binHigh = 400.0;
    const int binCount = 16;
    const int nXfBins = 17;
    const double xfBinEdges[nXfBins + 1] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85
    };

    if (ibin < 0 || ibin >= nXfBins) {
        std::cerr << "Error: Bin index " << ibin << " is out of range." << std::endl;
        return;
    }

    // --- Setup ---
    TCut baseCcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && DYCut_MC;
    TCut baseMcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp &&
                    occCuts_2111v42 && DYCut_MC;

    TFile* cleanFile = TFile::Open(cleanFileName, "READ");
    TFile* messyFile = TFile::Open(messyFileName, "READ");
    if (!cleanFile || !messyFile || cleanFile->IsZombie() || messyFile->IsZombie()) {
        std::cerr << "Error: Could not open input files for bin " << ibin << std::endl;
        return;
    }
    TTree* cleanTree = (TTree*)cleanFile->Get("Tree");
    TTree* messyTree = (TTree*)messyFile->Get("Tree");
    if (!cleanTree || !messyTree) {
        std::cerr << "Error: Could not retrieve TTrees for bin " << ibin << std::endl;
        return;
    }
    
    gSystem->mkdir("D2_occ", kTRUE);

    // --- Processing for the single specified bin ---
    double xfLow = xfBinEdges[ibin];
    double xfHigh = xfBinEdges[ibin + 1];

    std::cout << "🚀 Processing bin " << ibin << ": " << xfLow << " <= xF < " << xfHigh << std::endl;

    TCut xfcut = TString::Format("xF >= %f && xF < %f", xfLow, xfHigh).Data();
    TCut CcutsTemp = baseCcut && xfcut;
    TCut McutsTemp = baseMcut && xfcut;

    TString plotTitle = TString::Format("D2_Efficiency_%.2f <= xF < %.2f", xfLow, xfHigh);
    TString fileTag   = TString::Format("D2_Efficiency_xF_%.2f_to_%.2f", xfLow, xfHigh);
    TString hCname = TString::Format("hClean_%d", ibin);
    TString hMname = TString::Format("hMessy_%d", ibin);
    
    TH1D* hClean = new TH1D(hCname, plotTitle + ";D2;Events", binCount, binLow, binHigh);
    TH1D* hMessy = new TH1D(hMname, plotTitle + ";D2;Events", binCount, binLow, binHigh);

    messyTree->Draw(binVar + ">>" + hMname, "ReWeight*(" + TString(McutsTemp) + ")", "goff");
    cleanTree->Draw(binVar + ">>" + hCname, "ReWeight*(" + TString(CcutsTemp) + ")", "goff");
    
    getErrorsAndSave(hMessy, hClean, cleanTree, CcutsTemp, fileTag, plotTitle);
    
    delete hClean;
    delete hMessy;

    cleanFile->Close();
    messyFile->Close();
    delete cleanFile;
    delete messyFile;

    std::cout << "✅ Finished bin " << ibin << "." << std::endl;
}


/**
 * @brief Calculates, plots, and saves the efficiency for a single bin.
 * (This function remains unchanged from the previous version)
 */
void getErrorsAndSave(const TH1D* hMessy, const TH1D* hClean, TTree* cleanTree, TCut cutsTemp, const TString& fileTag, const TString& plotTitle) {
    if (!hMessy || !hClean || !cleanTree) return;
    int nBins = hClean->GetNbinsX();
    auto gr = new TGraphAsymmErrors(nBins);
    gr->SetMinimum(0);
    gr->SetMaximum(1.2);
    gr->SetMarkerColor(kAzure + 7);
    gr->SetLineColor(kAzure + 7);
    gr->SetMarkerStyle(20);
    gr->SetTitle(plotTitle + ";D2;Efficiency");
    TH1D* hWeightTemp = new TH1D("hWeightTemp", "hWeightTemp", 100, 0, 1e10);
    for (int i = 1; i <= nBins; ++i) {
        double lowEdge = hClean->GetBinLowEdge(i);
        double highEdge = lowEdge + hClean->GetBinWidth(i);
        TCut intensityBin = TString::Format("D2 > %f && D2 < %f", lowEdge, highEdge).Data();
        cleanTree->Draw("ReWeight>>hWeightTemp", cutsTemp && intensityBin, "goff");
        double sumOfWeights = hWeightTemp->GetMean() * hWeightTemp->GetEntries();
        double sumOfWeights_squared = sumOfWeights * sumOfWeights;
        cleanTree->Draw("ReWeight*ReWeight>>hWeightTemp", cutsTemp && intensityBin, "goff");
        double sumOfSquaredWeights = hWeightTemp->GetMean() * hWeightTemp->GetEntries();
        double N_messy = hMessy->GetBinContent(i);
        double N_clean = hClean->GetBinContent(i);
        double x = hMessy->GetBinCenter(i);
        double y = (N_clean != 0) ? N_messy / N_clean : 0.0;
        double w = (sumOfWeights_squared != 0) ? sumOfSquaredWeights / sumOfWeights_squared : 0.0;
        double eyHigh = 0, eyLow = 0;
        if (!(y == 0 && w == 0)) {
            double midWilson = (y + 0.5 * w) / (1 + w);
            double term = TMath::Sqrt(y * (1 - y) * w + 0.25 * w * w) / (1 + w);
            eyHigh = midWilson + term - y;
            eyLow  = y - (midWilson - term);
        }
        gr->SetPoint(i - 1, x, y);
        gr->SetPointError(i - 1, 0, 0, eyLow, eyHigh);
    }
    delete hWeightTemp;
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    gr->Draw("AP");
    TString saveName = "D2_occ/" + fileTag;
    c1->SaveAs(saveName + ".png");
    c1->SaveAs(saveName + ".eps");
    TFile* outFile = new TFile(saveName + ".root", "RECREATE");
    gr->Write("eff_D2_vs_xF");
    outFile->Close();
    delete outFile;
    delete c1;
    delete gr;
}