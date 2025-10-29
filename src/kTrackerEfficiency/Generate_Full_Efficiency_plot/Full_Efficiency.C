/**
 * @file Full_Efficiency.cpp
 * @brief Calculates and plots the D2-dependent occupancy efficiency.
 *
 * This script compares a "clean" (no occupancy effects) Monte Carlo sample
 * with a "messy" (with occupancy effects) Monte Carlo sample to determine
 * the reconstruction efficiency as a function of the D2 variable.
 *
 * It calculates the efficiency as (Weighted Messy Events) / (Weighted Clean Events)
 * and computes asymmetric error bars using the weighted Wilson score interval method.
 * The calculation is optimized to fill histograms once, rather than looping
 * TTree::Draw() calls.
 */

#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCut.h>
#include <iostream>
#include "chuckcuts.h" // Assumed to contain TCut definitions

// --- Using namespace ---
using namespace std;

// ---------------- Input files ----------------
/** @var messyFileName Path to the "messy" MC file (with occupancy effects). */
TString messyFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_messy_occ_pTxFweight_v2.root";
/** @var cleanFileName Path to the "clean" MC file (no occupancy effects). */
TString cleanFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S002_clean_occ_pTxFweight_v2.root";

// ---------------- Binning ----------------
/** @var binVar The variable to bin against (e.g., "D2"). */
TString binVar = "D2";
/** @var binLow The lower edge of the first bin. */
double binLow = 0.0;
/** @var binHigh The upper edge of the last bin. */
double binHigh = 400.0;
/** @var binCount The total number of bins. */
int binCount = 16;

// ---------------- Mass bins (Originals, currently unused in loop) ----------------
/** @var NBINS The number of mass bins defined. */
const int NBINS = 11;
/** @var edges The mass bin edges. */
Double_t edges[NBINS + 1] = {4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7};

/**
 * @brief Gets the label for the integrated mass bin.
 * @return TString A formatted label for plots.
 */
TString get_mass_label() {
    return "4.2 #leq M < 8.7 and 0.0 #leq x_{F} < 0.85";
}

/**
 * @brief Gets the TCut for the integrated mass bin.
 * @return TCut A TCut object for the mass and xF range.
 */
TCut get_mass_cut() {
    return "mass >= 4.2 && mass < 8.7 && xF >= 0.0 && xF < 0.85";
}


/**
 * @brief Calculates efficiency, computes errors, plots, and saves the result.
 *
 * This function replaces the original plotAbsolute and getErrors.
 * It fills three histograms:
 * 1. hClean_sumW: Denominator (clean MC), weighted by "ReWeight" (Sum of Weights).
 * 2. hClean_sumW2: Denominator (clean MC), weighted by "ReWeight*ReWeight" (Sum of Weights^2).
 * 3. hMessy_sumW: Numerator (messy MC), weighted by "ReWeight" (Sum of Weights).
 *
 * It then loops over bins to calculate the efficiency (y = messy/clean) and
 * the weighted Wilson score errors, which depend on y and w_eff = (SumW2_clean) / (SumW_clean)^2.
 *
 * @param cleanTree Pointer to the TTree from the "clean" file.
 * @param messyTree Pointer to the TTree from the "messy" file.
 * @param CcutsTemp The selection cuts to apply to the clean tree.
 * @param McutsTemp The selection cuts to apply to the messy tree.
 * @param outTag A string tag (based on mass bin) for output file and plot titles.
 */
void calculateAndPlotEfficiency(TTree* cleanTree, TTree* messyTree, TCut CcutsTemp, TCut McutsTemp, TString outTag) {
    
    // --- 1. Define Histograms ---
    // We need three histograms:
    // 1. Clean (denominator) sum of weights (SumW)
    // 2. Clean (denominator) sum of weights squared (SumW2)
    // 3. Messy (numerator) sum of weights (SumW)
    
    TH1D* hClean_sumW = new TH1D("hClean_sumW", "Clean Events (SumW)", binCount, binLow, binHigh);
    TH1D* hClean_sumW2 = new TH1D("hClean_sumW2", "Clean Events (SumW2)", binCount, binLow, binHigh);
    TH1D* hMessy_sumW = new TH1D("hMessy_sumW", "Messy Events (SumW)", binCount, binLow, binHigh);

    // --- 2. Fill Histograms ---
    // This is the optimized part: only 3 Draw calls total.
    TString cleanWeightCut = "ReWeight * (" + TString(CcutsTemp) + ")";
    TString cleanWeightSqCut = "ReWeight*ReWeight * (" + TString(CcutsTemp) + ")";
    TString messyWeightCut = "ReWeight * (" + TString(McutsTemp) + ")";

    cleanTree->Draw(binVar + ">>hClean_sumW", cleanWeightCut, "goff");
    cleanTree->Draw(binVar + ">>hClean_sumW2", cleanWeightSqCut, "goff");
    messyTree->Draw(binVar + ">>hMessy_sumW", messyWeightCut, "goff");

    // --- 3. Create Graph and Calculate Bin-by-Bin Efficiency & Errors ---
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(binCount);
    gr->SetMinimum(0);
    gr->SetMaximum(1.2);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    
    for (int i = 1; i <= binCount; i++) {
        double sumW_messy = hMessy_sumW->GetBinContent(i);   // Numerator
        double sumW_clean = hClean_sumW->GetBinContent(i);   // Denominator
        double sumW2_clean = hClean_sumW2->GetBinContent(i); // Denominator (for error)

        double x = hMessy_sumW->GetBinCenter(i);
        double y = (sumW_clean != 0) ? sumW_messy / sumW_clean : 0.0; // Efficiency

        // Effective weight parameter: w = sum(w_i^2) / (sum(w_i))^2
        // This 'w' is used in the weighted Wilson score formula.
        double w_eff = (sumW_clean != 0) ? sumW2_clean / (sumW_clean * sumW_clean) : 0.0;

        // Calculate weighted Wilson score interval
        double midWilson = (y + 0.5 * w_eff) / (1.0 + w_eff);
        double eyHigh = 0.0;
        double eyLow = 0.0;

        // Check to avoid NaNs if y=0 and w_eff=0
        if (y > 0.0 || w_eff > 0.0) {
            double term = TMath::Sqrt(y * (1.0 - y) * w_eff + 0.25 * w_eff * w_eff) / (1.0 + w_eff);
            eyHigh = (midWilson + term) - y;
            eyLow  = y - (midWilson - term);
        }

        // Clean up potential floating point inaccuracies (e.g., tiny negative errors)
        if (eyLow < 0.0) eyLow = 0.0;
        if (eyHigh < 0.0) eyHigh = 0.0;

        gr->SetPoint(i - 1, x, y);
        gr->SetPointError(i - 1, 0, 0, eyLow, eyHigh);
    }

    // --- 4. Build interpolation curve (as in original) ---
    TGraphAsymmErrors* grInterp = new TGraphAsymmErrors();
    int n = gr->GetN();
    for (int i = 0; i < n - 1; i++) {
        double x1, y1, x2, y2;
        gr->GetPoint(i, x1, y1);
        gr->GetPoint(i + 1, x2, y2);

        double eyl1 = gr->GetErrorYlow(i);
        double eyh1 = gr->GetErrorYhigh(i);
        double eyl2 = gr->GetErrorYlow(i+1);
        double eyh2 = gr->GetErrorYhigh(i+1);

        double xm = 0.5 * (x1 + x2);
        double ym = 0.5 * (y1 + y2);

        double eym_low  = 0.5 * TMath::Sqrt(eyl1 * eyl1 + eyl2 * eyl2);
        double eym_high = 0.5 * TMath::Sqrt(eyh1 * eyh1 + eyh2 * eyh2);

        grInterp->SetPoint(i, xm, ym);
        grInterp->SetPointError(i, 0, 0, eym_low, eym_high);
    }

    grInterp->SetMarkerColor(kRed);
    grInterp->SetMarkerStyle(24);
    grInterp->SetLineColor(kRed);
    grInterp->SetLineWidth(2);

    // --- 5. Draw results ---
    TCanvas* c1 = new TCanvas("c1", "Efficiency Canvas", 800, 600);
    TString title = "D2_Efficiency_" + outTag;
    gr->SetTitle(title + ";D2;Efficiency");
    gr->Draw("AP");
    // grInterp->Draw("LP SAME"); // Kept commented as in original

    c1->SaveAs("D2_occ/Efficiency_D2_Complete.pdf");

    // --- 6. Save to file ---
    TFile* outFile = new TFile("D2_occ/Efficiency_D2_Complete.root", "RECREATE");
    if (outFile->IsOpen()) {
        gr->Write("eff_original");
        // grInterp->Write("eff_interpolated"); // Kept commented
        outFile->Close();
    } else {
        cerr << "Error: Could not create output ROOT file." << endl;
    }

    // --- 7. Clean up memory ---
    delete hClean_sumW;
    delete hClean_sumW2;
    delete hMessy_sumW;
    delete gr;
    delete grInterp;
    delete c1;
    delete outFile;
}


/**
 * @brief Main driver function for the efficiency calculation.
 *
 * This function defines the base cuts, opens the input TFiles and TTrees,
 * sets up the specific cuts for the analysis (currently one integrated bin),
 * and calls calculateAndPlotEfficiency() to perform the main work.
 */
void Full_Efficiency() {
    // --- Base cuts ---
    TCut baseCcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && DYCut_MC;

    TCut baseMcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp &&
                    occCuts_2111v42 && DYCut_MC;

    // --- Open files and get Trees ---
    TFile* cleanFile = TFile::Open(cleanFileName, "READ");
    TFile* messyFile = TFile::Open(messyFileName, "READ");

    if (!cleanFile || cleanFile->IsZombie() || !messyFile || messyFile->IsZombie()) {
        cerr << "Error: Could not open one or both input files." << endl;
        if (cleanFile) cleanFile->Close();
        if (messyFile) messyFile->Close();
        delete cleanFile;
        delete messyFile;
        return;
    }

    TTree* cleanTree = (TTree*)cleanFile->Get("Tree");
    TTree* messyTree = (TTree*)messyFile->Get("Tree");

    if (!cleanTree || !messyTree) {
        cerr << "Error: Could not retrieve TTrees from files." << endl;
        cleanFile->Close();
        messyFile->Close();
        delete cleanFile;
        delete messyFile;
        return;
    }

    // --- Loop over mass bins (currently running for one integrated bin) ---
    // for (int imass=0; imass<NBINS; imass++) {
        // TCut masscut = get_mass_cut(); // Original script logic
        // TString masslabel = get_mass_label(); // Original script logic

        // Combined cuts (mass/xF cut removed as in original)
        TCut CcutsTemp = baseCcut;
        TCut McutsTemp = baseMcut;

        // Output tag (naming updated)
        TString outTag = get_mass_label();

        // Call the simplified calculation and plotting function
        calculateAndPlotEfficiency(cleanTree, messyTree, CcutsTemp, McutsTemp, outTag);
    // }

    // --- Clean up files ---
    cleanFile->Close();
    messyFile->Close();
    delete cleanFile;
    delete messyFile;
}