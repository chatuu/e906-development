#include "cuts.h"
#include <cmath>
#include <iostream>
#include <vector>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCut.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

/**
 * @file
 * @brief Main analysis script for calculating the double-differential cross-section.
 */

/**
 * @brief Calculates and plots the double-differential cross-section for a given xF bin.
 *
 * This function performs a comprehensive analysis to extract the Drell-Yan double-differential
 * cross-section, $M^3 \frac{d^2\sigma}{dMdx_F}$, for a specific Feynman-x ($x_F$) bin from
 * the E906/SeaQuest experiment's liquid hydrogen (LH2) target data.
 *
 * The process involves several key steps:
 * 1.  Loading raw dimuon event data from the main target and an empty "flask" target,
 * including both real and mixed-event background pairs.
 * 2.  Subtracting the mixed-event background and the normalized empty flask contribution.
 * 3.  Applying crucial corrections for detector acceptance, trigger efficiency, and
 * k-factor efficiency.
 * 4.  Normalizing the result by luminosity, target properties, and bin width to obtain
 * the final cross-section.
 * 5.  Calculating and displaying systematic uncertainties from the applied corrections.
 * 6.  Plotting the final cross-section data along with theoretical predictions from
 * NNPDF4.0 and CT18 parton distribution functions.
 * 7.  Saving the resulting plot and the final data histogram to output files.
 *
 * @param xFbinNum The index of the $x_F$ bin to be analyzed (from 0 to 15).
 */
void diffCross_acc2d_H(Int_t xFbinNum) {
    // =========================================================================
    // ## Binning and Constant Definitions
    // =========================================================================
    const int nMassBins = 11;
    const int nXFBins = 16;

    // Mass bins in GeV: [4.2, 4.5), [4.5, 4.8), ..., [7.5, 8.7)
    const Double_t mass_bins[nMassBins + 1] = {
        4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7
    };

    // xF bins: [0.0, 0.05), [0.05, 0.1), ..., [0.8, 0.85)
    const Double_t xf_bins[nXFBins + 1] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8
    };

    // Physics and normalization constants
    const Double_t PROTONS_ON_TARGET_LH2 = 1.57319e+17;
    const Double_t PROTONS_ON_TARGET_FLASK = 3.57904e+16;
    const Double_t flaskNorm = PROTONS_ON_TARGET_LH2 / PROTONS_ON_TARGET_FLASK;
    const Double_t LH2_TARGET_DENSITY_MOL_CM2 = 3.597; // mol/cm^2
    const Double_t AVOGADRO_CONSTANT = 6.022e23; // atoms/mol
    const Double_t PROTONS_PER_NUCLEON_LH2 = 1.008; // g/mol
    const Double_t XF_BIN_WIDTH = 0.05;
    const Double_t BEAM_ATTENUATION = (((52.0) / ((0.0708) * (50.8))) * (1 - exp(-(50.8 * 0.0708) / (52.0))));
    
    // Final cross-section scaling factor (units: nb * GeV^2)
    const Double_t constant = (PROTONS_PER_NUCLEON_LH2 * 1e33) / 
                              (LH2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LH2 * BEAM_ATTENUATION * XF_BIN_WIDTH);


    // =========================================================================
    // ## File and Histogram Loading
    // =========================================================================
    TString binNumStr = TString::Format("%d", xFbinNum);
    TCut xFCut = get_xf_cut(xFbinNum);

    auto *acceptanceFile = new TFile("acceptance_mass_xF.root", "READ");
    auto *kEffFile = new TFile("kEff_xF.root", "READ");
    auto *trigEffFile = new TFile("Trig_Eff_Correction.root", "READ");
    
    auto *hAcceptance = (TH1F*)acceptanceFile->Get("h_ratio_LH2_xF_bin" + binNumStr);
    auto *hKEff = (TH1F*)kEffFile->Get("kEff_xF_" + binNumStr);
    auto *hTrigEff = (TH1D*)trigEffFile->Get("h_eff");
    
    TFile* dataFile = new TFile("/root/github/e906-development/ROOTFiles/MixedEvents/merged_RS67_3089LH2.root", "READ");
    TTree* dataTree = (TTree*)dataFile->Get("result");
    TTree* mixTree = (TTree*)dataFile->Get("result_mix");

    TFile* flaskFile = new TFile("/root/github/e906-development/ROOTFiles/MixedEvents/merged_RS67_3089flask.root", "READ");
    TTree* flaskTree = (TTree*)flaskFile->Get("result");
    TTree* flaskMixTree = (TTree*)flaskFile->Get("result_mix");


    // =========================================================================
    // ## Data Extraction and Filling
    // =========================================================================
    // Create histograms to store data from TTrees, weighted by M^3
    auto *hData = new TH1D("hData", "Signal Dimuons", nMassBins, mass_bins);
    hData->Sumw2();
    auto *hDataMix = new TH1D("hDataMix", "Mixed-Event Background", nMassBins, mass_bins);
    hDataMix->Sumw2();
    auto *hFlask = new TH1D("hFlask", "Empty Flask Dimuons", nMassBins, mass_bins);
    hFlask->Sumw2();
    auto *hFlaskMix = new TH1D("hFlaskMix", "Empty Flask Mixed-Event Background", nMassBins, mass_bins);
    hFlaskMix->Sumw2();

    const TCut weight = "mass^3";
    dataTree->Draw("mass>>hData", weight * (dataCut && xFCut));
    mixTree->Draw("mass>>hDataMix", weight * (dataCut && xFCut));
    flaskTree->Draw("mass>>hFlask", weight * (dataCut && xFCut));
    flaskMixTree->Draw("mass>>hFlaskMix", weight * (dataCut && xFCut));


    // =========================================================================
    // ## Background Subtraction
    // =========================================================================
    hData->Add(hDataMix, -1.0);      // Subtract combinatorial background from target data
    
    hFlask->Scale(flaskNorm);        // Scale flask data by proton-on-target ratio
    hFlaskMix->Scale(flaskNorm);     // Scale flask mixed-event data
    
    hData->Add(hFlask, -1.0);        // Subtract normalized empty flask contribution
    hData->Add(hFlaskMix, 1.0);      // Add back normalized flask combinatorial background


    // =========================================================================
    // ## Corrections and Normalization
    // =========================================================================
    hData->Divide(hAcceptance);      // Correct for detector acceptance
    hData->Multiply(hKEff);          // Apply k-factor efficiency correction
    hData->Multiply(hTrigEff);       // Correct for trigger efficiency
    hData->Scale(constant);          // Scale by physics constants to get cross-section

    // Normalize by bin width
    for (int k = 0; k < nMassBins; ++k) {
        double binWidth = mass_bins[k+1] - mass_bins[k];
        hData->SetBinContent(k + 1, hData->GetBinContent(k + 1) / binWidth);
        hData->SetBinError(k + 1, hData->GetBinError(k + 1) / binWidth);
    }


    // =========================================================================
    // ## Systematic Uncertainty Calculation
    // =========================================================================
    std::vector<double> totalSystematicUncertainty(nMassBins, 0.0);
    for (int i = 0; i < nMassBins; ++i) {
        double dataContent = hData->GetBinContent(i + 1);
        if (dataContent == 0) continue;

        double sysAcceptance = (hAcceptance->GetBinError(i + 1) / hAcceptance->GetBinContent(i + 1)) * dataContent;
        double sysTrigEff = (hTrigEff->GetBinError(i + 1) / hTrigEff->GetBinContent(i + 1)) * dataContent;
        double syskEff = (hKEff->GetBinError(i + 1) / hKEff->GetBinContent(i + 1)) * dataContent;

        totalSystematicUncertainty[i] = sqrt(pow(sysAcceptance, 2) + pow(sysTrigEff, 2) + pow(syskEff, 2));
    }
    
    // Create a histogram clone to display systematic uncertainties
    auto* hDataOverlay = (TH1D*)hData->Clone("hDataOverlay");
    for (int k = 0; k < nMassBins; ++k) {
        hDataOverlay->SetBinError(k + 1, totalSystematicUncertainty[k]);
    }


    // =========================================================================
    // ## Plotting and Styling
    // =========================================================================
    auto *canvas = new TCanvas("canvas", "Cross-Section", 800, 600);
    canvas->SetLogy();
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Style main data histogram (statistical error)
    hData->GetXaxis()->SetRangeUser(4.2, 8.5);
    hData->GetYaxis()->SetRangeUser(0.001, 1.0);
    hData->SetMarkerColor(kRed);
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1);
    hData->SetLineColor(kRed);
    hData->GetXaxis()->SetTitle("M (GeV)");
    hData->GetXaxis()->CenterTitle(true);
    hData->GetYaxis()->SetTitle("M^{3} d^{2}#sigma/dMdx_{F} (nb GeV^{2})");
    hData->GetYaxis()->CenterTitle(true);
    hData->SetTitle(Form("Double differential Cross-Section for %.2f #leq x_{F} < %.2f", xf_bins[xFbinNum], xf_bins[xFbinNum + 1]));
    hData->Draw("E1");

    // Load and style theoretical prediction: NNPDF 4.0
    TFile *nnpdf4_file = new TFile("/root/github/e906-development/ROOTFiles/Hugo/NNPDF40_xFnew_p.root", "READ");
    auto *nnpdf4 = (TGraphAsymmErrors*)nnpdf4_file->Get("gr_xFbin" + binNumStr);
    nnpdf4->SetLineColor(kBlue + 2);
    nnpdf4->SetFillColorAlpha(kAzure + 1, 0.5);
    nnpdf4->SetFillStyle(3002);
    nnpdf4->Draw("L3 same"); 
    
    // Load and style theoretical prediction: CT18
    TFile *ct18_file = new TFile("/root/github/e906-development/ROOTFiles/Hugo/CT18_xFnew_p.root", "READ");
    auto *ct18 = (TGraphAsymmErrors*)ct18_file->Get("gr_xFbin" + binNumStr);
    ct18->SetLineColor(kGreen + 2);
    ct18->SetFillColorAlpha(kGreen - 5, 0.5);
    ct18->SetFillStyle(3002);
    ct18->Draw("L3 same");

    // Style and draw systematic error overlay
    hDataOverlay->SetMarkerSize(0);
    hDataOverlay->SetLineColor(kRed);
    hDataOverlay->SetFillColorAlpha(kPink - 9, 0.5);
    hDataOverlay->Draw("E2 same");

    // Legend
    auto legend = new TLegend(0.6, 0.7, 0.85, 0.85);
    legend->SetBorderSize(0);
    legend->AddEntry(hData, "E906 (stat. only)", "lep");
    legend->AddEntry(hDataOverlay, "E906 (syst. only)", "f");
    legend->AddEntry(nnpdf4, "NNPDF 4.0", "lf");
    legend->AddEntry(ct18, "CT18", "lf");
    legend->Draw();


    // =========================================================================
    // ## Output and Cleanup
    // =========================================================================
    canvas->SaveAs("plots/LH2_" + binNumStr + "_roofit.pdf");
    canvas->SaveAs("plots/LH2_" + binNumStr + "_roofit.root");

    // Save final histogram to a ROOT file
    auto *outputFile = new TFile("result_rootfiles/LH2_" + binNumStr + "_updatedPoT.root", "RECREATE");
    hData->Write("hAccCor");
    outputFile->Close();

    // Clean up dynamically allocated memory
    delete outputFile;
    delete nnpdf4_file;
    delete ct18_file;
    delete hData;
    delete hDataMix;
    delete hFlask;
    delete hFlaskMix;
    delete hDataOverlay;
    delete acceptanceFile;
    delete kEffFile;
    delete trigEffFile;
    delete dataFile;
    delete flaskFile;
    delete canvas;
    delete legend;
}