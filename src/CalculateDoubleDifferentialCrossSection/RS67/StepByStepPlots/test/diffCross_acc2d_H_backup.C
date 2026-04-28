#include "cuts.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>   // For writing to files
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::upper_bound

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
#include "TGraphErrors.h" 

/**
 * @brief Helper function to find the mass bin index for a given mass value.
 * @param mass The mass value.
 * @param nMassBins The total number of mass bins.
 * @param mass_bins The array defining the bin edges.
 * @return The index of the mass bin, or -1 if not in any bin.
 */
int find_mass_bin(double mass, const int nMassBins, const Double_t* mass_bins) {
    auto it = std::upper_bound(mass_bins, mass_bins + nMassBins + 1, mass);
    int index = std::distance(mass_bins, it) - 1;
    if (index >= 0 && index < nMassBins) {
        return index;
    }
    return -1; 
}

/**
 * @brief Helper function to populate a 2D vector with mass values from a TTree.
 */
void FillMassVectors(TTree* tree, TCut cuts, const int nMassBins, const Double_t* mass_bins, std::vector<std::vector<double>>& mass_vectors) {
    if (!tree) return;

    tree->Draw(">>elist", cuts);
    TEventList* elist = (TEventList*)gDirectory->Get("elist");
    if (!elist) return;

    Float_t mass;
    tree->SetBranchAddress("mass", &mass);

    for (Long64_t i = 0; i < elist->GetN(); ++i) {
        Long64_t entry = elist->GetEntry(i);
        tree->GetEntry(entry);
        
        int bin_idx = find_mass_bin(mass, nMassBins, mass_bins);
        if (bin_idx != -1) {
            mass_vectors[bin_idx].push_back(mass);
        }
    }
    delete elist; 
}

void diffCross_acc2d_H_backup(Int_t xFbinNum) {
    // =========================================================================
    // ## Binning and Constant Definitions
    // =========================================================================
    const int nMassBins = 11;
    const int nXFBins = 16;

    const Double_t mass_bins[nMassBins + 1] = {
        4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7
    };

    const Double_t xf_bins[nXFBins + 1] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8
    };

    const Double_t PROTONS_ON_TARGET_LH2 = 1.57319e+17;
    const Double_t PROTONS_ON_TARGET_FLASK = 3.57904e+16;
    const Double_t flaskNorm = PROTONS_ON_TARGET_LH2 / PROTONS_ON_TARGET_FLASK;
    const Double_t LH2_TARGET_DENSITY_MOL_CM2 = 3.597;
    const Double_t AVOGADRO_CONSTANT = 6.022e23;
    const Double_t PROTONS_PER_NUCLEON_LH2 = 1.008;
    const Double_t XF_BIN_WIDTH = 0.05;
    const Double_t BEAM_ATTENUATION = (((52.0) / ((0.0708) * (50.8))) * (1 - exp(-(50.8 * 0.0708) / (52.0))));
    
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
    // ## Event-by-Event Mass Extraction for Bin Average Calculation
    // =========================================================================
    std::cout << "Extracting event-by-event mass values for xF bin " << xFbinNum << "..." << std::endl;
    std::vector<std::vector<double>> data_masses(nMassBins);
    std::vector<std::vector<double>> mix_masses(nMassBins);
    std::vector<std::vector<double>> flask_masses(nMassBins);
    std::vector<std::vector<double>> flaskmix_masses(nMassBins);

    TCut totalCut = dataCut && xFCut;
    FillMassVectors(dataTree, totalCut, nMassBins, mass_bins, data_masses);
    FillMassVectors(mixTree, totalCut, nMassBins, mass_bins, mix_masses);
    FillMassVectors(flaskTree, totalCut, nMassBins, mass_bins, flask_masses);
    FillMassVectors(flaskMixTree, totalCut, nMassBins, mass_bins, flaskmix_masses);
    
    std::vector<double> avg_masses(nMassBins);
    
    TString csv_filename = "average_masses_xF" + binNumStr + ".csv";
    std::ofstream outfile(csv_filename.Data());
    outfile << "MassBinIndex,MassBinCenter,AverageMass,N_Data,N_Mix,N_Flask,N_FlaskMix,Numerator,Denominator\n";

    for(int i = 0; i < nMassBins; ++i) {
        double sum_data = std::accumulate(data_masses[i].begin(), data_masses[i].end(), 0.0);
        double sum_mix = std::accumulate(mix_masses[i].begin(), mix_masses[i].end(), 0.0);
        double sum_flask = std::accumulate(flask_masses[i].begin(), flask_masses[i].end(), 0.0);
        double sum_flaskmix = std::accumulate(flaskmix_masses[i].begin(), flaskmix_masses[i].end(), 0.0);

        long n_data = data_masses[i].size();
        long n_mix = mix_masses[i].size();
        long n_flask = flask_masses[i].size();
        long n_flaskmix = flaskmix_masses[i].size();

        double numerator = sum_data - sum_mix - flaskNorm * (sum_flask - sum_flaskmix);
        double denominator = n_data - n_mix - flaskNorm * (n_flask - n_flaskmix);

        if (denominator != 0) {
            avg_masses[i] = numerator / denominator;
        } else {
            avg_masses[i] = (mass_bins[i] + mass_bins[i+1]) / 2.0; 
        }
        
        double bin_center = (mass_bins[i] + mass_bins[i+1]) / 2.0;
        outfile << i << "," << bin_center << "," << avg_masses[i] << "," << n_data << "," << n_mix << "," << n_flask << "," << n_flaskmix << "," << numerator << "," << denominator << "\n";
    }
    outfile.close();
    std::cout << "Average mass calculation complete. Results saved to " << csv_filename.Data() << std::endl;

    // =========================================================================
    // ## Data Extraction and Filling (for Cross-Section Calculation)
    // =========================================================================
    auto *hData_stat = new TH1D("hData_stat", "Signal Dimuons (Stat Errors)", nMassBins, mass_bins);
    hData_stat->Sumw2();
    auto *hDataMix = new TH1D("hDataMix", "Mixed-Event Background", nMassBins, mass_bins);
    hDataMix->Sumw2();
    auto *hFlask = new TH1D("hFlask", "Empty Flask Dimuons", nMassBins, mass_bins);
    hFlask->Sumw2();
    auto *hFlaskMix = new TH1D("hFlaskMix", "Empty Flask Mixed-Event Background", nMassBins, mass_bins);
    hFlaskMix->Sumw2();

    const TCut weight = "mass^3";
    dataTree->Draw("mass>>hData_stat", weight * totalCut);
    mixTree->Draw("mass>>hDataMix", weight * totalCut);
    flaskTree->Draw("mass>>hFlask", weight * totalCut);
    flaskMixTree->Draw("mass>>hFlaskMix", weight * totalCut);

    // =========================================================================
    // ## Background Subtraction
    // =========================================================================
    hData_stat->Add(hDataMix, -1.0);
    hFlask->Scale(flaskNorm);
    hFlaskMix->Scale(flaskNorm);
    hData_stat->Add(hFlask, -1.0);
    hData_stat->Add(hFlaskMix, 1.0);

    // =========================================================================
    // ## Corrections and Normalization
    // =========================================================================
    hData_stat->Divide(hAcceptance);
    hData_stat->Multiply(hKEff);
    hData_stat->Multiply(hTrigEff);
    hData_stat->Scale(constant);

    for (int k = 0; k < nMassBins; ++k) {
        double binWidth = mass_bins[k+1] - mass_bins[k];
        hData_stat->SetBinContent(k + 1, hData_stat->GetBinContent(k + 1) / binWidth);
        hData_stat->SetBinError(k + 1, hData_stat->GetBinError(k + 1) / binWidth);
    }

    // =========================================================================
    // ## Systematic Uncertainty Calculation & TH1D Creation
    // =========================================================================
    auto* hData_syst = (TH1D*)hData_stat->Clone("hData_syst");
    hData_syst->SetTitle("Signal Dimuons (Syst Errors)");

    std::vector<double> totalSystematicUncertainty(nMassBins, 0.0);
    for (int i = 0; i < nMassBins; ++i) {
        double dataContent = hData_stat->GetBinContent(i + 1);
        if (dataContent == 0) continue;

        double sysAcceptance = (hAcceptance->GetBinError(i + 1) / hAcceptance->GetBinContent(i + 1)) * dataContent;
        double sysTrigEff = (hTrigEff->GetBinError(i + 1) / hTrigEff->GetBinContent(i + 1)) * dataContent;
        double syskEff = (hKEff->GetBinError(i + 1) / hKEff->GetBinContent(i + 1)) * dataContent;

        totalSystematicUncertainty[i] = sqrt(pow(sysAcceptance, 2) + pow(sysTrigEff, 2) + pow(syskEff, 2));
        hData_syst->SetBinError(i + 1, totalSystematicUncertainty[i]);
    }
    
    // =========================================================================
    // ## Create TGraphs for Plotting at Average Mass (Stat & Syst separately)
    // =========================================================================
    auto* gData_stat = new TGraphErrors(nMassBins);
    gData_stat->SetName("gData_avgMass_stat");
    
    auto* gData_syst = new TGraphErrors(nMassBins);
    gData_syst->SetName("gData_avgMass_syst");

    for(int i = 0; i < nMassBins; ++i) {
        double crossSection = hData_stat->GetBinContent(i + 1);
        double statError = hData_stat->GetBinError(i + 1);
        double systError = totalSystematicUncertainty[i];

        if (i == 0) {
            gData_stat->SetPoint(i, avg_masses[i], 0.0);
            gData_stat->SetPointError(i, 0, 0.0);
            
            gData_syst->SetPoint(i, avg_masses[i], 0.0);
            gData_syst->SetPointError(i, 0, 0.0);
        } else {
            gData_stat->SetPoint(i, avg_masses[i], crossSection);
            gData_stat->SetPointError(i, 0, statError); // Outer error
            
            gData_syst->SetPoint(i, avg_masses[i], crossSection);
            gData_syst->SetPointError(i, 0, systError); // Inner error
        }
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

    // Draw the histogram invisibly to set up the axes and title
    hData_stat->GetXaxis()->SetRangeUser(4.2, 8.7);
    hData_stat->GetYaxis()->SetRangeUser(0.001, 3.0);
    hData_stat->SetMarkerSize(0);
    hData_stat->SetLineColor(kWhite);
    hData_stat->GetXaxis()->SetTitle("M (GeV)");
    hData_stat->GetXaxis()->CenterTitle(true);
    hData_stat->GetYaxis()->SetTitle("M^{3} d^{2}#sigma/dMdx_{F} (nb GeV^{2})");
    hData_stat->GetYaxis()->CenterTitle(true);
    hData_stat->SetTitle(Form("Double differential Cross-Section for %.2f #leq x_{F} < %.2f", xf_bins[xFbinNum], xf_bins[xFbinNum + 1]));
    hData_stat->Draw("E1");

    // Theoretical Predictions
    TFile *nnpdf4_file = new TFile("/root/github/e906-development/ROOTFiles/Hugo/NNPDF40_xFnew_p.root", "READ");
    auto *nnpdf4 = (TGraphAsymmErrors*)nnpdf4_file->Get("gr_xFbin" + binNumStr);
    nnpdf4->SetLineColor(kBlue + 2);
    nnpdf4->SetFillColorAlpha(kAzure + 1, 0.5);
    nnpdf4->SetFillStyle(3002);
    nnpdf4->Draw("L3 same"); 
    
    TFile *ct18_file = new TFile("CT18_xFnew_p_1sigma.root", "READ");
    auto *ct18 = (TGraphAsymmErrors*)ct18_file->Get("gr_xFbin" + binNumStr);
    ct18->SetLineColor(kGreen + 2);
    ct18->SetFillColorAlpha(kGreen - 5, 0.5);
    ct18->SetFillStyle(3002);
    ct18->Draw("L3 same");

    // Draw inner systematic error bar (Thick line, no marker)
    gData_syst->SetLineColor(kRed);
    gData_syst->SetLineWidth(3); 
    gData_syst->Draw("Z same"); // "Z" prevents horizontal caps from drawing

    // Draw outer statistical error bar (Standard line, with marker)
    gData_stat->SetMarkerColor(kRed);
    gData_stat->SetMarkerStyle(20);
    gData_stat->SetMarkerSize(1);
    gData_stat->SetLineColor(kRed);
    gData_stat->SetLineWidth(1);
    gData_stat->Draw("P Z same");

    // Legend
    auto legend = new TLegend(0.6, 0.7, 0.85, 0.85);
    legend->SetBorderSize(0);
    legend->AddEntry(gData_stat, "E906 (stat. err)", "lep");
    legend->AddEntry(gData_syst, "E906 (syst. err)", "l"); 
    legend->AddEntry(nnpdf4, "NNPDF 4.0", "lf");
    legend->AddEntry(ct18, "CT18", "lf");
    legend->Draw();

    TLatex *prelim = new TLatex();
    prelim->SetNDC();
    prelim->SetTextColor(kBlue);
    prelim->SetTextAlign(31);
    prelim->SetTextSize(0.04);
    prelim->DrawLatex(0.88, 0.91, "SeaQuest Preliminary");

    TImage *img = TImage::Open("SeaQuestLogo.png");
    if (img && img->IsValid()) {
        TPad* logoPad = new TPad("logoPad", "logoPad", 0.25, 0.71, 0.35, 0.79);
        logoPad->SetFillStyle(0);
        logoPad->SetBorderSize(0);
        logoPad->Draw();
        logoPad->cd();
        img->SetConstRatio(kTRUE);
        img->Draw("");
        canvas->cd();
    } else {
        std::cout << "Warning: Could not open 'SeaQuestLogo.png'. Make sure the file is in the correct directory." << std::endl;
    }

    TLatex *lumi_note = new TLatex();
    lumi_note->SetNDC(kFALSE);
    lumi_note->SetTextColor(kBlack);
    lumi_note->SetTextAlign(13);
    lumi_note->SetTextSize(0.025);
    lumi_note->DrawLatex(4.3, 1.5e-3, "10% global uncertainty due to the integrated luminosity is not included in the error bands");

    // =========================================================================
    // ## Output and Cleanup
    // =========================================================================
    canvas->SaveAs("plots/LH2_" + binNumStr + "_roofit.pdf");
    canvas->SaveAs("plots/LH2_" + binNumStr + "_avgmass.root");

    // Save the new TH1Ds and TGraphs to a distinct ROOT file
    TString outFilename = TString::Format("result_rootfiles/CrossSections_xF%s.root", binNumStr.Data());
    auto *outputFile = new TFile(outFilename, "RECREATE");
    
    hData_stat->Write("hCrossSection_StatErr");
    hData_syst->Write("hCrossSection_SystErr");
    gData_stat->Write("gCrossSection_avgMass_StatErr");
    gData_syst->Write("gCrossSection_avgMass_SystErr");

    outputFile->Close();
    delete outputFile;

    delete nnpdf4_file;
    delete ct18_file;
    delete hData_stat;
    delete hData_syst;
    delete hDataMix;
    delete hFlask;
    delete hFlaskMix;
    delete gData_stat;
    delete gData_syst;
    delete acceptanceFile;
    delete kEffFile;
    delete trigEffFile;
    delete dataFile;
    delete flaskFile;
    delete canvas;
    delete legend;
}