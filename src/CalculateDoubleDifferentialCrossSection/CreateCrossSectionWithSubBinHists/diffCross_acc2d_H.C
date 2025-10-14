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
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

/**
 * @brief Helper function to find the mass bin index for a given mass value.
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
    if (!elist) {
        std::cerr << "Error: Could not create event list for tree " << tree->GetName() << std::endl;
        return;
    }
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

void diffCross_acc2d_H(Int_t xFbinNum) {
    TH1::AddDirectory(kFALSE);

    // =========================================================================
    // ## Binning and Constant Definitions
    // =========================================================================
    const int nMassBins = 11;
    const int nXFBins = 16;
    const int nSubBins = 10;

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
    // ## Event-by-Event Mass Extraction
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
    
    std::vector<TH1D*> hData_sub(nMassBins), hMix_sub(nMassBins), hFlask_sub(nMassBins), hFlaskMix_sub(nMassBins);
    
    std::vector<std::vector<double>> data_sums(nMassBins, std::vector<double>(nSubBins, 0.0));
    std::vector<std::vector<double>> mix_sums(nMassBins, std::vector<double>(nSubBins, 0.0));
    std::vector<std::vector<double>> flask_sums(nMassBins, std::vector<double>(nSubBins, 0.0));
    std::vector<std::vector<double>> flaskmix_sums(nMassBins, std::vector<double>(nSubBins, 0.0));

    gStyle->SetOptStat(0);

    for(int i = 0; i < nMassBins; ++i) {
        double sum_data = std::accumulate(data_masses[i].begin(), data_masses[i].end(), 0.0);
        double sum_mix = std::accumulate(mix_masses[i].begin(), mix_masses[i].end(), 0.0);
        double sum_flask = std::accumulate(flask_masses[i].begin(), flask_masses[i].end(), 0.0);
        double sum_flaskmix = std::accumulate(flaskmix_masses[i].begin(), flaskmix_masses[i].end(), 0.0);
        long n_data = data_masses[i].size(); long n_mix = mix_masses[i].size(); long n_flask = flask_masses[i].size(); long n_flaskmix = flaskmix_masses[i].size();
        double numerator = sum_data - sum_mix - flaskNorm * (sum_flask - sum_flaskmix);
        double denominator = n_data - n_mix - flaskNorm * (n_flask - n_flaskmix);
        avg_masses[i] = (denominator != 0) ? (numerator / denominator) : ((mass_bins[i] + mass_bins[i+1]) / 2.0);
        double bin_center = (mass_bins[i] + mass_bins[i+1]) / 2.0;
        outfile << i << "," << bin_center << "," << avg_masses[i] << "," << n_data << "," << n_mix << "," << n_flask << "," << n_flaskmix << "," << numerator << "," << denominator << "\n";
        
        hData_sub[i] = new TH1D(TString::Format("hData_sub_bin%d_xF%d", i, xFbinNum), "", nSubBins, mass_bins[i], mass_bins[i+1]);
        for (double val : data_masses[i]) { hData_sub[i]->Fill(val); int sub_bin = hData_sub[i]->FindBin(val); if(sub_bin > 0 && sub_bin <= nSubBins) data_sums[i][sub_bin-1] += val;}
        hData_sub[i]->SetLineColor(kRed); hData_sub[i]->SetLineWidth(2);

        hMix_sub[i] = new TH1D(TString::Format("hMix_sub_bin%d_xF%d", i, xFbinNum), "", nSubBins, mass_bins[i], mass_bins[i+1]);
        for (double val : mix_masses[i]) { hMix_sub[i]->Fill(val); int sub_bin = hMix_sub[i]->FindBin(val); if(sub_bin > 0 && sub_bin <= nSubBins) mix_sums[i][sub_bin-1] += val; }
        hMix_sub[i]->SetLineColor(kBlue); hMix_sub[i]->SetLineWidth(2);

        hFlask_sub[i] = new TH1D(TString::Format("hFlask_sub_bin%d_xF%d", i, xFbinNum), "", nSubBins, mass_bins[i], mass_bins[i+1]);
        for (double val : flask_masses[i]) { hFlask_sub[i]->Fill(val); int sub_bin = hFlask_sub[i]->FindBin(val); if(sub_bin > 0 && sub_bin <= nSubBins) flask_sums[i][sub_bin-1] += val;}
        hFlask_sub[i]->SetLineColor(kGreen + 2); hFlask_sub[i]->SetLineWidth(2);

        hFlaskMix_sub[i] = new TH1D(TString::Format("hFlaskMix_sub_bin%d_xF%d", i, xFbinNum), "", nSubBins, mass_bins[i], mass_bins[i+1]);
        for (double val : flaskmix_masses[i]) { hFlaskMix_sub[i]->Fill(val); int sub_bin = hFlaskMix_sub[i]->FindBin(val); if(sub_bin > 0 && sub_bin <= nSubBins) flaskmix_sums[i][sub_bin-1] += val;}
        hFlaskMix_sub[i]->SetLineColor(kOrange + 7); hFlaskMix_sub[i]->SetLineWidth(2);
    }
    outfile.close();
    std::cout << "Average mass calculation and sub-histogram filling complete." << std::endl;

    // =========================================================================
    // ## Plotting of Sub-Histograms with Sum(mass) and Counts Text
    // =========================================================================
    std::cout << "Generating sub-histogram plots..." << std::endl;
    for (int i = 0; i < nMassBins; ++i) {
        auto *subCanvas = new TCanvas(TString::Format("subCanvas_bin%d_xF%d", i, xFbinNum), "", 800, 600);
        subCanvas->SetTickx(1); subCanvas->SetTicky(1);
        double max_val = 0;
        if (hData_sub[i]->GetMaximum() > max_val) max_val = hData_sub[i]->GetMaximum();
        if (hMix_sub[i]->GetMaximum() > max_val) max_val = hMix_sub[i]->GetMaximum();
        if (hFlask_sub[i]->GetMaximum() > max_val) max_val = hFlask_sub[i]->GetMaximum();
        if (hFlaskMix_sub[i]->GetMaximum() > max_val) max_val = hFlaskMix_sub[i]->GetMaximum();
        
        double y_max = max_val > 0 ? max_val * 1.5 : 1.0; // Increased padding for multi-line text
        hData_sub[i]->SetMaximum(y_max);
        hData_sub[i]->GetXaxis()->SetTitle("M (GeV)"); hData_sub[i]->GetXaxis()->CenterTitle(true);
        hData_sub[i]->GetYaxis()->SetTitle("Counts / Bin"); hData_sub[i]->GetYaxis()->CenterTitle(true);
        hData_sub[i]->SetTitle(TString::Format("Mass Dist. in Bin %d (%.2f-%.2f GeV) for x_{F} Bin %d", i, mass_bins[i], mass_bins[i+1], xFbinNum));
        
        hData_sub[i]->Draw("hist");
        if (hMix_sub[i]->GetEntries() > 0) hMix_sub[i]->Draw("hist SAME");
        if (hFlask_sub[i]->GetEntries() > 0) hFlask_sub[i]->Draw("hist SAME");
        if (hFlaskMix_sub[i]->GetEntries() > 0) hFlaskMix_sub[i]->Draw("hist SAME");

        // --- MODIFICATION START ---
        TLatex* latex = new TLatex();
        latex->SetTextSize(0.025);
        latex->SetTextAlign(23); // Align to Center-Top of the coordinate
        double y_offset = y_max * 0.025; // Small vertical offset from the top of the bin

        // Data
        latex->SetTextColor(kRed);
        for(int j=1; j <= nSubBins; ++j) {
            double count = hData_sub[i]->GetBinContent(j);
            if(count > 0) {
                double sum = data_sums[i][j-1];
                latex->DrawLatex(hData_sub[i]->GetBinCenter(j), count + y_offset, TString::Format("#splitline{%.3f}{(%.0f)}", sum, count));
            }
        }
        // Mix
        latex->SetTextColor(kBlue);
        for(int j=1; j <= nSubBins; ++j) {
            double count = hMix_sub[i]->GetBinContent(j);
            if(count > 0) {
                double sum = mix_sums[i][j-1];
                latex->DrawLatex(hMix_sub[i]->GetBinCenter(j), count + y_offset, TString::Format("#splitline{%.3f}{(%.0f)}", sum, count));
            }
        }
        // Flask
        latex->SetTextColor(kGreen + 2);
        for(int j=1; j <= nSubBins; ++j) {
            double count = hFlask_sub[i]->GetBinContent(j);
            if(count > 0) {
                double sum = flask_sums[i][j-1];
                latex->DrawLatex(hFlask_sub[i]->GetBinCenter(j), count + y_offset, TString::Format("#splitline{%.3f}{(%.1f)}", sum * flaskNorm, count * flaskNorm));
            }
        }
        // FlaskMix
        latex->SetTextColor(kOrange + 7);
        for(int j=1; j <= nSubBins; ++j) {
            double count = hFlaskMix_sub[i]->GetBinContent(j);
            if(count > 0) {
                double sum = flaskmix_sums[i][j-1];
                latex->DrawLatex(hFlaskMix_sub[i]->GetBinCenter(j), count + y_offset, TString::Format("#splitline{%.3f}{(%.1f)}", sum * flaskNorm, count * flaskNorm));
            }
        }
        delete latex;
        // --- MODIFICATION END ---

        auto subLegend = new TLegend(0.6, 0.75, 0.88, 0.88);
        subLegend->SetBorderSize(0);
        subLegend->AddEntry(hData_sub[i], "Data", "l");
        subLegend->AddEntry(hMix_sub[i], "Mix", "l");
        subLegend->AddEntry(hFlask_sub[i], "Flask", "l");
        subLegend->AddEntry(hFlaskMix_sub[i], "Empty Flask Mix", "l");
        subLegend->Draw();
        subCanvas->SaveAs(TString::Format("plots/sub_mass_dist_xF%d_massBin%d.pdf", xFbinNum, i));
        delete subLegend;
        delete subCanvas;
    }
    std::cout << "Sub-histogram plots generated." << std::endl;
    
    // =========================================================================
    // ## Data Extraction and Filling (for Cross-Section Calculation)
    // =========================================================================
    auto *hData = new TH1D("hData", "Signal Dimuons", nMassBins, mass_bins);
    hData->Sumw2();
    auto *hDataMix = new TH1D("hDataMix", "Mixed-Event Background", nMassBins, mass_bins);
    hDataMix->Sumw2();
    auto *hFlask = new TH1D("hFlask", "Empty Flask Dimuons", nMassBins, mass_bins);
    hFlask->Sumw2();
    auto *hFlaskMix = new TH1D("hFlaskMix", "Empty Flask Mixed-Event Background", nMassBins, mass_bins);
    hFlaskMix->Sumw2();
    const TCut weight = "mass^3";
    dataTree->Draw("mass>>hData", weight * totalCut);
    mixTree->Draw("mass>>hDataMix", weight * totalCut);
    flaskTree->Draw("mass>>hFlask", weight * totalCut);
    flaskMixTree->Draw("mass>>hFlaskMix", weight * totalCut);

    // =========================================================================
    // ## Background Subtraction and Corrections
    // =========================================================================
    hData->Add(hDataMix, -1.0);
    hFlask->Scale(flaskNorm);
    hFlaskMix->Scale(flaskNorm);
    hData->Add(hFlask, -1.0);
    hData->Add(hFlaskMix, 1.0);
    hData->Divide(hAcceptance);
    hData->Multiply(hKEff);
    hData->Multiply(hTrigEff);
    hData->Scale(constant);
    for (int k = 0; k < nMassBins; ++k) {
        double binWidth = mass_bins[k+1] - mass_bins[k];
        hData->SetBinContent(k + 1, hData->GetBinContent(k + 1) / binWidth);
        hData->SetBinError(k + 1, hData->GetBinError(k + 1) / binWidth);
    }

    // =========================================================================
    // ## Systematic Uncertainty Calculation and Final Plotting
    // =========================================================================
    auto* hDataOverlay = (TH1D*)hData->Clone("hDataOverlay");
    for (int i = 0; i < nMassBins; ++i) {
        double dataContent = hData->GetBinContent(i + 1);
        double sysError = 0.0;
        if (dataContent != 0 && hAcceptance->GetBinContent(i+1) != 0 && hTrigEff->GetBinContent(i+1) != 0 && hKEff->GetBinContent(i+1) != 0) {
            double sysAcceptance = (hAcceptance->GetBinError(i + 1) / hAcceptance->GetBinContent(i + 1));
            double sysTrigEff = (hTrigEff->GetBinError(i + 1) / hTrigEff->GetBinContent(i + 1));
            double syskEff = (hKEff->GetBinError(i + 1) / hKEff->GetBinContent(i + 1));
            sysError = dataContent * sqrt(pow(sysAcceptance, 2) + pow(sysTrigEff, 2) + pow(syskEff, 2));
        }
        hDataOverlay->SetBinError(i + 1, sysError);
    }
    
    auto* gData = new TGraphErrors(nMassBins);
    for(int i = 0; i < nMassBins; ++i) {
        gData->SetPoint(i, avg_masses[i], hData->GetBinContent(i + 1));
        gData->SetPointError(i, 0, hData->GetBinError(i + 1));
    }
    if (nMassBins > 0) { gData->SetPoint(0, avg_masses[0], 0.0); gData->SetPointError(0, 0, 0.0); }
    
    auto *canvas = new TCanvas("canvas", "Cross-Section", 800, 600);
    canvas->SetLogy(); canvas->SetTickx(1); canvas->SetTicky(1);

    hData->GetXaxis()->SetRangeUser(4.2, 8.7); hData->GetYaxis()->SetRangeUser(0.001, 3.0);
    hData->SetMarkerSize(0); hData->SetLineColor(kWhite);
    hData->GetXaxis()->SetTitle("M (GeV)"); hData->GetXaxis()->CenterTitle(true);
    hData->GetYaxis()->SetTitle("M^{3} d^{2}#sigma/dMdx_{F} (nb GeV^{2})"); hData->GetYaxis()->CenterTitle(true);
    hData->SetTitle(Form("Double differential Cross-Section for %.2f #leq x_{F} < %.2f", xf_bins[xFbinNum], xf_bins[xFbinNum + 1]));
    hData->Draw("E1");

    TFile *nnpdf4_file = new TFile("/root/github/e906-development/ROOTFiles/Hugo/NNPDF40_xFnew_p.root", "READ");
    auto *nnpdf4 = (TGraphAsymmErrors*)nnpdf4_file->Get("gr_xFbin" + binNumStr);
    nnpdf4->SetLineColor(kBlue + 2); nnpdf4->SetFillColorAlpha(kAzure + 1, 0.5); nnpdf4->SetFillStyle(3002);
    nnpdf4->Draw("L3 same"); 
    
    TFile *ct18_file = new TFile("CT18_xFnew_p_1sigma.root", "READ");
    auto *ct18 = (TGraphAsymmErrors*)ct18_file->Get("gr_xFbin" + binNumStr);
    ct18->SetLineColor(kGreen + 2); ct18->SetFillColorAlpha(kGreen - 5, 0.5); ct18->SetFillStyle(3002);
    ct18->Draw("L3 same");

    hDataOverlay->SetMarkerSize(0); hDataOverlay->SetLineColor(kRed); hDataOverlay->SetFillColorAlpha(kPink - 9, 0.5);
    if(nMassBins > 0) {hDataOverlay->SetBinContent(1, 0.0); hDataOverlay->SetBinError(1, 0.0);}
    hDataOverlay->Draw("E2 same");

    gData->SetMarkerColor(kRed); gData->SetMarkerStyle(20); gData->SetMarkerSize(1); gData->SetLineColor(kRed);
    gData->Draw("P same");

    auto legend = new TLegend(0.6, 0.7, 0.85, 0.85);
    legend->SetBorderSize(0);
    legend->AddEntry(gData, "E906 (stat. only)", "lep");
    legend->AddEntry(hDataOverlay, "E906 (syst. only)", "f");
    legend->AddEntry(nnpdf4, "NNPDF 4.0", "lf");
    legend->AddEntry(ct18, "CT18", "lf");
    legend->Draw();

    // =========================================================================
    // ## Output and Cleanup
    // =========================================================================
    canvas->SaveAs("plots/LH2_" + binNumStr + "_roofit.pdf");
    canvas->SaveAs("plots/LH2_" + binNumStr + "_avgmass.root");

    delete nnpdf4_file; delete ct18_file;
    delete hData; delete hDataMix; delete hFlask; delete hFlaskMix;
    delete hDataOverlay; delete gData;
    delete acceptanceFile; delete kEffFile; delete trigEffFile;
    delete dataFile; delete flaskFile;
    delete canvas; delete legend;
    for(int i = 0; i < nMassBins; ++i) {
        delete hData_sub[i]; delete hMix_sub[i];
        delete hFlask_sub[i]; delete hFlaskMix_sub[i];
    }
}