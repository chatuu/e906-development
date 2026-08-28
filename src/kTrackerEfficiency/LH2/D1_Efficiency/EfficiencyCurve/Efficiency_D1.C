#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>  // Required for writing to files
#include <iomanip>  // Required for formatting numbers (setprecision)
#include "chuckcuts.h"

using namespace std;

// ---------------- Input files ----------------
TString messyFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root";
TString cleanFileName = "~/github/e906-development/ROOTFiles/Hugo/mc_drellyan_LH2_M027_S001_clean_occ_pTxFweight_v2.root";

// ---------------- Binning ----------------
TString binVar = "D1";
double binLow = 0.0;
double binHigh = 600.0;
int binCount = 24;

// Trees & histograms
TTree *cleanTree = nullptr;
TTree *messyTree = nullptr;
TH1D *hClean = nullptr;
TH1D *hMessy = nullptr;
TGraphAsymmErrors *gr = nullptr;        // original efficiency curve
TGraphAsymmErrors *grInterp = nullptr;  // interpolated curve

// ---------------- Forward declarations ----------------
void plotAbsolute(TCut CcutsTemp, TCut McutsTemp, TString outTag);
void getErrors(TCut cutsTemp, TString outTag);

// ---------------- Mass bins ----------------
const int NBINS = 11;
Double_t edges[NBINS + 1] = {4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7};

TString get_mass_label(int iMass) {
    return Form("M%.1fto%.1f", edges[iMass], edges[iMass+1]);
}

TCut get_mass_cut(int iMass) {
    return Form("mass > %f && mass < %f", edges[iMass], edges[iMass+1]);
}

// ---------------- Main driver ----------------
void Efficiency_D1() {
    // Base cuts
    TCut baseCcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && DYCut_MC;

    TCut baseMcut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp &&
                    chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp &&
                    occCuts_2111v42 && DYCut_MC;

    // Open files
    TFile* cleanFile = new TFile(cleanFileName, "READ");
    cleanTree = (TTree*)cleanFile->Get("Tree");

    TFile* messyFile = new TFile(messyFileName, "READ");
    messyTree = (TTree*)messyFile->Get("Tree");

    // Combined cuts (xF cut removed)
    TCut CcutsTemp = baseCcut;
    TCut McutsTemp = baseMcut;

    // Histograms
    TString hCname = Form("hClean_mass_D1");
    TString hMname = Form("hMessy_mass_D1");
    hClean = new TH1D(hCname, hCname, binCount, binLow, binHigh);
    hMessy = new TH1D(hMname, hMname, binCount, binLow, binHigh);

    // Output tag
    TString outTag = "D1";

    plotAbsolute(CcutsTemp, McutsTemp, outTag);
}

// ---------------- Efficiency plotting ----------------
void plotAbsolute(TCut CcutsTemp, TCut McutsTemp, TString outTag) {
    // Fill histograms
    cleanTree->Draw(binVar + ">>" + hClean->GetName(), "ReWeight*" + TString(CcutsTemp));
    messyTree->Draw(binVar + ">>" + hMessy->GetName(), "ReWeight*" + TString(McutsTemp));

    // ---------------------------------------------------------
    // 1. Save Table to CSV
    // ---------------------------------------------------------
    TString csvName = "D1_occ/D1_Table_Combined_" + outTag + ".csv";
    std::ofstream csvFile(csvName.Data());

    if (csvFile.is_open()) {
        csvFile << "Bin Index,D1 Bin Center,hMessy Bin Content,hClean Bin Content,Ratio (Messy/Clean)\n";
        for (int i = 1; i <= binCount; i++) {
            double binCenter = hMessy->GetBinCenter(i);
            double contentMessy = hMessy->GetBinContent(i);
            double contentClean = hClean->GetBinContent(i);
            double ratio = (contentClean != 0) ? contentMessy / contentClean : 0.0;

            csvFile << i << "," << binCenter << "," << contentMessy << "," << contentClean << "," << ratio << "\n";
        }
        csvFile.close();
        std::cout << "CSV Table saved to: " << csvName << std::endl;
    } else {
        std::cerr << "Error writing CSV: " << csvName << std::endl;
    }

    // ---------------------------------------------------------
    // 2. Save Table to LaTeX
    // ---------------------------------------------------------
    TString texName = "D1_occ/D1_Table_Combined_" + outTag + ".tex";
    std::ofstream texFile(texName.Data());

    if (texFile.is_open()) {
        // Write LaTeX Preamble
        texFile << "\\documentclass{article}\n";
        texFile << "\\usepackage[utf8]{inputenc}\n";
        texFile << "\\usepackage{geometry}\n";
        texFile << "\\geometry{margin=1in}\n";
        texFile << "\\usepackage{booktabs}\n"; // For prettier tables
        texFile << "\\usepackage{longtable}\n"; // To handle tables splitting across pages
        texFile << "\\begin{document}\n";
        texFile << "\\begin{center}\n";
        
        // Begin Table
        texFile << "\\begin{longtable}{c c c c c}\n";
        texFile << "\\caption{Efficiency Data for " << outTag << "} \\\\\n";
        
        // Header
        texFile << "\\toprule\n";
        texFile << "Bin Index & Bin Center & Messy Content & Clean Content & Ratio (M/C) \\\\\n";
        texFile << "\\midrule\n";
        texFile << "\\endfirsthead\n";
        
        // Header for subsequent pages (if needed)
        texFile << "\\toprule\n";
        texFile << "Bin Index & Bin Center & Messy Content & Clean Content & Ratio (M/C) \\\\\n";
        texFile << "\\midrule\n";
        texFile << "\\endhead\n";
        
        // Footer for table end
        texFile << "\\bottomrule\n";
        texFile << "\\endfoot\n";

        // Loop over bins and write rows
        for (int i = 1; i <= binCount; i++) {
            double binCenter = hMessy->GetBinCenter(i);
            double contentMessy = hMessy->GetBinContent(i);
            double contentClean = hClean->GetBinContent(i);
            double ratio = (contentClean != 0) ? contentMessy / contentClean : 0.0;

            texFile << i << " & " 
                    << fixed << setprecision(1) << binCenter << " & " 
                    << setprecision(2) << contentMessy << " & " 
                    << setprecision(2) << contentClean << " & " 
                    << setprecision(4) << ratio << " \\\\\n";
        }

        texFile << "\\end{longtable}\n";
        texFile << "\\end{center}\n";
        texFile << "\\end{document}\n";

        texFile.close();
        std::cout << "LaTeX Table saved to: " << texName << std::endl;
    } else {
        std::cerr << "Error writing LaTeX: " << texName << std::endl;
    }
    // ---------------------------------------------------------

    getErrors(CcutsTemp, outTag);
}

void getErrors(TCut cutsTemp, TString outTag) {
    double weightSqClean = 0;
    double weightSumSqClean = 0;

    Double_t lowIntensityEdge, highIntensityEdge;
    TString intensityBin;

    TH1D* hWeightTemp = new TH1D("hWeightTemp", "hWeightTemp", 100, 0, 1e10);

    gr = new TGraphAsymmErrors(binCount);
    gr->SetMinimum(0);
    gr->SetMaximum(1.2);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->SetTitle(" ");

    for (int i = 1; i <= binCount; i++) {
        lowIntensityEdge = hClean->GetBinLowEdge(i);
        highIntensityEdge = lowIntensityEdge + hClean->GetBinWidth(i);

        intensityBin = TString::Format("D1 > %f && D1 < %f", lowIntensityEdge, highIntensityEdge);

        // Compute weights
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

    // --- Build interpolation curve ---
    int n = gr->GetN();
    grInterp = new TGraphAsymmErrors();
    for (int i = 0; i < n - 1; i++) {
        double x1,y1,x2,y2;
        gr->GetPoint(i, x1,y1);
        gr->GetPoint(i+1, x2,y2);

        double eyl1 = gr->GetErrorYlow(i);
        double eyh1 = gr->GetErrorYhigh(i);
        double eyl2 = gr->GetErrorYlow(i+1);
        double eyh2 = gr->GetErrorYhigh(i+1);

        double xm = 0.5*(x1+x2);
        double ym = 0.5*(y1+y2);

        double eym_low  = 0.5*sqrt(eyl1*eyl1 + eyl2*eyl2);
        double eym_high = 0.5*sqrt(eyh1*eyh1 + eyh2*eyh2);

        grInterp->SetPoint(i, xm, ym);
        grInterp->SetPointError(i, 0,0, eym_low, eym_high);
    }

    grInterp->SetMarkerColor(kRed);
    grInterp->SetMarkerStyle(24);
    grInterp->SetLineColor(kRed);
    grInterp->SetLineWidth(2);

    // --- Draw results ---
    TCanvas* c1 = new TCanvas();
    TString title = "D1_Efficiency_" + outTag;
    gr->SetTitle(title + ";D1;Efficiency");
    gr->Draw("AP");
    grInterp->Draw("LP SAME");

    c1->SaveAs("D1_occ/"+title+".pdf");

    // Save to file
    TFile* outFile = new TFile("D1_occ/"+title+".root","recreate");
    gr->Write("eff_original");
    grInterp->Write("eff_interpolated");
    outFile->Close();
}