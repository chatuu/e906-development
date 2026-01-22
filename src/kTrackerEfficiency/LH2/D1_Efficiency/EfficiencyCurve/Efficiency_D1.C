#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <iostream>
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



    // Loop over mass bins
//     for (int imass=0; imass<NBINS; imass++) {
//         TCut masscut = get_mass_cut(imass);
//         TString masslabel = get_mass_label(imass);

//         // Combined cuts (xF cut removed)
         TCut CcutsTemp = baseCcut;
         TCut McutsTemp = baseMcut;

//         // Histograms (naming updated)
         TString hCname = Form("hClean_mass_D1");
         TString hMname = Form("hMessy_mass_D1");
         hClean = new TH1D(hCname, hCname, binCount, binLow, binHigh);
         hMessy = new TH1D(hMname, hMname, binCount, binLow, binHigh);

//         // Output tag (naming updated)
         TString outTag = "D1";

         plotAbsolute(CcutsTemp, McutsTemp, outTag);
//     }
 }

// ---------------- Efficiency plotting ----------------
void plotAbsolute(TCut CcutsTemp, TCut McutsTemp, TString outTag) {
    cleanTree->Draw(binVar + ">>" + hClean->GetName(), "ReWeight*" + TString(CcutsTemp));
    messyTree->Draw(binVar + ">>" + hMessy->GetName(), "ReWeight*" + TString(McutsTemp));
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