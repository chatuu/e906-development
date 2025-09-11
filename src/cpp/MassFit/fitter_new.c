#include "include/cuts/chuckcuts.h"
#include "include/liveproton_v2.h"
#include "massfit.h"

#include <iostream>
#include <fstream>
#include <vector>
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"

int color[] = {kPink, kViolet+2, kOrange-3, kSpring-6, kBlack};


void fitter_new(int rInt, int target){
	set_PoT();
	gStyle->SetOptStat(0);
	TString s1Name = "step1/"+ roads[rInt]+"_"+targetName[target]+"_mass.root";
	TString s11Name = "step11/"+ roads[3]+"_"+targetName[1]+"_mass.root";
	TString outName =	"step2/"+ roads[rInt]+"_"+targetName[target]+"";
	TString fitTitle = roads[rInt]+" "+ targetName[target]+" mass fit;mass (GeV);count/(0.20GeV)";
	TString pic = ".pdf";
	int flaskPos = 2;
	if(target>=4){
		flaskPos=4;
		component[4]="none";
	}
	Double_t sFlask = liveproton[rInt][target]/liveproton[rInt][flaskPos];

	TFile* s1File = new TFile(s1Name);
	TFile* s11File = new TFile(s11Name);
	TH1D* hData = (TH1D*) s1File->Get("h_data");
	TH1D* hMC[3];
	TH1D* wMC[3];
	TH1D* hMCw[3];
	for(int i=0; i<3;i++){
		hMC[i]=(TH1D*) s1File->Get("n_"+component[i]);
		wMC[i]=(TH1D*) s1File->Get("w_"+component[i]);
		hMCw[i]=(TH1D*) s1File->Get("h_"+component[i]);
	}
	TH1D* hMix = (TH1D*) s1File->Get("h_mix");
	//TH1D* hMix = (TH1D*) s11File->Get("h_mix");
	TH1D* hFlask = (TH1D*) s1File->Get("h_"+component[4]);

	Double_t eData=0;
	Double_t nData = hData->IntegralAndError(1,nMass, eData);
	Double_t n0[5], e0[5];
	for(int i=0; i<3; i++){
		n0[i] = hMCw[i]->IntegralAndError(1,nMass,e0[i]);
	}
	n0[3] = hMix->IntegralAndError(1,nMass,e0[3]);
	n0[4] = hFlask->IntegralAndError(1,nMass,e0[4]);

	TObjArray *mc = new TObjArray(5);
	for(int i=0; i<3; i++){
		mc->Add(hMC[i]);
	}
	mc->Add(hMix);
	mc->Add(hFlask);

	TFractionFitter* fit = new TFractionFitter(hData,mc);

	for (int i=0; i<3; i++){
		fit->SetWeight(i, wMC[i]);
	}
	double fFlask = sFlask * n0[4]/nData;
	//double fFlask = 4.96843 * n0[4]/nData;
	double uncertainty = 0.001;
	fit->Constrain(4, fFlask * (1.0 - uncertainty), fFlask * (1.0 + uncertainty));
	//fit->Constrain(4, fFlask, fFlask);
	//fit->Constrain(4, fFlask, fFlask);

	TCanvas* c1 = new TCanvas();
	TLegend* leg = new TLegend(0.8,0.7,1,1);

	TFitResultPtr fit_result = fit->Fit();
	if (!fit_result->IsValid()) {
		std::cout << "Fit failed!" << std::endl;
		return;
	}

	//std::cout << "Fit Result Chi2/NDF: " <<  fit_result->Chisquare() << "/" << fit_result->Ndf() << std::endl;
	std::cout << "Fit Result Chi2/NDF: " << fit->GetChisquare() << "/" << fit->GetNDF() << std::endl;


	double f[5], ef[5];
	double s0[5], es0[5];
	double n1[5], en1[5];

	TFile* outGr = new TFile(outName+".root", "recreate");
	TH1D* mc_new[5];

	hData->SetMarkerStyle(21);
	hData->SetTitle(fitTitle);
	hData->DrawCopy("Ep");
	outGr->cd();
	hData->Write("data");
	leg->AddEntry(hData, "data","p");

	TH1D* result = (TH1D*)fit->GetPlot();
	result->SetLineWidth(3);
	leg->AddEntry(result, "Total", "l");
	result->Draw("same");
	TPaveText *pt = new TPaveText(.8,.4,1,.7,"nbNDC");
	pt->SetBorderSize(1);
	pt->AddText(TString::Format( "Total = %.0f", nData));

	ofstream outFile;
	outFile.open(outName+".csv");
	outFile<<"comp/C:scale/D:eScale/D"<<endl;
	for (int i=0; i<5; i++){
		fit->GetResult(i,f[i],ef[i]);
		n1[i]=f[i]*nData;
		en1[i]=ef[i]*nData;
		s0[i]=n1[i]/n0[i];
		es0[i]=en1[i]/n0[i];
		outFile<<component[i]<<", "<<s0[i]<<", "<<es0[i]<<endl;
		pt->AddText( component[i]+TString::Format( " = %.0f+/-%.0f", n1[i],en1[i]) );
		mc_new[i] = (TH1D*)fit->GetMCPrediction(i)->Clone();
		if(i<3){
			mc_new[i]->Multiply(wMC[i]);
		}
		double etmp;
		double tmp=mc_new[i]->IntegralAndError(1,nMass,etmp);
		if (tmp > 0) mc_new[i]->Scale(n1[i]/tmp);
		mc_new[i]->SetLineColor(color[i]);
		mc_new[i]->SetLineWidth(2);
		mc_new[i]->SetLineStyle(9);
		mc_new[i]->Draw("hist same");
		leg->AddEntry(mc_new[i],component[i],"l");
		outGr->cd();
		mc_new[i]->Write(component[i]);
	}

	const TMatrixDSym& cov_mat = fit_result->GetCovarianceMatrix();
	int n_par = fit_result->NPar();

	std::vector<int> free_params;
	for(int i = 0; i < n_par; ++i) {
		if (!fit_result->IsParameterFixed(i)) {
			free_params.push_back(i);
		}
	}
	int n_free_par = free_params.size();
	TMatrixD cor(n_free_par, n_free_par);
	
    // --- START: Manual Chi2 and NDF Calculation ---
    
    std::cout << "--- Manual Calculation ---" << std::endl;

    // 1. Calculate NDF
    int n_bins_used = 0;
    for (int i = 1; i <= hData->GetNbinsX(); ++i) {
        // Only count bins that have data to be included in the fit
        if (hData->GetBinContent(i) > 0) {
            n_bins_used++;
        }
    }
    int manual_NDF = n_bins_used - n_free_par;
    std::cout << "Number of bins with data: " << n_bins_used << std::endl;
    std::cout << "Number of free parameters: " << n_free_par << std::endl;
    std::cout << "Manual NDF = " << manual_NDF << std::endl;

    // 2. Calculate Pearson Chi2
    double manual_chi2 = 0;
    // 'result' is the histogram with the total fit prediction
    for (int i = 1; i <= hData->GetNbinsX(); ++i) {
        double data_val = hData->GetBinContent(i);
        // Only include bins with data in the chi2 calculation
        if (data_val > 0) {
            double fit_val = result->GetBinContent(i);
            double error_sq = hData->GetBinError(i) * hData->GetBinError(i);
            if (error_sq > 0) {
                manual_chi2 += TMath::Power(data_val - fit_val, 2) / error_sq;
            }
        }
    }
    std::cout << "Manual Chi2 = " << manual_chi2 << std::endl;
    std::cout << "Manual Chi2/NDF = " << (manual_NDF > 0 ? manual_chi2 / manual_NDF : 0) << std::endl;
    std::cout << "-------------------------" << std::endl;
    
    // --- END: Manual Chi2 and NDF Calculation ---

	outGr->cd();
	cor.Write("cor_mat"); // This was missing in the original script
	outFile.close();
    
    // Use the MANUALLY calculated values on the plot
	pt->AddText(TString::Format( "#chi^{2}/NDF = %.1f/%d", manual_chi2, manual_NDF));
	leg->Draw();
	pt->Draw();
	c1->SetTickx(1);
	c1->SetTicky(1);
	c1->SaveAs(outName+pic);
	c1->SetLogy();
	c1->SaveAs(outName+"_log"+pic);
}
