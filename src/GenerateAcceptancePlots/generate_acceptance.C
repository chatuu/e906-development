/**
 * @file generate_acceptance.C
 * @brief A ROOT script to calculate and plot acceptance corrections for Drell-Yan analysis.
 * @details This script reads Monte Carlo events generated for LH2 and LD2 targets,
 * applies physics cuts, and generates acceptance histograms as a function of
 * dimuon mass for different bins of xF. It also calculates the ratio of
 * acceptances between the two targets.
 */

#include "include/cuts/chuckcuts.h"
#include "include/cuts/otherCuts.c"
#include <TMath.h>
#include "abs_cs.h"

/** @brief Array of colors for plotting histograms (kRed, kBlue, kBlack). */
int color[3]={kRed, kBlue, kBlack};

/**
 * @brief Calculates the binomial error.
 * @param n The total number of trials.
 * @param ns The number of successful trials.
 * @return The binomial error, calculated as sqrt(p*(1-p)/n). Returns 0 if n is 0.
 */
double biError(double n, double ns){
	if (n==0){
		return 0;
	}
	double p = ns/n;
	return TMath::Sqrt(p*(1-p)/n);
}

/**
 * @brief Calculates and sets the binomial error for a ratio histogram with weighted events.
 * @details This function iterates through each bin and calculates the error
 * for the ratio of accepted to thrown events, accounting for event weights.
 * The formula used is: error = (w_sqr / w_sum) * sqrt(p * (1 - p)), where p is the efficiency.
 * @param thrown A pointer to the histogram of thrown events. The bin content should be the sum of weights,
 * and the bin error should be the sum of squared weights.
 * @param accept A pointer to the histogram of accepted events (sum of weights).
 * @param ratio A pointer to the ratio histogram (accept/thrown) where the calculated bin errors will be set.
 */
void biErrorWeight(TH1F* thrown, TH1F* accept, TH1F* ratio){
	int nbins =  ratio->GetNbinsX();
	for (int i=0; i<2+nbins;i++){
		cout<<ratio->GetBinContent(i)<<endl;
		if (ratio->GetBinContent(i)==0){
			continue;
		}
		double w_sum=thrown->GetBinContent(i);
		double p=accept->GetBinContent(i)/w_sum;

		// Assumes the error in the thrown histogram is the sum of weights squared
		double w_sqr = thrown->GetBinError(i); 
		double error = w_sqr/w_sum*TMath::Sqrt(p*(1-p));

		ratio->SetBinError(i,error);
	}
}

/**
 * @brief Main function to generate and analyze acceptance correction histograms.
 * @details This function performs the following steps:
 * 1. Sets global ROOT styling options.
 * 2. Defines file paths for Monte Carlo (MC) simulation data for LH2 and LD2 targets.
 * 3. Opens the ROOT files and retrieves TTrees for "thrown" and "accepted" events.
 * 4. Defines the physics and MC cuts to be applied.
 * 5. Creates an output ROOT file to store resulting histograms.
 * 6. Loops over kinematic bins in xF (Feynman-x).
 * 7. For each xF bin:
 * a. Creates and fills histograms for thrown and accepted events vs. mass for LH2 and LD2.
 * b. Calculates the acceptance (accepted/thrown ratio) for each target.
 * c. Calculates the correct binomial errors for the weighted acceptance histograms using biErrorWeight().
 * d. Draws the LH2, LD2, and combined acceptance histograms.
 * e. Calculates and fits the ratio of LH2 to LD2 acceptance.
 * f. Writes all generated histograms to the output file.
 */
void generate_acceptance(){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	TString MCdir = "/data2/users/chleung/pT_ReWeight/";
	TString LH2_thrown_name = MCdir + "mc_drellyan_LH2_M027_S001_4pi_pTxFweight_v2.root";
	TString LH2_clean_name = MCdir + "mc_drellyan_LH2_M027_S001_clean_occ_pTxFweight_v2.root";
	TString LH2_messy_name = MCdir + "mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root";

	TFile* LH2_thrownFile = new TFile(LH2_thrown_name);
	TFile* LH2_acceptFile = new TFile(LH2_clean_name);

	TTree* LH2_thrownTree = (TTree*) LH2_thrownFile->Get("Tree");
	TTree* LH2_acceptTree = (TTree*) LH2_acceptFile->Get("Tree");

	TString LD2_thrown_name = MCdir + "mc_drellyan_LD2_M027_S001_4pi_pTxFweight_v2.root";
	TString LD2_clean_name = MCdir + "mc_drellyan_LD2_M027_S001_clean_occ_pTxFweight_v2.root";
	TString LD2_messy_name = MCdir + "mc_drellyan_LD2_M027_S001_messy_occ_pTxFweight_v2.root";

	TFile* LD2_thrownFile = new TFile(LD2_thrown_name);
	TFile* LD2_acceptFile = new TFile(LD2_clean_name);

	TTree* LD2_thrownTree = (TTree*) LD2_thrownFile->Get("Tree");
	TTree* LD2_acceptTree = (TTree*) LD2_acceptFile->Get("Tree");


	TCut allCut= chuckCutsPositive_2111v42 && chuckCutsNegative_2111v42 && chuckCutsDimuon_2111v42 &&  physicsCuts_2111v42;
	TCut MCcut= "1";
	//TCut MCcut= "mass > 4.5 && mass < 8.8 && xF < 0.95 && xF > -0.1";

	TFile * outFile =new TFile(TString::Format("acceptance_mass_xF.root"),"recreate");
	for (int i=0; i <nxF_CS; i++){
		TString title = TString::Format("%.2f<=xF<%.2f;mass",xFEdge[i],xFEdge[i+1]);
		TString binName = TString::Format("xF_bin%i",i);

		TCut binCut = TString::Format("%f<=xF && xF<%f",xFEdge[i], xFEdge[i+1]).Data();

		TCanvas* c1 = new TCanvas();
		TH1F* h_LH2_thrown = new TH1F("h_LH2_thrown_"+binName,title,nMass_CS,massEdge);
		TH1F* h_LH2_accept = new TH1F("h_LH2_accept_"+binName,title,nMass_CS,massEdge);

		TH1F* h_LD2_thrown = new TH1F("h_LD2_thrown_"+binName,title,nMass_CS,massEdge);
		TH1F* h_LD2_accept = new TH1F("h_LD2_accept_"+binName,title,nMass_CS,massEdge);

		h_LH2_thrown->Sumw2();
		h_LH2_accept->Sumw2();

		LH2_thrownTree->Draw("mass>>h_LH2_thrown_"+binName,"ReWeight"*(binCut&&MCcut), "goff");
		LH2_acceptTree->Draw("0.99*mass>>h_LH2_accept_"+binName,"ReWeight"*(allCut&&binCut), "goff");

		TH1F* h_ratio_LH2 = (TH1F*) h_LH2_accept->Clone("h_ratio_LH2_"+binName);
		h_ratio_LH2->Divide(h_LH2_thrown);
		biErrorWeight(h_LH2_thrown, h_LH2_accept, h_ratio_LH2);
		h_ratio_LH2->Draw();
		h_ratio_LH2->Write();

		h_LD2_thrown->Sumw2();
		h_LD2_accept->Sumw2();

		LD2_thrownTree->Draw("mass>>h_LD2_thrown_"+binName,"ReWeight"*(binCut&&MCcut), "goff");
		LD2_acceptTree->Draw("0.99*mass>>h_LD2_accept_"+binName,"ReWeight"*(allCut&&binCut), "goff");

		TH1F* h_ratio_LD2 = (TH1F*) h_LD2_accept->Clone("h_ratio_LD2_"+binName);
		h_ratio_LD2->Divide(h_LD2_thrown);
		biErrorWeight(h_LD2_thrown, h_LD2_accept, h_ratio_LD2);
		h_ratio_LD2->SetLineColor(kRed);
		h_ratio_LD2->Draw("same");
		h_ratio_LD2->Write();

		TH1F* h_ratio_combine = (TH1F*) h_ratio_LH2->Clone("h_ratio_combine_"+binName);
		h_ratio_combine->Add(h_ratio_LD2);
		h_ratio_combine->Scale(0.5);
		h_ratio_combine->SetLineColor(kBlack);
		h_ratio_combine->Draw("same");
		h_ratio_combine->Write();

		TLegend* leg = new TLegend(0.8,0.8,1,1);
		leg->AddEntry(h_ratio_LD2,"LD2","l");
		leg->AddEntry(h_ratio_LH2,"LH2","l");
		leg->AddEntry(h_ratio_combine,"LH2+LD2","l");
		leg->Draw();

		TCanvas* c2 = new TCanvas();

		TH1F* h_ratio_acceptance = (TH1F*) h_ratio_LH2->Clone("h_ratio_acceptance_"+binName);
		h_ratio_acceptance->Divide(h_ratio_LD2);
		h_ratio_acceptance->SetTitle("LH2/LD2 "+title);
		h_ratio_acceptance->Draw();
		h_ratio_acceptance->Fit("pol0");
		h_ratio_acceptance->Write();
	}
}

