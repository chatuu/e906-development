#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include "include/cuts/chuckcuts.h"

TString messyFileName = "/seaquest/users/chleung/pT_ReWeight/mc_jpsi_LD2_M027_S002_messy_occ_pTxFweight_v2.root";
//TString messyFileName = "/seaquest/users/chleung/pT_ReWeight/drellyan_LH2_messy_nhits.root";
TString cleanFileName = "/seaquest/users/chleung/pT_ReWeight/mc_jpsi_LD2_M027_S002_clean_occ_pTxFweight_v2.root";
//TString cleanFileName = "/seaquest/users/chleung/pT_ReWeight/drellyan_LH2_clean_nhits.root";
//TString messyFileName = "/data2/users/shivangi/work/analysis/R008/rootfilesMC/LH2_M026_S002_messy_v2_occ_mTables_pTxFweight1fl_v42.root";
//TString cleanFileName = "/data2/users/shivangi/work/analysis/R008/rootfilesMC/LH2_M026_S002_clean_v2_occ_mTables_pTxFweight1fl_v42.root";
//D2

TString binVar = "D2";
//TString binVar = "TriggerIntP";
double binLow = 0.0;
double binHigh = 400.0;
int binCount = 16;//16

/*
//Intensity
//TString binVar = "Intensity_p";
TString binVar = "(RFp00-32.73)*G2SEM/(QIEsum-369000*588*32.73)";
double binLow = 0.0;
double binHigh = 100000.0;
int binCount = 20;
*/
TTree *cleanTree;
TTree *messyTree;

void Efficiency(){
	/*TString cutsTemp = chuckCutsPositive && chuckCutsNegative && chuckCutsDimuon && physics && YC && occ;
	  TCut dataCuts = "sigWeight"*cutsTemp;
	  TString Cuts = "";
	  TCut messyCuts = "sigWeight"*Cuts;
	  */
	TString CcutsTemp = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && jPsiCut_MC ; 
	TString McutsTemp = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && occCuts_2111v42_Run56 && jPsiCut_MC; 

	TCut Ccuts = "ReWeight"*CcutsTemp;
	TCut Mcuts = "ReWeight"*McutsTemp;
	TFile* cleanFile = new TFile(cleanFileName, "READ");
	cleanTree = (TTree*)cleanFile->Get("Tree");

	TFile* messyFile = new TFile(messyFileName, "READ");
	messyTree = (TTree*)messyFile->Get("Tree");

	hClean = new TH1D("hClean","hClean",binCount,binLow,binHigh);
	hMessy = new TH1D("hMessy","hMessy",binCount,binLow,binHigh);
	/*	TCanvas *c1 = new TCanvas();
		cleanTree->Draw(binVar + ">>hClean", Ccuts);
		TCanvas *c2 = new TCanvas();
		messyTree->Draw(binVar + ">>hMessy", Mcuts);
		*/
	//	hClean->Sumw2();
	//	hMessy->Sumw2();
	//	TCanvas *c = new TCanvas();
	//        TH1D* hRatio = new TH1D("hRatio","hRatio",binCount,binLow,binHigh);
	//        hRatio->Add(hMessy);
	/*	TCanvas *c3 = new TCanvas();
		TH1D *hMessyCopy = (TH1D*)hMessy->Clone("hMessyCopy");
		hMessyCopy->Divide(hClean);
		hMessyCopy->Draw();
		*/
	/*	TGraphAsymmErrors* gr = new TGraphAsymmErrors();
		gr->BayesDivide(hMessy, hClean,"w");
		gr->Draw("AP");
		*/
	//	getErrors(cutsTemp);
	plotAbsolute(CcutsTemp, McutsTemp);
}
void plotAbsolute(TString CcutsTemp, TString McutsTemp)
{

	//cout<<cutsTemp<<endl;
	TCut Ccuts = "ReWeight"*CcutsTemp;
	TCut Mcuts = "ReWeight"*McutsTemp;
	TCanvas *c1 = new TCanvas();
	cleanTree->Draw(binVar + ">>hClean","ReWeight"*CcutsTemp);
	//      messyTree->Draw(binVar + ">>hMessy","gmcWeight"*cutsTemp);
	//        messyTree->Draw(binVar + ">>hMessy","sigWeight"*cutsTemp);
	TCanvas *c2 = new TCanvas();
	messyTree->Draw(binVar + ">>hMessy", "ReWeight"*McutsTemp);


	getErrors(CcutsTemp);
	//        fitRatio();
}

void getErrors(TString cutsTemp)
{

	//Cuts already includes weighting by gmcWeights
	//so divide the cuts by gmcWeight

	double weightSqClean;
	//      double weightSqMessy;
	double weightSqMessy;

	double weightSumSqClean;
	//      double weightSumSqMessy;
	double weightSumSqMessy;


	Double_t lowIntensityEdge;
	Double_t highIntensityEdge;


	TString intensityBin;

	TCanvas *c3 = new TCanvas();
	gStyle->SetStatY(0.9);                
	// Set y-position (fraction of pad size)
	gStyle->SetStatX(0.9);                
	// Set x-position (fraction of pad size)
	gStyle->SetStatW(0.2);                
	// Set width of stat-box (fraction of pad size)
	gStyle->SetStatH(0.2);                
	// Set height of stat-box (fraction of pad size)
	TH1D* hWeightTemp = new TH1D("hWeightTemp","hWeightTemp",100,0,10e9);
	gr = new TGraphAsymmErrors(binCount);
	gr->SetMinimum(0);
	gr->SetMaximum(1.2);
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->SetTitle(" ");
	//        gr->Fit("pol1","F EX0");

	for(int i = 1; i < binCount+1; i++)
	{


		lowIntensityEdge = hClean->GetBinLowEdge(i);
		highIntensityEdge = lowIntensityEdge + hClean->GetBinWidth(i);



		//intensityBin = TString::Format("Intensity_p > %f && Intensity_p < %f",lowIntensityEdge,highIntensityEdge);
		//intensityBin = TString::Format("(RFp00-32.73)*G2SEM/(QIEsum-369000*588*32.73) > %f && (RFp00-32.73)*G2SEM/(QIEsum-369000*588*32.73) < %f",lowIntensityEdge,highIntensityEdge);
		//intensityBin = TString::Format(binVar+" > %f && "+binVar+" < %f",lowIntensityEdge,highIntensityEdge);
		cout << intensityBin << endl;
		intensityBin = TString::Format("D2 > %f && D2 < %f",lowIntensityEdge,highIntensityEdge);


		cleanTree->Draw("ReWeight>>hWeightTemp",cutsTemp+"&&"+intensityBin);
		weightSumSqClean = hWeightTemp->GetMean()*hWeightTemp->GetEntries();
		weightSumSqClean = weightSumSqClean * weightSumSqClean;
		//      cout << weightSumSqClean << endl;
		cleanTree->Draw("ReWeight*ReWeight>>hWeightTemp",cutsTemp+"&&"+intensityBin);
		weightSqClean = hWeightTemp->GetMean()*hWeightTemp->GetEntries();
		//      cout << weightSqClean<< endl;

		/*
		   messyTree->Draw("gmcWeight>>hWeightTemp",cutsTemp);
		   weightSumSqMessy = hWeightTemp->GetMean()*hWeightTemp->GetEntries();
		   weightSumSqMessy = weightSumSqMessy * weightSumSqMessy;

		   messyTree->Draw("gmcWeight*gmcWeight>>hWeightTemp",cutsTemp);
		   weightSqMessy = hWeightTemp->GetMean()*hWeightTemp->GetEntries();

		   overlapTree->Draw("gmcWeight>>hWeightTemp",cutsTemp);
		   weightSumSqMessy = hWeightTemp->GetMean()*hWeightTemp->GetEntries();
		   weightSumSqMessy = weightSumSqMessy * weightSumSqMessy;

		   overlapTree->Draw("gmcWeight*gmcWeight>>hWeightTemp",cutsTemp);
		   weightSqOverlap = hWeightTemp->GetMean()*hWeightTemp->GetEntries();

*/


		Double_t x;
		Double_t y;
		Double_t exLow;
		Double_t exHigh;
		Double_t eyLow;
		Double_t eyHigh;

		Double_t midWilson;     //middle of the wilson score interval
		if(weightSumSqClean != 0)
		{
			Double_t w = weightSqClean/weightSumSqClean;    //x defined in BZ's note on weighted Wilson scores #8080
		} else
		{
			Double_t w = 0;
		}

		x = hMessy->GetBinCenter(i);
		exLow = 0;
		exHigh = 0;
		//                cout << x << endl;
		if(hClean->GetBinContent(i) != 0)
		{
			y = hMessy->GetBinContent(i)/hClean->GetBinContent(i);
		} else
		{
			y = 0;
		}
		//        if(y>1) continue;
		midWilson = (y + 0.5 * w)/(1 + w);
		//              cout << y << ";"<<midWilson << endl;
		if(y == 0 && w == 0)
		{
			eyHigh = 0;
			eyLow = 0;
		} else

		{
			eyHigh = midWilson + TMath::Sqrt(y*(1-y) * w + 0.25 * w*w)/(1 + w)  - y;
			eyLow = y - (midWilson - TMath::Sqrt(y*(1-y) * w + 0.25 * w*w)/(1 + w));
		}
		//         cout<< eyLow<< ":" << eyHigh << endl;

		gr->SetPoint(i-1,x,y);
		gr->SetPointError(i-1,exLow,exHigh,eyLow,eyHigh);
		//              cout<<"y: "<<y-eyLow<<" "<<y<<" "<<y+eyHigh<<endl;
		cout<<x<<" "<<y<<" "<<eyLow<<" "<<eyHigh<<endl;


	}
	TCanvas* c1 = new TCanvas();
	gr->Draw("AP");
	gr->SetTitle("Run 5-6 LD2;D2;Efficiency");
	TF1 *f1 = new TF1("f1","[0]+[1]*x",0,400);
	//TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x",0,400);
	gr->Fit("f1","R");
	gStyle->SetOptFit(1);
	c1->SaveAs("D2_occ/run5-6_jpsi_LD2.eps");
	TFile* outFile = new TFile("D2_occ/run5-6_jpsi_LD2.root","recreate");
	gr->Write("eff");
}

