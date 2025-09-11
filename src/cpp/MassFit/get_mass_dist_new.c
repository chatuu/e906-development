#include "include/cuts/chuckcuts.h"
#include "include/cuts/otherCuts.c"
#include "include/liveproton_v2.h"
#include "massfit.h"
#include "TFile.h"
#include "TTree.h"

void get_mass_dist_new(int rInt, int target){
	TString outName = "step11/"+ roads[rInt]+"_"+targetName[target]+"_mass.root";
	TFile* outFile = new TFile(outName, "recreate");

	int flaskPos = 2;
	if(target>=4){
		flaskPos=4;	
		component[4]="none";
	}
	TCut occCut ="1";
	if (rInt>=6){
		occCut = occCuts_2111v42_Run56;	
	}else{
		occCut = occCuts_2111v42;	
		
	}

	//TCut dataCut = physicsCuts_noMassCut_2111v42_tmp && occCut && intensityCuts_2111v42 && "mass<8.8";

	//Modified to apply on NMSU files:
	TCut dataCut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && occCut && "mass<8.8";

	//TCut dataCut = physicsCuts_noMassCut_2111v42_tmp && occCut && intensityCuts_2111v42_run56 && "mass<8.8"&&badRun6;
	TCut MCCut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && occCut && "0.99*mass<8.8";
	//TCut mixCut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && "mass<8.8" && TString::Format("targetPos==%i",target); 
	//TCut mixCut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && "mass<8.8 && (targetPos==1||targetPos==3)"; 

	//Modified to apply on NMSU files:
	TCut mixCut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_noMassCut_2111v42_tmp && occCut && "mass<8.8"; 
	//TCut mixCut = physicsCuts_noMassCut_2111v42_tmp  && "mass<8.8" && TString::Format("targetPos==%i",target); //FPGA1mix Run 2-3 

	//TFile* dataFile = new TFile(inName[rInt]);
	// Added this line to include NMSU RS-67 data file
	auto dataFile = new TFile("/seaquest/users/apun/e906_projects/rs67_merged_files/merged_RS67_3089LH2.root","OPEN");
	auto dataFlaskFile = new TFile("/seaquest/users/apun/e906_projects/rs67_merged_files/merged_RS67_3089flask.root","OPEN");

	// Added this line to include NMSU RS-67 data TTree
	auto dataTree = (TTree*) dataFile->Get("result");
	

	// Added this line to include NMSU RS-67 Mix TTree
	auto dataMixTree = (TTree*) dataFile->Get("result_mix");

	// Added this line to include NMSU RS-67 Flask TTree
	auto dataFlaskTree = (TTree*) dataFlaskFile->Get("result");

	// Added this line to include NMSU RS-67 Flask Mix TTree
	auto dataFlaskMixTree = (TTree*) dataFlaskFile->Get("result_mix");



	TH1D* hData = new TH1D("h_data", roads[rInt]+" "+targetName[target]+" data", nMass,massLo, massHi);
	TH1D* hFlask = new TH1D("h_flask", roads[rInt]+" "+component[4], nMass,massLo, massHi);
	TH1D* hFlaskMix = new TH1D("h_flask_mix", roads[rInt]+" "+component[4], nMass,massLo, massHi);
	
	hData->Sumw2();
	hFlask->Sumw2();
	hFlaskMix->Sumw2();
	
	//dataTree->Draw("mass>>h_data",dataCut&&TString::Format("targetPos==%i",target), "goff");
	dataTree->Draw("mass>>h_data",dataCut, "goff");
	//dataTree->Draw("mass>>h_flask",dataCut&&TString::Format("targetPos==%i",flaskPos),"goff");
	dataFlaskTree->Draw("mass>>h_flask",mixCut,"goff");
	dataFlaskMixTree->Draw("mass>>h_flask_mix",mixCut,"goff");
	
	outFile->cd();
	hData->Write();
	hFlask->Write();
	hFlaskMix->Write();
	dataFile->Close();

	for (int i=0; i<3 ;i++){
		TString nName = "n_"+component[i];
		TString hName = "h_"+component[i];
		TString wName = "w_"+component[i];
	
		cout<<MCname[target][i]<<endl;
		
		TFile* mcFile = new TFile(MCname[target][i]);
		auto mcTree = (TTree*) mcFile->Get("Tree");
		
		TH1D* hMC= new TH1D(nName, component[i], nMass, massLo, massHi);	
		TH1D* hMCw= new TH1D(hName, component[i], nMass, massLo, massHi);	
		
		hMC->Sumw2();
		hMCw->Sumw2();
		
		mcTree->Draw("0.99*mass>>"+nName, MCCut, "goff");
		mcTree->Draw("0.99*mass>>"+hName, "ReWeight"*MCCut, "goff");
		
		TH1D* hWeight = (TH1D*) hMCw->Clone(wName);
		hWeight->Divide(hMC);
		for(int j=1; j<nMass+1;j++){
			if (hWeight->GetBinContent(j)==0){
				hWeight->SetBinContent(j, 1)	;
			}	
		}
	
		outFile->cd();
		hMC->Write();
		hMCw->Write();
		hWeight->Write();
		mcFile->Close();
	}

	TFile* mixFile = new TFile(mixName);
	auto mixTree = (TTree*) mixFile->Get("Tree");
	
	TH1D* hmix = new TH1D("h_mix", "FPGA1 mix", nMass,massLo, massHi);
	//TH1D* hmix = new TH1D("h_mix", "liquidFPGA4 mix", nMass,massLo, massHi);
	
	hmix->Sumw2();
	
	mixTree->Draw("mass>>h_mix",mixCut, "goff");
	outFile->cd();
	hmix->Write();
	mixFile->Close();
	outFile->Close();
}

