//#include "cuts.h"
#include "./norm_massfit.h"
//#include "dataAndMC_rootfiles.h"
//#include "get_xf_cut.cxx"
void diffCross_acc2d_H(int binNum){
Double_t kEff[16]={1.77909,1.74145,1.62481,1.63713,1.57874,1.58371,1.47446,1.4607,1.42335,1.43693,1.36176,1.35986,1.37244,1.34141,1.33321,1.30802};	
const Int_t NBINS = 11;
Double_t edges[NBINS + 1] = {4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7};
	TString binNumStr = TString::Format("%d",binNum);
	TCut xFCut = get_xf_cut(binNum);
	TCut gmcTempxF = gmcTemp && xFCut;

	TFile *saveAcc = new TFile("acceptance_h.root","READ");
	TH1F *h2 = (TH1F*)saveAcc->Get("LH2_acc_M_xF_"+binNumStr);
	
	TCanvas *c = new TCanvas();
	c->SetLogy();
	gStyle->SetOptStat(0);
	TFile* dataFile = new TFile(rootFile, "READ");
	TTree* dataTree = (TTree*)dataFile->Get("Tree");
	TFile* mixFile = new TFile(mixFileName, "READ");
	TTree* mixTree = (TTree*)mixFile->Get("Tree");
	TFile* jpFile = new TFile(jpsiFileName, "READ");
	TTree* jpTree = (TTree*)jpFile->Get("Tree");
	TFile* ppFile = new TFile(psipFileName, "READ");
	TTree* ppTree = (TTree*)ppFile->Get("Tree");

	TH1F *hData= new TH1F("hData","hData",NBINS, edges);
	hData->Sumw2();
	TH1F *hDataF= new TH1F("hDataF","hDataF",NBINS, edges);//0.3 bin size 18 bins, 0.5 bin size 11 bins, 0.1 bin size 55 bins
	hDataF->Sumw2();
	TH1F *hDataMix= new TH1F("hDataMix","hDataMix",NBINS, edges);
	hDataMix->Sumw2();
	TH1F *hMC_jp= new TH1F("hMC_jp","hMC_jp",NBINS, edges);
	hMC_jp->Sumw2();
	TH1F *hMC_pp= new TH1F("hMC_pp","hMC_pp",NBINS, edges);
	hMC_pp->Sumw2();
	
	dataTree->Draw("mass>>hData","mass^3"*(gmcTemp_charm && physics && target1 && xFCut && occ && tInt));
	//dataTree->Draw("mass>>hData","mass^3*(abs(costh)<0.5 && xF>=0.45 && xF < 0.5 && xT>0.05 && xT<.55 && mass>4.2 && mass<8.8 && targetPos==3)");
	dataTree->Draw("mass>>hDataF","mass^3"*(gmcTemp_charm && physics && xFCut && flask && occ && tInt));
	//dataTree->Draw("mass>>hDataF","mass^3*(abs(costh)<0.5 && xF>=0.45 && xF < 0.5 && xT>0.05 && xT<.55 && mass>4.2 && mass<8.8 && targetPos==2)");
	mixTree->Draw("mass>>hDataMix","mass^3"*(gmcTemp_charm && physics && xFCut && (target1 || target3)));
	jpTree->Draw("0.99*mass>>hMC_jp","0.99^3*mass^3*sigWeight"*(gmcTemp_charm && mass99 && xFCut));
	ppTree->Draw("0.99*mass>>hMC_pp","0.99^3*mass^3*sigWeight"*(gmcTemp_charm && mass99 && xFCut));
	// scaling and normalization from mass fits in norm_massfit.h
	for(int k=0;k<NBINS;k++){
	cout<<hData->GetBinContent(k+1)<<","<<hDataF->GetBinContent(k+1)<<","<<hDataMix->GetBinContent(k+1)<<","<<hMC_jp->GetBinContent(k+1)<<","<<hMC_pp->GetBinContent(k+1)<<","<<h2->GetBinContent(k+1)<<endl;
	}	
	hDataF->Scale(flaskNorm); //flasknorm is ratio of liveP_targ/liveP_fl
	hData->Add(hDataF,-1); // subtract flask contribution
	hDataMix->Scale(mixNorm); // mix background normalization from massfit
	hData->Add(hDataMix,-1);  // subtract mix bkg. contribution
	hMC_jp->Scale(jpNorm); //scale jpsi by normalization from massfit
	hMC_pp->Scale(ppNorm); //scale psi' similarly
	hData->Add(hMC_jp,-1); // subtract jpsi contribution
	hData->Add(hMC_pp,-1); // subtract psi' contribution
	hData->Divide(h2);
	hData->Scale(kEff[binNum]); //scale by 1/kEfficiency
	cout<<kEff[binNum]<<endl;
	//hData->Scale(3.35607e-8);// luminosity (includes delta xF but not  delta M- delta M is not uniform)
	hData->Scale(3.48489e-8);// updated live PoT from Hugo Doc 10679 ;luminosity (includes delta xF but not delta M- delta M is not uniform)
	hData->Scale(1.102); // hodo eff
	for(int k=0;k<NBINS;k++){
	hData->SetBinContent(k+1,hData->GetBinContent(k+1)/(edges[k+1]-edges[k]));
	hData->SetBinError(k+1,hData->GetBinError(k+1)/(edges[k+1]-edges[k]));
	}
	for(int k=0;k<NBINS;k++){
	cout<<hData->GetBinContent(k+1)<<endl;
	}
	for(int k=0;k<NBINS;k++){
	if(hData->GetBinError(k+1)/hData->GetBinContent(k+1)>0.99){
		hData->SetBinContent(k+1,0.0);
		hData->SetBinError(k+1,0.0);
			continue;
		}
	}

	hData->SetBinContent(1,0.0);
	hData->Draw();
	hData->GetXaxis()->SetRangeUser(4.2,8.5);
	//hData->GetYaxis()->SetRangeUser(.001,3);
	hData->SetMarkerColor(kRed);
	hData->SetMarkerStyle(20);
	hData->SetMarkerSize(1);
//	hData->SetLineColor(kRed);

	hData->GetXaxis()->SetTitle("M (GeV)");
	hData->GetYaxis()->SetTitle("M^{3} d^{2}#sigma/dMdx_{F} (nb GeV^{2})");
	hData->SetTitle("");
	TFile *nnpdf4_file = new TFile("/seaquest/users/chleung/theory/new/v2/E906/NNPDF40_xFnew_p.root","READ");
	TGraphAsymmErrors *nnpdf4 = (TGraphAsymmErrors*)nnpdf4_file->Get("gr_xFbin"+binNumStr);//"G_MaxMin");
	nnpdf4->SetMarkerSize(0.75);
	nnpdf4->SetMarkerStyle(21);
	nnpdf4->SetMarkerColor(kBlue+2);
	nnpdf4->SetLineColor(kBlue+2);
	nnpdf4->SetFillColorAlpha(38, 0.5);
	nnpdf4->SetFillStyle(3002);
	nnpdf4->Draw("L3");
	TFile *ct18_file = new TFile("/seaquest/users/chleung/theory/new/v2/E906/CT18_xFnew_p.root","READ");
	TGraphAsymmErrors *ct18 = (TGraphAsymmErrors*)ct18_file->Get("gr_xFbin"+binNumStr);//"G_MaxMin");
	ct18->SetMarkerSize(0.75);
	ct18->SetMarkerStyle(21);
	ct18->SetMarkerColor(kGreen+2);
	ct18->SetLineColor(kGreen+2);
	ct18->SetFillColorAlpha(30, 0.5);
	ct18->SetFillStyle(3002);
	ct18->Draw("L3");
	hData->Draw("same");
	auto legend = new TLegend(0.6,0.7,0.85,0.85);
	legend->SetBorderSize(0);
	legend->AddEntry(hData,"E906","lep");
	legend->AddEntry(nnpdf4,"NNPDF 4.0","lf");
	legend->AddEntry(ct18,"CT18","lf");
	legend->Draw();
//	c->SaveAs("plots/LH2_"+binNumStr+"_roofit.eps");
	TH1F *hAccCor = (TH1F*)hData->Clone("hAccCor");
//	hAccCor->Scale(1.1183e-7);
//	hAccCor->Draw();
//
	
	TFile *savehist = new TFile("result_rootfiles/LH2_"+binNumStr+"_updatedPoT.root","RECREATE");
	hAccCor->Write();
	savehist->Close();
	delete savehist;
	
}
