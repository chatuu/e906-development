#include <iostream>
#include <string.h>
#include "Riostream.h"
#include "TEfficiency.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TArrow.h"
#include "TBox.h"
#include "TVirtualPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TVector.h"
#include "TVectorD.h"
#include "TClass.h"
#include "Math/QuantFuncMathCore.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "./dataAndMC_rootfiles.h"
#include "./cuts.h"
#include "./get_xf_cut.cxx"
#include "./diffCross_acc2d_H.cxx"
//#include <diffCross_acc2d_H_syst.cxx>
void go_diffCross_H(){
int i=15;//0-15
const Int_t NBINS = 11;
Double_t edges[NBINS + 1] = {4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7};
diffCross_acc2d_H(i);
/*
for(i=0;i<16;i++){
cout<<i<<endl;
diffCross_acc2d_H(i);
}
*/
//diffCross_acc2d_H_syst(i, true, false, false, false);
//diffCross_acc2d_H_syst(i, false, true, false, false);
//diffCross_acc2d_H_syst(i, false, false, true, true);
}
