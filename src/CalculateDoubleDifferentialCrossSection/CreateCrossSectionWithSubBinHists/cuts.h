#ifndef CUTS_H
#define CUTS_H

// --- Add necessary ROOT headers ---
#include "TCut.h"
#include "TString.h"
#include "Rtypes.h"

// Define constants and cuts
Double_t beamOffset = 1.604; // in cm, for 2011v42

// Function to generate the xF cut for a specific bin
TCut get_xf_cut(int binNumber) {
    const Double_t xf_bins[] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8
    };
    if (binNumber < 0 || binNumber >= 16) {
        return ""; // Return an empty cut for invalid bin numbers
    }
    // --- FIX: Added .Data() to convert TString to const char* ---
    return TString::Format("xF >= %f && xF < %f", xf_bins[binNumber], xf_bins[binNumber+1]).Data();
}

// Define the individual cuts
TCut chuckCutsPositive_2111v42_tmp = TString::Format("chisq1_target < 15 && pz1_st1 > 9 && pz1_st1 < 75 && nHits1 > 13 && x1_t*x1_t + (y1_t-%f)*(y1_t-%f) < 320 && x1_d*x1_d + (y1_d-%f)*(y1_d-%f) < 1100 && x1_d*x1_d + (y1_d-%f)*(y1_d-%f) > 8 && chisq1_target < 1.5* chisq1_upstream && chisq1_target < 1.5*chisq1_dump && z1_v < -5 && z1_v > -320 && chisq1/(nHits1-5) < 12 && y1_st1/y1_st3 < 1 && abs( abs(px1_st1-px1_st3)-0.416) < 0.008 && abs(py1_st1-py1_st3) < 0.008 && abs(pz1_st1-pz1_st3) < 0.08 && y1_st1*y1_st3 > 0 && abs(py1_st1)>0.02", beamOffset,beamOffset,beamOffset,beamOffset,beamOffset,beamOffset).Data();
TCut chuckCutsNegative_2111v42_tmp = TString::Format("chisq2_target < 15 && pz2_st1 > 9 && pz2_st1 < 75 && nHits2 > 13 && x2_t*x2_t + (y2_t-%f)*(y2_t-%f) < 320 && x2_d*x2_d + (y2_d-%f)*(y2_d-%f) < 1100 && x2_d*x2_d + (y2_d-%f)*(y2_d-%f) > 8 && chisq2_target < 1.5* chisq2_upstream && chisq2_target < 1.5*chisq2_dump && z2_v < -5 && z2_v > -320 && chisq2/(nHits2-5) < 12 && y2_st1/y2_st3 < 1 && abs(abs(px2_st1-px2_st3)-0.416) < 0.008 && abs(py2_st1-py2_st3) < 0.008 && abs(pz2_st1-pz2_st3) < 0.08 && y2_st1*y2_st3 > 0 && abs(py2_st1)>0.02",beamOffset,beamOffset, beamOffset,beamOffset,beamOffset,beamOffset).Data();
TCut chuckCutsDimuon_2111v42 = TString::Format("abs(dx) < 0.25 && abs(dy-%f)< 0.22 && dz > -280 && dz < -5 && abs(dpx) < 1.8 && abs(dpy) < 2 && dpx*dpx + dpy*dpy < 5 && dpz > 38 && dpz < 116 && dx*dx + (dy-%f)*(dy-%f) < 0.06 && abs(trackSeparation) < 270 && chisq_dimuon < 18 && abs(chisq1_target+ chisq2_target-chisq_dimuon) < 2 && y1_st3*y2_st3 < 0 && nHits1 + nHits2 > 29 && nHits1St1 + nHits2St1 > 8 && abs(x1_st1+x2_st1)<42 ",beamOffset,beamOffset,beamOffset).Data();
TCut physicsCuts_2111v42 = "mass > 4.2 && mass < 8.8 && xF < 0.95 && xF > -0.1 && xT > 0.05 && xT < 0.55 && abs(costh) < 0.5";
TCut physicsCutsNoMass_2111v42 = "mass > 4.2 && mass < 8.8 && xF < 0.95 && xF > -0.1 && xT > 0.05 && xT < 0.55 && abs(costh) < 0.5";
TCut occCuts_2111v42 = "D1 < 400 && D2 < 400 && D3 < 400 && D1+D2+D3<1000";

// Final combined data cut used in the analysis
TCut dataCut = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCuts_2111v42 && occCuts_2111v42 && "mass<8.8";
TCut dataCutNoMass = chuckCutsPositive_2111v42_tmp && chuckCutsNegative_2111v42_tmp && chuckCutsDimuon_2111v42 && physicsCutsNoMass_2111v42 && occCuts_2111v42 && "mass<8.8";

#endif // CUTS_H