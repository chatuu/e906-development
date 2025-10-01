double beamOffset = 1.6;

TCut chuckCutsPositive_2111v42 = TString::Format("chisq1_target < 15 && pz1_st1 > 9 && pz1_st1 < 75 && nHits1 > 13 && x1_t*x1_t + (y1_t-%f)*(y1_t-%f) < 320 && x1_d*x1_d + (y1_d-%f)*(y1_d-%f) < 1100 && x1_d*x1_d + (y1_d-%f)*(y1_d-%f) > 16 && chisq1_target < 1.5* chisq1_upstream && chisq1_target < 1.5*chisq1_dump && z1_v < -5 && z1_v > -320 && chisq1/(nHits1-5) < 12 && y1_st1/y1_st3 < 1 && abs( abs(px1_st1-px1_st3)-0.416) < 0.008 && abs(py1_st1-py1_st3) < 0.008 && abs(pz1_st1-pz1_st3) < 0.08 && y1_st1*y1_st3 > 0 && abs(py1_st1)>0.02", beamOffset,beamOffset,beamOffset,beamOffset,beamOffset,beamOffset).Data();

TCut chuckCutsNegative_2111v42 = TString::Format("chisq2_target < 15 && pz2_st1 > 9 && pz2_st1 < 75 && nHits2 > 13 && x2_t*x2_t + (y2_t-%f)*(y2_t-%f) < 320 && x2_d*x2_d + (y2_d-%f)*(y2_d-%f) < 1100 && x2_d*x2_d + (y2_d-%f)*(y2_d-%f) > 16 && chisq2_target < 1.5* chisq2_upstream && chisq2_target < 1.5*chisq2_dump && z2_v < -5 && z2_v > -320 && chisq2/(nHits2-5) < 12 && y2_st1/y2_st3 < 1 && abs(abs(px2_st1-px2_st3)-0.416) < 0.008 && abs(py2_st1-py2_st3) < 0.008 && abs(pz2_st1-pz2_st3) < 0.08 && y2_st1*y2_st3 > 0 && abs(py2_st1)>0.02",beamOffset,beamOffset, beamOffset,beamOffset,beamOffset,beamOffset).Data();

TCut chuckCutsDimuon_2111v42 = TString::Format("abs(dx) < 0.25 && abs(dy-%f)< 0.22 && dz > -280 && dz < -5 && abs(dpx) < 1.8 && abs(dpy) < 2 && dpx*dpx + dpy*  dpy < 5 && dpz > 38 && dpz < 116 && dx*dx + (dy-%f)*(dy-%f) < 0.06 && abs(trackSeparation) < 270 && chisq_dimuon < 18 && abs(chisq1_target+ chisq2_target-chisq_dimuon) < 2 && y1_st3*y2_st3 < 0 && nHits1 + nHits2 > 29 &&  nHits1St1 + nHits2St1 > 8 && abs(x1_st1+x2_st1)<42 ",beamOffset,beamOffset,beamOffset).Data();      

TCut physicsCuts_2111v42 = "mass > 4.2 && mass < 8.8 && xF < 0.95 && xF > -0.1 && xT > 0.05 && xT < 0.55 && abs(costh) < 0.5";                                  
TCut physicsCuts_noxT_2111v42 = "mass > 4.2 && mass < 8.8 && xF < 0.95 && xF > -0.1 && abs(costh) < 0.5";                                  
TCut physicsCuts_noMassCut_2111v42 = " xF < 0.95 && xF > -0.1 && xT > 0.05 && xT < 0.55 && abs(costh) < 0.5";                                  

TCut occCuts_2111v42 = "D1 < 400 && D2 < 400 && D3 < 400 && D1+D2+D3<1000";
TCut occCuts_2111v42_Run56 = "D0 < 400 && D2 < 400 && D3 < 400 && D0+D2+D3<1000";

//TCut intensityCuts_2111v42 = "(RFp00-34)*PotPerQie >0 && (RFp00-34)*PotPerQie<80000";
TCut intensityCuts_2111v42 = "(RFp00-34)*G2SEM/(QIEsum-369000*588*34) >0 && (RFp00-34)*G2SEM/(QIEsum-369000*588*34)<80000";
TCut intensityCuts_2111v42_run56 = "tInt >0 && tInt<80000";

TCut jPsiCut ="mass>2.7 && mass<3.5 && xF>0 && xF<0.95 && xT>0.05 && xT<0.15";
TCut psipCut ="mass>3.3 && mass<4.1 && xF>0 && xF<0.95 && xT>0.05 && xT<0.2";
TCut jPsiCut_MC ="0.99*mass>2.7 && 0.99*mass<3.5 && xF>0 && xF<0.95 && xT>0.05 && xT<0.15";
//TCut jPsiCut ="mass>2.8 && mass<3.4 && xF>0 && xF<0.95 && xT>0.05 && xT<0.15";
TCut DYCut ="mass>4.5 && mass<8.8 &&xB<0.8";
TCut DYCut_MC ="0.99*mass>4.5 && 0.99*mass<8.8 &&xB<0.8";

TCut chuckCutsPositive_2111v42_tmp = TString::Format("chisq1_target < 15 && pz1_st1 > 9 && pz1_st1 < 75 && nHits1 > 13 && x1_t*x1_t + (y1_t-%f)*(y1_t-%f) < 320 && x1_d*x1_d + (y1_d-%f)*(y1_d-%f) < 1100 && x1_d*x1_d + (y1_d-%f)*(y1_d-%f) > 8 && chisq1_target < 1.5* chisq1_upstream && chisq1_target < 1.5*chisq1_dump && z1_v < -5 && z1_v > -320 && chisq1/(nHits1-5) < 12 && y1_st1/y1_st3 < 1 && abs( abs(px1_st1-px1_st3)-0.416) < 0.008 && abs(py1_st1-py1_st3) < 0.008 && abs(pz1_st1-pz1_st3) < 0.08 && y1_st1*y1_st3 > 0 && abs(py1_st1)>0.02", beamOffset,beamOffset,beamOffset,beamOffset,beamOffset,beamOffset).Data();

TCut chuckCutsNegative_2111v42_tmp = TString::Format("chisq2_target < 15 && pz2_st1 > 9 && pz2_st1 < 75 && nHits2 > 13 && x2_t*x2_t + (y2_t-%f)*(y2_t-%f) < 320 && x2_d*x2_d + (y2_d-%f)*(y2_d-%f) < 1100 && x2_d*x2_d + (y2_d-%f)*(y2_d-%f) > 8 && chisq2_target < 1.5* chisq2_upstream && chisq2_target < 1.5*chisq2_dump && z2_v < -5 && z2_v > -320 && chisq2/(nHits2-5) < 12 && y2_st1/y2_st3 < 1 && abs(abs(px2_st1-px2_st3)-0.416) < 0.008 && abs(py2_st1-py2_st3) < 0.008 && abs(pz2_st1-pz2_st3) < 0.08 && y2_st1*y2_st3 > 0 && abs(py2_st1)>0.02",beamOffset,beamOffset, beamOffset,beamOffset,beamOffset,beamOffset).Data();

TCut physicsCuts_noMassCut_2111v42_tmp = "  xF < 0.95 && xF > -0.1  && xT>0.005 && xT < 0.55 && abs(costh) < 0.5";                                  

TCut chuckCutsTrack_2111v42_tmp = "(((x1_t*x1_t + (y1_t-0.400000)*(y1_t-0.400000) < 320 && x1_d*x1_d + (y1_d-0.400000)*(y1_d-0.400000) < 1100 && x1_d*x1_d + (y1_d-0.400000)*(y1_d-0.400000) > 16 && runID < 11000) || (x1_t*x1_t + (y1_t-1.600000)*(y1_t-1.600000) < 320 && x1_d*x1_d + (y1_d-1.600000)*(y1_d-1.600000) < 1100 && x1_d*x1_d + (y1_d-1.600000)*(y1_d-1.600000) > 16 && runID > 11000))) && (((x2_t*x2_t + (y2_t-0.400000)*(y2_t-0.400000) < 320 && x2_d*x2_d + (y2_d-0.400000)*(y2_d-0.400000) < 1100 && x2_d*x2_d + (y2_d-0.400000)*(y2_d-0.400000) > 16 && runID < 11000) || (x2_t*x2_t + (y2_t-1.600000)*(y2_t-1.600000) < 320 && x2_d*x2_d + (y2_d-1.600000)*(y2_d-1.600000) < 1100 && x2_d*x2_d + (y2_d-1.600000)*(y2_d-1.600000) > 16 && runID > 11000)))";

