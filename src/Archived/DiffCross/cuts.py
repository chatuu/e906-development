# cuts.py
# v42 cuts

chuckCutsPositive = "chisq1_target<15 && pz1_st1 > 9 && pz1_st1 < 75 && nHits1 > 13 && ((x1_t*x1_t + (y1_t-1.6)*(y1_t-1.6) < 320 && x1_d*x1_d + (y1_d-1.6)*(y1_d-1.6) > 16 && x1_d*x1_d + (y1_d-1.6)*(y1_d-1.6) < 1100 && runID> 11000) || (x1_t*x1_t + (y1_t-0.4)*(y1_t-0.4) < 320 && x1_d*x1_d + (y1_d-0.4)*(y1_d-0.4) > 16 && x1_d*x1_d + (y1_d-0.4)*(y1_d-0.4) < 1100 && runID < 11000)) && chisq1_target < 1.5*chisq1_upstream && chisq1_target < 1.5*chisq1_dump && z1_v < -5 && z1_v > -320 && chisq1/(nHits1-5) < 12 && abs(abs(px1_st1-px1_st3) - 0.416) < 0.008 && abs(py1_st1-py1_st3) < 0.008 && abs(pz1_st1-pz1_st3) < 0.08 && y1_st1/y1_st3<1 && y1_st1*y1_st3>0 && abs(py1_st1)>0.02"

chuckCutsNegative = "chisq2_target<15 && pz2_st1 > 9 && pz2_st1 < 75 && nHits2 > 13 && ((x2_t*x2_t + (y2_t-1.6)*(y2_t-1.6) < 320 && x2_d*x2_d + (y2_d-1.6)*(y2_d-1.6) > 16 && x2_d*x2_d + (y2_d-1.6)*(y2_d-1.6) < 1100 && runID> 11000) || (x2_t*x2_t + (y2_t-0.4)*(y2_t-0.4) < 320 && x2_d*x2_d + (y2_d-0.4)*(y2_d-0.4) > 16 && x2_d*x2_d + (y2_d-0.4)*(y2_d-0.4) < 1100 && runID < 11000)) && chisq2_target < 1.5*chisq2_upstream && chisq2_target < 1.5*chisq2_dump && z2_v < -5 && z2_v > -320 && chisq2/(nHits2-5) < 12 && abs(abs(px2_st1-px2_st3) - 0.416) < 0.008 && abs(py2_st1-py2_st3) < 0.008 && abs(pz2_st1-pz2_st3) < 0.08 && y2_st1/y2_st3<1 && y2_st1*y2_st3>0 && abs(py2_st1)>0.02"

chuckCutsDimuon = "abs(dx) < 0.25 && ((abs(dy-1.6) < 0.22 && runID > 11000) || (abs(dy-0.4) < 0.22 && runID < 11000)) && dz > -280 && dz < -5 && abs(dpx) < 1.8 && abs(dpy) < 2 && dpz > 38 && dpz < 116 && dpx*dpx+dpy*dpy<5 && ((dx*dx+(dy-1.6)*(dy-1.6)<0.06 && runID > 11000) || (dx*dx+(dy-0.4)*(dy-0.4)<0.06 && runID < 11000)) && abs(trackSeparation) < 270 && chisq_dimuon < 18 && abs(chisq1_target+chisq2_target-chisq_dimuon) < 2 && y1_st3*y2_st3<0 && nHits1+nHits2>29 && abs(x1_st1+x2_st1) < 42 && nHits1St1+nHits2St1>8"

physics = "abs(costh)<0.5 && xT>0.05 && xT<.55 && mass>4.2 && mass<9.1"
physics_charm = "abs(costh)<0.5 && xT>0.05 && xT<.55"
# physics = "abs(costh)<0.5 && xT>0.05 && xT<.55 && mass>4.2 && mass<8.8" # Original commented out
# xFCut = "xF>0.0 && xF < 0.05"; # Original commented out

mass99 = "0.99*mass>4.2"
target1 = "targetPos==1"
target3 = "targetPos==3"
flask = "targetPos==2"
occ = "D1<400 && D2<400 && D3<400 && D1+D2+D3<1000"
tInt = "((RFp00-34)*G2SEM/(QIEsum-369000*588*34))>0 && (RFp00-34)*G2SEM/(QIEsum-369000*588*34)<80000"

# Combining the cuts using f-strings for clarity
gmcTemp = f"({chuckCutsPositive}) && ({chuckCutsNegative}) && ({chuckCutsDimuon}) && ({physics})"
gmcTemp_charm = f"({chuckCutsPositive}) && ({chuckCutsNegative}) && ({chuckCutsDimuon}) && ({physics_charm})"

# The TString gmc4pi was commented out and relied on an undefined xFCut in the C++ version.
# If you need it, you'll have to define how xFCut is determined in Python.
# For example, if xFCut were defined elsewhere or passed as a parameter:
# xFCut_example = "xF > 0.1 && xF < 0.2" # Example
# gmc4pi = xFCut_example