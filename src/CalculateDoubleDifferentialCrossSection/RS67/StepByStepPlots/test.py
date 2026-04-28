(chisq1_target < 15.) & (pz1_st1 > 9.) & (pz1_st1 < 75.) &
(nHits1 > 13) & (x1_t**2 + (y1_t - beam_offset)**2 < 320.) &
(x1_d**2 + (y1_d - beam_offset)**2 < 1100.) &
(x1_d**2 + (y1_d - beam_offset)**2 > 16.) &
(chisq1_target < 1.5 * chisq1_upstream) &
(chisq1_target < 1.5 * chisq1_dump) &
(z1_v < -5.) & (z1_v > -320.) &
(chisq1/(nHits1 - 5) < 12) & ((y1_st1)/(y1_st3 ) < 1.) & 
(np.abs(np.abs(px1_st1 - px1_st3) - 0.416) < 0.008) & 
(np.abs(py1_st1 - py1_st3) < 0.008) &
(np.abs(pz1_st1 - pz1_st3) < 0.08) &
((y1_st1) * (y1_st3) > 0.) & (np.abs(py1_st1) > 0.02)