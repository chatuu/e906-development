#!/bin/bash
printf "\n Removing existing .pdf files...\n\n"
rm *.pdf

printf "\n Generating plots and cross-sections...\n\n"
python3 Generate2DHists_BinCentrioid.py

printf "\n Adding SeaQuest Logo to the cross-section plot...\n\n"
ls -1 CrossSection_LD2*.pdf | xargs -I {} python3 addlogo.py {} ./SeaQuestLogo.png

printf "\n Copying LH2 and LD2 cross section plots to Tech Note directory...\n\n"
cp CrossSection_LH2_xF_*.pdf ~/github/e906-development/docs/TechNote_LD2_GPS_2026/XSecPlotsBinCentroid/LH2/
cp CrossSection_LD2_xF_*.pdf ~/github/e906-development/docs/TechNote_LD2_GPS_2026/XSecPlotsBinCentroid/LD2/

printf "\n Copying Yield distribution plots to Tech Note directory...\n\n"
cp Y_corrected_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/CorrectedYieldsDists/
cp Y_corrected_Subtracted_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/CorrectedSubtractedYieldsDists/

printf "\n Copying Latex scripts to Tech Note directory...\n\n"
cp Table_PsiP_Contamination.tex /root/github/e906-development/docs/TechNote_LD2_GPS_2026/
cp Appendix_MassCentroids.tex /root/github/e906-development/docs/TechNote_LD2_GPS_2026/



printf "\n All done!\n\n"