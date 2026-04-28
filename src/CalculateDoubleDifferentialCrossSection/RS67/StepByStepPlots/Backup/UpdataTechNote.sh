#!/bin/bash

printf "\n Updating Yield Distributions:\n"
cp ./Y_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/YieldDists/

printf "\n Updating Reconstruction Efficiency Distributions:\n"
cp ./E_*_reco_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/RecoEffDists/

printf "\n Updating Hodoscope Efficiency Distributions:\n"
cp ./E_*_hodo_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/HodoEffDists/

printf "\n Updating Total Efficiency Distributions:\n"
cp ./E_*_final_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/FinalEffDists/

printf "\n Updating Signal Efficiency Distributions:\n"
cp ./E_final_signal_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/SigEffDists/

printf "\n Updating Corrected Yield Distributions:\n"
cp ./Y_corrected_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/CorrectedYieldsDists/

printf "\n Updating Corrected Subtracted Yield Distributions:\n"
cp ./Y_corrected_Subtracted_*.pdf /root/github/e906-development/docs/TechNote_LD2_GPS_2026/CorrectedSubtractedYieldsDists/


