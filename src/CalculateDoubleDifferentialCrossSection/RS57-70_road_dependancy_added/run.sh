#!/bin/bash

rm *.pdf *.csv *.tex

python3 main.py
ls -1 ls -1 CrossSection_L*.pdf | xargs -I {} python3 addlogo.py {} ./SeaQuestLogo.png
