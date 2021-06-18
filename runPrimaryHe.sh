#!/bin/bash
FILE="/home/fgrosa/Lb2He3/SimulatePrimary3He_cc_ACLiC_dict_rdict.pcm"
if test -f "$FILE"; then
    rm $FILE
fi
parallel -j100 "root -b -q -l 'SimulatePrimary3He.cc+(0, 100000, -2., 2., \"/home/fgrosa/Lb2He3/outputs/PrimaryHe3_pp13TeV_y2_HepMC2/AnalysisResults_"{}".root\", \"/home/fgrosa/Lb2He3/outputs/PrimaryHe3_pp13TeV_y2_HepMC2/AnalysisResults_"{}".hepmc\", "{}")'" ::: {1..100}
