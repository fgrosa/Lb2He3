#!/bin/bash
FILE="/home/fgrosa/Lb2He3/SimulateLbto3HeDecays_cc_ACLiC_dict_rdict.pcm"
if test -f "$FILE"; then
    rm $FILE
fi
FILE="/home/fgrosa/Lb2He3/SimulateLbto3HeDecays_cc.so"
if test -f "$FILE"; then
    rm $FILE
fi
FILE="/home/fgrosa/Lb2He3/SimulateLbto3HeDecays_cc.d"
if test -f "$FILE"; then
    rm $FILE
fi
parallel -j10 "root -b -q -l 'SimulateLbto3HeDecays.cc+(\"/home/fgrosa/Lb2He3/config/config_file_LbtoHe3_CRmode2.yml\", 1000000000, \"/home/fgrosa/Lb2He3/outputs/Lb2He3_pp13TeV_y5_CRmode2_HepMC2/AnalysisResults_"{}".root\", \"/home/fgrosa/Lb2He3/outputs/Lb2He3_pp13TeV_y5_CRmode2_HepMC2/AnalysisResults_"{}".hepmc\", "{}")'" ::: {1..10}
