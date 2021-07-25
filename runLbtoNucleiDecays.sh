#!/bin/bash
parallel -j10 "root -b -q -l 'SimulateLbtoNucleiDecays.cc+(\"/home/fgrosa/Lb2He3/config/config_file_LbtoHe3_CRmode2.yml\", 1000000000, \"/home/fgrosa/Lb2He3/outputs/Lb2He3_pp13TeV_y5_CRmode2_HepMC2_FONLLmax/AnalysisResults_"{}".root\", \"/home/fgrosa/Lb2He3/outputs/Lb2He3_pp13TeV_y5_CRmode2_HepMC2_FONLLmax/AnalysisResults_"{}".hepmc\", "{}")'" ::: {0..10}
