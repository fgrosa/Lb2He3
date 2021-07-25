#!/bin/bash
parallel -j10 "root -b -q -l 'SimulatePrimary3He.cc+(0, 25000, -0.8, 0.8, \"/home/fgrosa/Lb2He3/outputs/PrimaryHe3_pp13TeV_y08_HepMC2/AnalysisResults_"{}".root\", \"/home/fgrosa/Lb2He3/outputs/PrimaryHe3_pp13TeV_y08_HepMC2/AnalysisResults_"{}".hepmc\", true, "{}")'" ::: {1..400}
