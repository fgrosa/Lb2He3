#!/bin/bash
declare -a lumi_array=(1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 120 140 160 180 200)
parallel -j25 "python3 compute_Lb2He3_BR_CL.py ITS3/AnalysisResults_QAMieiPrimari.root ITS3/AnalysisResults_QASecondari.root CL_ITS3_part_"{}".root --batch --asymcalc --frac_hyper 0.05 --det ITS3 --lumi_totest "{}" " ::: ${lumi_array[@]}
hadd -f CL90_ITS3_pp13TeV_y08.root CL_ITS3_part_*.root
rm CL_ITS3_part_*.root