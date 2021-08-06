#!/bin/bash
declare -a lumi_array=(1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 120 140 160 180 200) # 300 500 700 1000 2000 3000 5000 7000 10000 20000 30000 50000)
parallel -j37 "python3 compute_Lb2He3_BR_CL.py ALICE3/AnalysisResults_QAMieiPrimari.root ALICE3/AnalysisResults_QASecondari.root CLmany_ALICE3_part_"{}".root --batch --asymcalc --frac_hyper 0.05 --det ALICE3 --lumi_totest "{}" " ::: ${lumi_array[@]}
hadd -f CL90_ALICE3_pp13TeV_y144_reduced.root CLmany_ALICE3_part_*.root
rm CLmany_ALICE3_part_*.root