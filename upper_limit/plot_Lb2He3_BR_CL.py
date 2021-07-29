#!/usr/bin/env python3

'''Script for the plot of the BR CL for Lb->3He
'''


import argparse
from ROOT import (
    TFile,
    gStyle,
    TCanvas,
    TLine,
    TLegend,
    TLatex,
    TColor,
    kRed,
    kAzure
)

# define custom colors to mimic transparency
kAzureMy = TColor.GetFreeColorIndex()
cAzureMy = TColor(kAzureMy, 159./255, 191./255, 223./255, 'kAzureMy', 1.0)
kRedMy = TColor.GetFreeColorIndex()
cRedMy = TColor(kRedMy, 250./255, 153./255, 153./255, 'kRedMy', 1.0)
kGreenMy = TColor.GetFreeColorIndex()
cGreenMy = TColor(kGreenMy, 179./255, 230./255, 179./255, 'kGreenMy', 1.0)
kOrangeMy = TColor.GetFreeColorIndex()
cOrangeMy = TColor(kOrangeMy, 255./255, 204./255, 128./255, 'kOrangeMy', 1.0)


def set_obj_style(obj, marker, color, markersize=0, linewidth=2, fillstyle=0, linestyle=1):
    '''
    Helper method to set style
    '''
    obj.SetLineColor(color)
    obj.SetLineWidth(linewidth)
    obj.SetLineStyle(linestyle)
    if not isinstance(obj, TLine):
        obj.SetMarkerStyle(marker)
        obj.SetMarkerSize(markersize)
        obj.SetMarkerColor(color)
        obj.SetFillColor(color)
        obj.SetFillStyle(fillstyle)


parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument(
    'fileALICE3',
    metavar='text',
    default='ALICE3.root',
    help='input root file with CL vs lumi for ALICE3'
)
parser.add_argument(
    '--fileITS3',
    metavar='text',
    default=None,
    help='input root file with CL vs lumi for ITS3'
)
parser.add_argument(
    '--fileITS2',
    metavar='text',
    default=None,
    help='input root file with CL vs lumi for ITS2'
)
parser.add_argument(
    '--fileITS1',
    metavar='text',
    default=None,
    help='input root file with CL vs lumi for ITS1'
)
args = parser.parse_args()

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.04, 'xy')
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleOffset(1.4, 'y')
gStyle.SetLegendBorderSize(0)

# not very precise, but close
pythia_pred = 2.6e-6
herwig_pred = 1.0e-9

file_names = {
    'ALICE3': args.fileALICE3,
    'ITS3': args.fileITS3,
    'ITS2': args.fileITS2,
    'ITS1': args.fileITS1
}
colors = {
    'ALICE3': kOrangeMy,
    'ITS3': kAzureMy,
    'ITS2': kGreenMy,
    'ITS1': kRedMy
}
gr_cl90 = {}

for det in file_names:
    if file_names[det]:
        infile = TFile.Open(file_names[det])
        gr_cl90[det] = infile.Get('gr_90cl')
        set_obj_style(gr_cl90[det], 0, colors[det], 0, 3, 1000, 1)

line_pythia = TLine(1., pythia_pred, 3.e3, pythia_pred)
set_obj_style(line_pythia, -1, kRed+1, 0, 3, 0, 1)
line_herwig = TLine(1., herwig_pred, 3.e3, herwig_pred)
set_obj_style(line_herwig, -1, kAzure+4, 0, 3, 0, 9)

leg = TLegend(0.65, 0.8-0.04*len(gr_cl90), 0.9, 0.85)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)
leg.SetHeader('90% CL')
for det in gr_cl90:
    leg.AddEntry(gr_cl90[det], det, 'f')

leg_models = TLegend(0.18, 0.15, 0.58, 0.25)
leg_models.SetFillStyle(0)
leg_models.SetTextSize(0.04)
leg_models.AddEntry(line_pythia, 'PYTHIA 8', 'l')
leg_models.AddEntry(line_herwig, 'HERWIG 7', 'l')

lat = TLatex()
lat.SetNDC()
lat.SetTextSize(0.045)
lat.SetTextFont(42)

cCL = TCanvas('cCL', '', 800, 800)
cCL.SetLogy()
cCL.SetLogx()
cCL.DrawFrame(
    1.,
    1.e-11,
    3.e3,
    1.e-3,
    ';#it{L}_{int} (pb^{#minus1}); BR(#Lambda_{b}^{0} #leftarrow ^{3}He + X)'
)
line_pythia.Draw()
line_herwig.Draw()
for det in gr_cl90:
    gr_cl90[det].Draw('3')
leg.Draw()
leg_models.Draw()
lat.DrawLatex(0.18, 0.87, 'ALICE upgrade projection')
cCL.Modified()
cCL.Update()

cCL.SaveAs('Lb2He3_ALICE3_proj_pp.pdf')

input('Press enter to exit')
