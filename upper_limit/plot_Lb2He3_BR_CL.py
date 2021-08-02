#!/usr/bin/env python3

'''Script for the plot of the BR CL for Lb->3He
'''


import argparse
import ctypes
import numpy as np
from ROOT import (
    TFile,
    gStyle,
    TCanvas,
    TLine,
    TLegend,
    TLatex,
    TColor,
    TGraphAsymmErrors,
    TGraph,
    TGraphSmooth,
    kRed,
    kAzure,
    kGray
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


def set_obj_style(obj, marker, color, markersize=0,
                  linewidth=2, fillstyle=0, linestyle=1, alpha=1.):
    '''
    Helper method to set style
    '''
    if alpha < 1.:
        obj.SetLineColorAlpha(color, alpha)
    else:
        obj.SetLineColor(color)
    obj.SetLineWidth(linewidth)
    obj.SetLineStyle(linestyle)
    if not isinstance(obj, TLine):
        obj.SetMarkerStyle(marker)
        obj.SetMarkerSize(markersize)
        obj.SetFillStyle(fillstyle)
        if alpha < 1.:
            obj.SetMarkerColorAlpha(color, alpha)
            obj.SetFillColorAlpha(color, alpha)
        else:
            obj.SetMarkerColor(color)
            obj.SetFillColor(color)


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
parser.add_argument(
    '--suffix',
    metavar='text',
    default='',
    help='suffix for output file'
)
args = parser.parse_args()

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.04, 'xy')
gStyle.SetPadBottomMargin(0.14)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleOffset(1.4, 'y')
gStyle.SetTitleOffset(1.2, 'x')
gStyle.SetLegendBorderSize(0)

# not very precise, but close
pythia_pred = 2.6e-6
herwig_pred = 1.0e-9

file_names = {
    'ITS1': args.fileITS1,
    'ITS2': args.fileITS2,
    'ITS3': args.fileITS3,
    'ALICE3': args.fileALICE3
}
colors = {
    'ALICE3': kOrangeMy,
    'ITS3': kAzureMy,
    'ITS2': kGreenMy,
    'ITS1': kRedMy
}
gr_cl90 = {}

lumi_to_test = [1., 10., 20., 50., 100., 200., 300.,
                400., 500., 1000., 2000., 5000., 10000., 50000]


line_90cl = TLine(1.e-8, 1.64, 1.e-4, 1.64)
set_obj_style(line_90cl, -1, kGray+1, 0, 1, 0, 9)

gr_signif, spl_signif, c_signif = {}, {}, {}
for lumi in lumi_to_test:
    c_signif[lumi] = TCanvas(f'c_signif_lumi{lumi}', '', 800, 800)
    h_frame = c_signif[lumi].DrawFrame(
        1.e-8, 0., 1.e-4, 5.,
        ';BR(#Lambda_{b}^{0} #rightarrow ^{3}He + X); significance'
    )
    h_frame.GetYaxis().SetDecimals()
    c_signif[lumi].SetLogx()
    line_90cl.Draw()

lat = TLatex()
lat.SetTextSize(0.045)
lat.SetTextFont(42)
lat.SetNDC()

gr_smoother = TGraphSmooth('smoother')
for det in file_names:
    gr_signif[det] = {}
    spl_signif[det] = {}
    is_first_det = False
    if file_names[det]:
        is_first_det = True
        infile = TFile.Open(file_names[det])
        gr_tosmooth = TGraph(0)
        gr_cl90[det] = TGraphAsymmErrors(0)
        set_obj_style(gr_cl90[det], 0, colors[det], 0, 3, 1000, 1, 0.5)
        br_90cl = {}
        for lumi in lumi_to_test:
            gr_signif[det][lumi] = []
            spl_signif[det][lumi] = []
            br_90cl[lumi] = [-1, -1, -1]
            for i_fonll in range(3):
                linestyle = 1
                if i_fonll != 0:
                    linestyle = 2
                gr_signif[det][lumi].append(
                    infile.Get(f'gr_nsigma_{lumi}_fonll{i_fonll}_cent')
                )
                spl_signif[det][lumi].append(
                    infile.Get(f'spl_nsigma_{lumi}_fonll{i_fonll}_cent')
                )
                if gr_signif[det][lumi][i_fonll] and spl_signif[det][lumi][i_fonll]:
                    min_br, max_br, sign = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
                    gr_signif[det][lumi][i_fonll].GetPoint(0, min_br, sign)
                    gr_signif[det][lumi][i_fonll].GetPoint(gr_signif[det][lumi][i_fonll].GetN()-1, max_br, sign)
                    for br in np.arange(min_br.value*2, max_br.value, min_br.value/50):
                        signif = spl_signif[det][lumi][i_fonll].Eval(br)
                        if abs(signif-1.64) < 1.e-2 and min_br.value < br < max_br.value:
                            br_90cl[lumi][i_fonll] = br
                            print(det, lumi, br)
                            break
                    set_obj_style(
                        gr_signif[det][lumi][i_fonll],
                        0, colors[det], 0, 3, 1000, linestyle, 1.
                    )
                    c_signif[lumi].cd()
                    gr_signif[det][lumi][i_fonll].Draw('l')
                    if i_fonll == 0 and is_first_det:
                        lat.DrawLatex(0.6, 0.85, f'#it{{L}}_{{int}} = {lumi:.0f} pb^{{#minus1}}')
                    c_signif[lumi].Modified()
                    c_signif[lumi].Update()

            if all(el > -1. for el in br_90cl[lumi]):
                gr_tosmooth.SetPoint(gr_tosmooth.GetN(), lumi, br_90cl[lumi][0])

        gr_smooth = gr_smoother.SmoothSuper(gr_tosmooth)
        for i_lumi in range(gr_smooth.GetN()):
            lumi, cl = ctypes.c_double(), ctypes.c_double()
            gr_smooth.GetPoint(i_lumi, lumi, cl)
            lumi = lumi.value
            cl = cl.value
            min_br = min(br_90cl[lumi])
            max_br = max(br_90cl[lumi])
            cent = (max_br + min_br) / 2
            if cl < 0.:
                cl = cent
            unc_low = (cent - min_br) / cent * cl
            unc_high = (max_br - cent) / cent * cl
            print(f'{det}, lumi = {lumi} ({gr_cl90[det].GetN()}), br = ({cl, unc_low, unc_high})')
            gr_cl90[det].SetPoint(gr_cl90[det].GetN(), lumi, cl)
            gr_cl90[det].SetPointError(
                gr_cl90[det].GetN()-1,
                0.,
                0.,
                unc_low,
                unc_high)

line_pythia = TLine(1., pythia_pred, 5.e4, pythia_pred)
set_obj_style(line_pythia, -1, kRed+1, 0, 3, 0, 1)
line_pythia_5dot6 = TLine(1., pythia_pred/5.6, 5.e4, pythia_pred/5.6)
set_obj_style(line_pythia_5dot6, -1, kRed+1, 0, 3, 0, 9)
line_pythia_17 = TLine(1., pythia_pred/17, 5.e4, pythia_pred/17)
set_obj_style(line_pythia_17, -1, kRed+1, 0, 3, 0, 2)
line_herwig = TLine(1., herwig_pred, 5.e4, herwig_pred)
set_obj_style(line_herwig, -1, kAzure+4, 0, 3, 0, 6)

leg = TLegend(0.7, 0.85-0.04*len(gr_cl90), 0.9, 0.9)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)
leg.SetHeader('90% CL')
for det in gr_cl90:
    leg.AddEntry(gr_cl90[det], det, 'f')

leg_pythia = TLegend(0.18, 0.24, 0.6, 0.40)
leg_pythia.SetFillStyle(0)
leg_pythia.SetTextSize(0.04)
leg_pythia.SetHeader('PYTHIA8')
leg_pythia.AddEntry(line_pythia, '#Lambda_{b}-tune', 'l')
leg_pythia.AddEntry(line_pythia_5dot6, '#Lambda_{b}-tune / 5.6', 'l')
leg_pythia.AddEntry(line_pythia_17, '#Lambda_{b}-tune / 17', 'l')

leg_herwig = TLegend(0.55, 0.24, 0.97, 0.28)
leg_herwig.SetFillStyle(0)
leg_herwig.SetTextSize(0.04)
leg_herwig.SetHeader('')
leg_herwig.AddEntry(line_herwig, 'HERWIG7', 'l')

lat = TLatex()
lat.SetNDC()
lat.SetTextSize(0.045)
lat.SetTextFont(42)

cCL = TCanvas('cCL', '', 800, 800)
cCL.SetLogy()
cCL.SetLogx()
cCL.DrawFrame(
    1.,
    3.e-10,
    5.e4,
    1.e-2,
    ';#it{L}_{int} (pb^{#minus1}); BR(#Lambda_{b}^{0} #rightarrow ^{3}He + X)'
)
for det in gr_cl90:
    gr_cl90[det].Draw('3')
line_pythia.Draw()
line_pythia_5dot6.Draw()
line_pythia_17.Draw()
line_herwig.Draw()
leg.Draw()
leg_pythia.Draw()
leg_herwig.Draw()
lat.DrawLatex(0.18, 0.89, 'ALICE upgrade projection')
lat.DrawLatex(0.18, 0.83, 'pp, #sqrt{#it{s}} = 14 TeV')
lat.DrawLatex(0.18, 0.77, '3 < #it{p}_{T}  < 4 GeV/#it{c}')
cCL.Modified()
cCL.Update()

cCL.SaveAs(f'Lb2He3_ALICE3_proj_pp{args.suffix}.pdf')

input('Press enter to exit')
