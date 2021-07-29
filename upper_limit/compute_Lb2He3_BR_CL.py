#!/usr/bin/env python3

'''Script for the evaluation of the expected 90% CL for the Lb->He3 decay
'''

import argparse
import numpy as np
from ROOT import (
    gROOT,
    gRandom,
    TFile,
    TCanvas,
    TGraph,
    TGraphAsymmErrors,
    TTree,
    TSpline3,
    kError,
    RooFit,
    RooDataSet,
    RooArgSet,
    RooRealVar,
    RooMsgService,
    RooWorkspace,
    RooStats,
    RooRandom,
    RooAbsReal,
    kRed,
    kAzure,
    kGreen,
    kBlack
)


def listo2roo(list_imppar, var, name='data'):

    name = var.GetName()
    x = np.zeros(1, dtype=np.float64)

    tree = TTree('tree', 'tree')
    tree.Branch(f'{name}', x, f'{name}/D')

    for imp_par in list_imppar:
        x[0] = imp_par
        tree.Fill()

    array_roo = RooDataSet(name, 'dataset', tree, RooArgSet(var))
    return array_roo


parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument(
    'infileprimary',
    metavar='text',
    help='AnalysisResults.root file with impact-parameter distribution of primary He3'
)
parser.add_argument(
    'infilefromLb',
    metavar='text',
    help='AnalysisResults.root file with impact-parameter distribution of He3 from Lb decays'
)
parser.add_argument(
    'outfile',
    metavar='text',
    help='root output file name'
)

parser.add_argument('--det', metavar='text', default='ALICE3', help='detector to load proper cross sections')
parser.add_argument('--ptmin', type=float, default=3., help='min pT')
parser.add_argument('--ptmax', type=float, default=4., help='max pT')
parser.add_argument('--batch', action='store_true', help='suppress video output')
parser.add_argument('--frac_hyper', type=float, default=0.01, help='fraction of 3He from hypertriton')
args = parser.parse_args()

if args.batch:
    gROOT.SetBatch(True)

h_crosssec_fromLb = []
for i_file, suffix in enumerate(['', '_FONLLmin', '_FONLLmax']):
    if args.det == 'ALICE3':
        file_name = f'input_xsec/Lb2He3_y144_CRmode2{suffix}.root'
    elif 'ITS' in args.det:
        file_name = f'input_xsec/Lb2He3_y08_CRmode2{suffix}.root'
    infile = TFile.Open()
    h_crosssec_fromLb.append(infile.Get('hHe3FromLb'))
    h_crosssec_fromLb[-1].SetDirectory(0)
    if i_file == 0:
        h_crosssec_prim = infile.Get('hHe3Primary')
        h_br = infile.Get('hBR')
        h_crosssec_prim.SetDirectory(0)
        h_br.SetDirectory(0)
    infile.Close()

input_br = h_br.GetBinContent(1) * h_br.GetBinContent(2) * h_br.GetBinContent(3)

names = ['bkg', 'fromLb', 'hyper']

cols = ['#e41a1c', '#377eb8', '#4daf4a']
h_imppar_vs_pt, htempl, ltempl = [], [], []
for i_file, (file_name, col) in enumerate(zip([args.infileprimary, args.infilefromLb], cols)):
    infile = TFile.Open(file_name)
    h_imppar_vs_pt.append(
        infile.Get('qa-tracking-resolution/impactParameter/impactParameterRPhiVsPt')
    )
    h_imppar_vs_pt[i_file].SetDirectory(0)
    infile.Close()

    ptbinmin = h_imppar_vs_pt[i_file].GetXaxis().FindBin(args.ptmin*1.0001)
    ptbinmax = h_imppar_vs_pt[i_file].GetXaxis().FindBin(args.ptmax*0.9999)

    htempl.append(h_imppar_vs_pt[i_file].ProjectionY(f'htempl{i_file}', ptbinmin, ptbinmax))

htempl.append(htempl[0].Clone('htempl_flat'))
htempl[-1].Reset()
for i_entry in range(int(htempl[0].Integral())):
    htempl[-1].Fill(gRandom.Uniform(-500, 500))  # flat in [-500, 500]

for h in htempl:
    ltempl.append([])
    for i_entry in range(int(h.Integral())):
        if i_entry > 20000:  # for time reasons
            continue
        ltempl[-1].append(h.GetRandom())

RooMsgService.instance().setGlobalKillBelow(RooFit.ERROR)
RooMsgService.instance().setSilentMode(True)
gErrorIgnoreLevel = kError

lumi_to_test = [1., 10., 20., 50., 100., 200.]
if args.det == 'ALICE3':
    lumi_to_test += [300., 400., 500., 1000., 2000., 5000., 10000., 50000]

br_mult_fact_to_test = [1.e-4, 2.e-4, 4.e-4, 6.e-4, 8.e-4,
                        1.e-3, 2.e-3, 4.e-3, 6.e-3, 8.e-3,
                        1.e-2, 2.e-2, 4.e-2, 6.e-2, 8.e-2,
                        1.e-1, 2.e-1, 4.e-1, 6.e-1, 8.e-1,
                        1.e0, 2.e0, 4.e0, 6.e0, 8.e0,
                        1.e+1, 2.e+1, 4.e+1, 6.e+1, 8.e+1]

outfile = TFile(args.outfile, 'recreate')
gr_nsigma, spl_nsigma = [], []
gr_90cl = TGraphAsymmErrors(0)
gr_90cl.SetName('gr_90cl')
cfit = None
for i_lumi, lumi in enumerate(lumi_to_test):
    gr_nsigma.append([])
    spl_nsigma.append([])
    for i_fonll in range(3):
        gr_nsigma[i_lumi].append([])
        spl_nsigma[i_lumi].append([])
        for i_case, case in enumerate(['cent', 'plusunc', 'minusunc']):
            gr_nsigma[i_lumi][i_fonll].append(TGraph(0))
            gr_nsigma[i_lumi][i_fonll][i_case].SetNameTitle(
                f'gr_nsigma_{lumi}_fonll{i_fonll}_{case}',
                ';BR (#Lambda_{b}^{0} #rightarrow ^{3}He + X); n#sigma'
            )

    ptbinmin = h_crosssec_prim.GetXaxis().FindBin(args.ptmin*1.0001)
    ptbinmax = h_crosssec_prim.GetXaxis().FindBin(args.ptmax*0.9999)
    cross_sec_prim = h_crosssec_prim.Integral(ptbinmin, ptbinmax, 'width')
    n_exp_prim = int(round(lumi * cross_sec_prim, 0))
    n_exp_hyper = int((n_exp_prim) * args.frac_hyper)
    cross_sec_fromLb, n_exp_fromLb = [], []
    for hist in h_crosssec_fromLb:
        cross_sec_fromLb.append(hist.Integral(ptbinmin, ptbinmax, 'width'))
        n_exp_fromLb.append(int(round(lumi * cross_sec_fromLb[-1], 0)))

    for i_point, br_mult_fact in enumerate(br_mult_fact_to_test):
        for i_fonll, n_exp_sec in enumerate(n_exp_fromLb):
            print(f'lumi = {lumi:4.0f}, fonll = {i_fonll}, br mult fact = {br_mult_fact:2.6f}')

            roo_ws = RooWorkspace("w")
            roo_imppar = RooRealVar("d0", "#it{d}_{0}", 0., -500., 500., "#mum")

            ktempl, bintempl = [], []
            for i, name in enumerate(names):
                ktempl.append(listo2roo(ltempl[i], roo_imppar, f'data_{name}'))
                getattr(roo_ws, 'import')(ktempl[i], RooFit.Rename(f'data_{name}'))
                # bintempl.append(ktempl[i].binnedClone())
                # getattr(roo_ws, 'import')(bintempl[i], RooFit.Rename(f'data_{name}'))

            n_exp_sec_brcorr = n_exp_sec * br_mult_fact

            roo_ws.factory(f"Poisson::pois(b[{n_exp_prim}, 0., {n_exp_prim*10}],meanb[{n_exp_prim}])")
            roo_ws.factory(f"KeysPdf::templ_prim(d0,data_{names[0]},NoMirror,1)")
            roo_ws.factory(f"KeysPdf::templ_fromLb(d0,data_{names[1]},NoMirror,1)")
            roo_ws.factory(f"KeysPdf::templ_hyper(d0,data_{names[2]},NoMirror,1)")
            roo_ws.factory(f"SUM:model(Nprim[{n_exp_prim}, 0., {n_exp_prim*10}]*templ_prim, Nhyper[{n_exp_hyper}, 0., {n_exp_hyper*10}]*templ_hyper, Nsig[{n_exp_sec_brcorr}, 0., {n_exp_sec_brcorr*1000}]*templ_fromLb)")
            roo_model = roo_ws.pdf("model")
            roo_model_prim = roo_ws.pdf("templ_prim")
            roo_model_hyper = roo_ws.pdf("templ_hyper")
            roo_model_fromLb = roo_ws.pdf("templ_fromLb")

            RooRandom.randomGenerator().SetSeed(42)
            roo_data = roo_model.generate(
                RooArgSet(roo_imppar),
                RooFit.NumEvents(n_exp_prim + n_exp_sec_brcorr + n_exp_hyper)
            )

            # propaganda plot
            if ((args.det != 'ALICE3' and lumi == 200) or
                (args.det == 'ALICE3'and lumi == 1000)) and i_fonll == 0 and br_mult_fact == 1.:
                roo_model.fitTo(roo_data)
                d0_frame = roo_imppar.frame()
                roo_data.plotOn(d0_frame)
                roo_model.plotOn(d0_frame, RooFit.LineColor(kBlack))
                roo_model_prim.plotOn(
                    d0_frame,
                    RooFit.Normalization(roo_ws.var("Nprim").getVal(), RooAbsReal.Raw),
                    RooFit.LineColor(kRed+1)
                )
                roo_model_hyper.plotOn(
                    d0_frame,
                    RooFit.Normalization(roo_ws.var("Nhyper").getVal(), RooAbsReal.Raw),
                    RooFit.LineColor(kGreen+2)
                )
                roo_model_fromLb.plotOn(
                    d0_frame,
                    RooFit.Normalization(roo_ws.var("Nsig").getVal(), RooAbsReal.Raw),
                    RooFit.LineColor(kAzure+4)
                )

                cfit = TCanvas('cfit', 'cfit', 800, 800)
                cfit.SetLogy()
                d0_frame.Draw()
                cfit.SaveAs(args.outfile.replace('.root', '.pdf'))

            roo_ws.defineSet("obs", "d0")
            roo_ws.defineSet("poi", "Nsig")
            roo_ws.defineSet("nuisParams","Nprim,Nhyper")

            sb_model = RooStats.ModelConfig("S+B_model", roo_ws)
            sb_model.SetPdf(roo_model)
            sb_model.SetObservables(roo_ws.set("obs"))
            sb_model.SetParametersOfInterest(roo_ws.set("poi"))
            sb_model.SetNuisanceParameters(roo_ws.set("nuisParams"))
            roo_ws.var("Nsig").setVal(n_exp_sec_brcorr)
            sb_model.SetSnapshot(roo_ws.set("poi"))

            b_model = RooStats.ModelConfig("B_model", roo_ws)
            b_model.SetPdf(roo_model)
            b_model.SetObservables(roo_ws.set("obs"))
            b_model.SetParametersOfInterest(roo_ws.set("poi"))
            b_model.SetNuisanceParameters(roo_ws.set("nuisParams"))
            roo_ws.var("Nsig").setVal(0.0)
            b_model.SetSnapshot(roo_ws.set("poi"))

            test_calc = RooStats.AsymptoticCalculator(roo_data, sb_model, b_model)
            RooStats.AsymptoticCalculator.SetPrintLevel(0)
            test_calc.SetOneSidedDiscovery(True)

            res = test_calc.GetHypoTest()
            signif = res.Significance()
            signif_unc = res.SignificanceError()

            if np.isposinf([signif])[0] or signif > 5:
                continue
            gr_nsigma[i_lumi][i_fonll][0].SetPoint(i_point, input_br * br_mult_fact, signif)
            gr_nsigma[i_lumi][i_fonll][1].SetPoint(i_point, input_br * br_mult_fact, signif + signif_unc)
            gr_nsigma[i_lumi][i_fonll][2].SetPoint(i_point, input_br * br_mult_fact, signif - signif_unc)

    br_90cl = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
    for i_fonll in range(3):
        for i_case, case in enumerate(['cent', 'plusunc', 'minusunc']):
            spl_nsigma[i_lumi][i_fonll].append(
                TSpline3(f'spl_nsigma_{lumi}_fonll{i_fonll}_{case}', gr_nsigma[i_lumi][i_fonll][i_case])
            )
            spl_nsigma[i_lumi][i_fonll][i_case].SetName(f'spl_nsigma_{lumi}_fonll{i_fonll}_{case}')
            outfile.cd()
            signif = -1.
            br = input_br * br_mult_fact_to_test[0]
            while abs(signif - 1.64) > 1.e-2 and br < input_br * br_mult_fact_to_test[-1]:
                signif = spl_nsigma[i_lumi][i_fonll][i_case].Eval(br)
                br += input_br * br_mult_fact_to_test[0]

            br_90cl[3*i_fonll + i_case] = br
            print(f'\nfor lumi = {lumi:4.0f}, fonll = {i_fonll} consider significance ({case}) {signif}, corresponding to br = {br}')

            gr_nsigma[i_lumi][i_fonll][i_case].Write()
            spl_nsigma[i_lumi][i_fonll][i_case].Write()

    br_90cl.sort()
    br_90cl_cent = br_90cl[-1] - br_90cl[0]
    gr_90cl.SetPoint(i_lumi, lumi, br_90cl_cent)
    gr_90cl.SetPointError(i_lumi, 0., 0., br_90cl_cent-br_90cl[0], br_90cl[-1]-br_90cl_cent)

outfile.cd()
gr_90cl.Write()
if cfit:
    cfit.Write()
outfile.Close()
