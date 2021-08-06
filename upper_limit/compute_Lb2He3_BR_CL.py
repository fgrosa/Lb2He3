#!/usr/bin/env python3

'''Script for the evaluation of the expected 90% CL for the Lb->He3 decay
'''

import sys
import argparse
import numpy as np
from ROOT import (
    gROOT,
    gRandom,
    gStyle,
    TH1F,
    TFile,
    TCanvas,
    TGraph,
    TGraphAsymmErrors,
    TTree,
    TSpline3,
    TLegend,
    TLatex,
    TGraphSmooth,
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

parser.add_argument('--lumi_totest', nargs='+')
parser.add_argument('--hybcalc', action='store_true', default=False, help='use HybridCalculator for CL calculation')
parser.add_argument('--asymcalc', action='store_true', default=False, help='use AsymptoticCalculator for CL calculation')
parser.add_argument('--freqcalc', action='store_true', default=False, help='use FrequentistCalculator for CL calculation')
parser.add_argument('--custom', action='store_true', default=False, help='use custom method for CL calculation')
parser.add_argument('--det', metavar='text', default='ALICE3', help='detector to load proper cross sections')
parser.add_argument('--ptmin', type=float, default=3., help='min pT')
parser.add_argument('--ptmax', type=float, default=4., help='max pT')
parser.add_argument('--batch', action='store_true', help='suppress video output')
parser.add_argument('--frac_hyper', type=float, default=0.01, help='fraction of 3He from hypertriton')
parser.add_argument('--flat_hyper', action='store_true', default=False, help='use flat distribution for hypertriton')
args = parser.parse_args()

if not args.hybcalc and not args.asymcalc and not args.freqcalc and not args.custom:
    print('ERROR: you should choose the way of computing the CL!')
    sys.exit()

if sum([args.hybcalc, args.asymcalc, args.freqcalc, args.custom]) >= 2:
    print('ERROR: you can choose only one way of computing the CL!')
    sys.exit()

if args.batch:
    gROOT.SetBatch(True)

gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.035)
gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.04, 'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)

h_crosssec_fromLb = []
for i_file, suffix in enumerate(['', '_FONLLmin', '_FONLLmax']):
    if args.det == 'ALICE3':
        file_name = f'input_xsec/Lb2He3_y144_CRmode2{suffix}.root'
    elif 'ITS' in args.det:
        file_name = f'input_xsec/Lb2He3_y08_CRmode2{suffix}.root'
    infile = TFile.Open(file_name)
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
    if i_file == 1 and not args.flat_hyper:
        ptbinmin = h_imppar_vs_pt[i_file].GetXaxis().FindBin(0.*1.0001) # consider hypertriton as low pT Lb
        ptbinmax = h_imppar_vs_pt[i_file].GetXaxis().FindBin(1.*0.9999)
        htempl.append(h_imppar_vs_pt[i_file].ProjectionY('htempl2', ptbinmin, ptbinmax))

if args.flat_hyper:
    htempl.append(htempl[0].Clone('htempl_hyper'))
    htempl[-1].Reset()
    for i_entry in range(int(htempl[0].Integral())):
        htempl[-1].Fill(gRandom.Uniform(-200, 200))  # flat in [-200, 200]

for h in htempl:
    ltempl.append([])
    for i_entry in range(int(h.Integral())):
        if i_entry > 20000:  # for time reasons
            continue
        ltempl[-1].append(h.GetRandom())

RooMsgService.instance().setGlobalKillBelow(RooFit.ERROR)
RooMsgService.instance().setSilentMode(True)
gErrorIgnoreLevel = kError

br_mult_fact_to_test = [1.e-3, 2.e-3, 4.e-3, 6.e-3, 8.e-3,
                        1.e-2, 2.e-2, 4.e-2, 6.e-2, 8.e-2,
                        1.e-1, 2.e-1, 4.e-1, 6.e-1, 8.e-1,
                        1.e0, 2.e0, 4.e0, 6.e0, 8.e0,
                        1.e+1, 2.e+1, 4.e+1, 6.e+1, 8.e+1,
                        1.e+2, 2.e+2, 4.e+2, 6.e+2, 8.e+2]

outfile = TFile(args.outfile, 'recreate')
gr_nsigma, spl_nsigma = [], []
gr_90cl = TGraphAsymmErrors(0)
gr_90cl.SetName('gr_90cl')
cfit = None

gr_smooth = TGraphSmooth('gr_smooth')
for i_lumi, lumi in enumerate(args.lumi_totest):
    lumi = float(lumi)
    signif_lim = [-1., -1., -1.]
    gr_nsigma.append([])
    spl_nsigma.append([])
    for i_fonll in range(3):
        gr_nsigma[i_lumi].append([])
        spl_nsigma[i_lumi].append([])
        for i_case, case in enumerate(['cent', 'plusunc', 'minusunc']):
            gr_nsigma[i_lumi][i_fonll].append(TGraph(0))

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
            print(f'\nlumi = {lumi:4.0f}, fonll = {i_fonll}, br mult fact = {br_mult_fact:2.6f}')

            if all(sign > 5 for sign in signif_lim):  # avoid to keep doing it if we already reached the needed precision
                continue

            if lumi < 100 and br_mult_fact < 1.e-2:  # avoid to test surely not significant cases
                continue

            roo_ws = RooWorkspace('w')
            roo_imppar = RooRealVar('d0', '#it{d}_{0}^{#it{xy}}', 0., -200., 200., '#mum')

            ktempl, bintempl = [], []
            for i, name in enumerate(names):
                ktempl.append(listo2roo(ltempl[i], roo_imppar, f'data_{name}'))
                getattr(roo_ws, 'import')(ktempl[i], RooFit.Rename(f'data_{name}'))
                # bintempl.append(ktempl[i].binnedClone())
                # getattr(roo_ws, 'import')(bintempl[i], RooFit.Rename(f'data_{name}'))

            n_exp_sec_brcorr = n_exp_sec * br_mult_fact

            roo_ws.factory(f'Poisson::pois(x[{n_exp_prim}, 0., {n_exp_prim*100}],mean[{n_exp_prim}, 0., {n_exp_prim*100}])')
            roo_ws.factory(f'KeysPdf::templ_prim(d0,data_{names[0]},NoMirror,2)')
            roo_ws.factory(f'KeysPdf::templ_fromLb(d0,data_{names[1]},NoMirror,2)')
            roo_ws.factory(f'KeysPdf::templ_hyper(d0,data_{names[2]},NoMirror,2)')
            roo_ws.factory(f'SUM:model(Nprim[{n_exp_prim}, 0., {n_exp_prim*100}]*templ_prim,'
                           f'Nhyper[{n_exp_hyper}, 0., {n_exp_hyper*100}]*templ_hyper,'
                           f'Nsig[{n_exp_sec_brcorr}, 0., {n_exp_sec_brcorr*100}]*templ_fromLb)')
            roo_model = roo_ws.pdf('model')
            roo_model_prim = roo_ws.pdf('templ_prim')
            roo_model_hyper = roo_ws.pdf('templ_hyper')
            roo_model_fromLb = roo_ws.pdf('templ_fromLb')

            RooRandom.randomGenerator().SetSeed(42)
            roo_data = roo_model.generate(
                RooArgSet(roo_imppar),
                RooFit.NumEvents(n_exp_prim + n_exp_sec_brcorr + n_exp_hyper)
            )

            # propaganda plot
            if ((args.det != 'ALICE3' and lumi == 200) or
                (args.det == 'ALICE3' and lumi == 1000)) and i_fonll == 0 and br_mult_fact == 1.:
                roo_model.fitTo(roo_data)
                d0_frame = roo_imppar.frame(
                    RooFit.Title('')
                )
                roo_data.plotOn(
                    d0_frame,
                    RooFit.LineWidth(2),
                    RooFit.Name('data')
                )
                roo_model.plotOn(
                    d0_frame,
                    RooFit.LineColor(kBlack),
                    RooFit.Name('total')
                )
                roo_model_prim.plotOn(
                    d0_frame,
                    RooFit.Normalization(roo_ws.var('Nprim').getVal(), RooAbsReal.Raw),
                    RooFit.LineColor(kRed+1),
                    RooFit.LineStyle(2),
                    RooFit.DrawOption('L'),
                    RooFit.Name('primary')
                )
                roo_model_hyper.plotOn(
                    d0_frame,
                    RooFit.Normalization(roo_ws.var('Nhyper').getVal(), RooAbsReal.Raw),
                    RooFit.LineColor(kGreen+2),
                    RooFit.LineStyle(9),
                    RooFit.DrawOption('L'),
                    RooFit.Name('hyper')
                )
                roo_model_fromLb.plotOn(
                    d0_frame,
                    RooFit.Normalization(roo_ws.var('Nsig').getVal(), RooAbsReal.Raw),
                    RooFit.LineColor(kAzure+4),
                    RooFit.LineStyle(6),
                    RooFit.DrawOption('L'),
                    RooFit.Name('Lb')
                )
                d0_frame.SetMinimum(0.9)
                d0_frame.SetMaximum(1.e6)

                cfit = TCanvas('cfit', 'cfit', 800, 800)
                cfit.SetLogy()
                d0_frame.Draw()

                leg = TLegend(0.6, 0.55, 0.9, 0.8)
                leg.SetTextSize(0.04)
                leg.SetFillStyle(0)
                leg.SetBorderSize(0)
                leg.AddEntry('data', 'simulated data', 'p')
                leg.AddEntry('total', 'total fit function', 'l')
                leg.AddEntry('primary', 'primary ^{3}#bar{He}', 'l')
                leg.AddEntry('hyper', '^{3}#bar{He} #leftarrow ^{3}_{#Lambda}#bar{H}^{+}', 'l')
                leg.AddEntry('Lb', '^{3}#bar{He} #leftarrow #bar{#Lambda}_{b}^{0}', 'l')
                leg.Draw()

                lat = TLatex()
                lat.SetNDC()
                lat.SetTextSize(0.045)
                lat.SetTextFont(42)

                if args.det == 'ALICE3':
                    lat.DrawLatex(0.18, 0.89, 'ALICE 3 simulation')
                    lat.DrawLatex(0.18, 0.84, f'pp, #sqrt{{#it{{s}}}} = 14 TeV, |#it{{y}}| < 1.44, #it{{L}}_{{int}} = {lumi/1000:.0f} fb^{{#minus1}}')
                else:
                    lat.DrawLatex(0.18, 0.89, f'ALICE {args.det} simulation')
                    lat.DrawLatex(0.18, 0.84, f'pp, #sqrt{{#it{{s}}}} = 14 TeV, |#it{{y}}| < 0.8, #it{{L}}_{{int}} = {lumi:.0f} pb^{{#minus1}}')
                lat.DrawLatex(0.18, 0.79, f'{args.ptmin:.0f} < #it{{p}}_{{T}} < {args.ptmax:.0f} GeV/#it{{c}}')

                cfit.SaveAs(args.outfile.replace('.root', '.pdf'))

            roo_ws.defineSet('obs', 'd0')
            roo_ws.defineSet('poi', 'Nsig')
            roo_ws.defineSet('nuisParams', 'Nprim,Nhyper')

            sb_model = RooStats.ModelConfig('S+B_model', roo_ws)
            sb_model.SetPdf(roo_model)
            sb_model.SetObservables(roo_ws.set('obs'))
            sb_model.SetParametersOfInterest(roo_ws.set('poi'))
            sb_model.SetNuisanceParameters(roo_ws.set('nuisParams'))
            roo_ws.var('Nsig').setVal(n_exp_sec_brcorr)
            sb_model.SetSnapshot(roo_ws.set('poi'))

            b_model = RooStats.ModelConfig('B_model', roo_ws)
            b_model.SetPdf(roo_model)
            b_model.SetObservables(roo_ws.set('obs'))
            b_model.SetParametersOfInterest(roo_ws.set('poi'))
            b_model.SetNuisanceParameters(roo_ws.set('nuisParams'))
            roo_ws.var('Nsig').setVal(0.0)
            b_model.SetSnapshot(roo_ws.set('poi'))

            if args.hybcalc:
                slrts = RooStats.SimpleLikelihoodRatioTestStat(sb_model.GetPdf(), b_model.GetPdf())
                slrts.SetNullParameters(b_model.GetSnapshot())
                slrts.SetAltParameters(sb_model.GetSnapshot())

                test_calc = RooStats.HybridCalculator(roo_data, sb_model, b_model)
                toy = test_calc.GetTestStatSampler()
                toy.SetTestStatistic(slrts)
                test_calc.SetToys(1000, 1000)
                test_calc.ForcePriorNuisanceAlt(roo_ws.pdf('pois'))
                test_calc.ForcePriorNuisanceNull(roo_ws.pdf('pois'))

            elif args.asymcalc:
                test_calc = RooStats.AsymptoticCalculator(roo_data, sb_model, b_model)
                RooStats.AsymptoticCalculator.SetPrintLevel(3)
                test_calc.SetOneSidedDiscovery(True)

            elif args.freqcalc:
                profll = RooStats.ProfileLikelihoodTestStat(sb_model.GetPdf())
                profll.SetOneSidedDiscovery(True)

                test_calc = RooStats.FrequentistCalculator(roo_data, sb_model, b_model)
                toy = test_calc.GetTestStatSampler()
                toy.SetTestStatistic(profll)
                test_calc.SetToys(1000, 1000)

            if not args.custom:
                res = test_calc.GetHypoTest()
                signif = res.Significance()
                signif_unc = res.SignificanceError()
                res.Print()
            else:
                hDistr = TH1F(f'hDistr_{lumi}_{i_fonll}_{br_mult_fact}', '', 1000, n_exp_sec_brcorr/5, n_exp_sec_brcorr*5)
                for i_gen in range(100):
                    # roo_ws.var('Nhyper').setConstant(True)
                    # roo_ws.var('Nprim').setConstant(True)
                    roo_data = roo_model.generate(
                        RooArgSet(roo_imppar),
                        RooFit.NumEvents(n_exp_prim + n_exp_sec_brcorr + n_exp_hyper)
                    )
                    roo_model.fitTo(roo_data)
                    hDistr.Fill(roo_ws.var('Nsig').getVal())
                if hDistr.GetRMS() > 0.:
                    signif = hDistr.GetMean()/hDistr.GetRMS()
                else:
                    signif = 0.
                signif_unc = 0.

            signif_lim[i_fonll] = signif

            if np.isposinf([signif])[0] or signif > 5 or np.isnan(signif):
                continue

            if gr_nsigma[i_lumi][i_fonll][0].GetN() == 0:
                min_br = input_br * br_mult_fact

            gr_nsigma[i_lumi][i_fonll][0].SetPoint(
                gr_nsigma[i_lumi][i_fonll][0].GetN(),
                input_br * br_mult_fact,
                signif
            )
            gr_nsigma[i_lumi][i_fonll][1].SetPoint(
                gr_nsigma[i_lumi][i_fonll][0].GetN(),
                input_br * br_mult_fact,
                signif + signif_unc
            )
            gr_nsigma[i_lumi][i_fonll][2].SetPoint(
                gr_nsigma[i_lumi][i_fonll][0].GetN(),
                input_br * br_mult_fact,
                signif - signif_unc
            )

    br_90cl = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
    for i_fonll in range(3):
        for i_case, case in enumerate(['cent', 'plusunc', 'minusunc']):
            gr_nsigma[i_lumi][i_fonll][i_case] = gr_smooth.SmoothSuper(gr_nsigma[i_lumi][i_fonll][i_case])
            gr_nsigma[i_lumi][i_fonll][i_case].SetNameTitle(
                f'gr_nsigma_{lumi}_fonll{i_fonll}_{case}',
                ';BR (#Lambda_{b}^{0} #rightarrow ^{3}He + X); n#sigma'
            )
            spl_nsigma[i_lumi][i_fonll].append(
                TSpline3(f'spl_nsigma_{lumi}_fonll{i_fonll}_{case}', gr_nsigma[i_lumi][i_fonll][i_case])
            )
            spl_nsigma[i_lumi][i_fonll][i_case].SetName(f'spl_nsigma_{lumi}_fonll{i_fonll}_{case}')
            outfile.cd()
            signif = -1.
            br = min_br
            while abs(signif - 1.64) > 1.e-3 and br < input_br * br_mult_fact_to_test[-1]:
                signif = spl_nsigma[i_lumi][i_fonll][i_case].Eval(br)
                br += input_br * br_mult_fact_to_test[0]

            br_90cl[3*i_fonll + i_case] = br
            print(f'\nfor lumi = {lumi:4.0f}, fonll = {i_fonll} consider significance ({case}) {signif}, corresponding to br = {br}')

            gr_nsigma[i_lumi][i_fonll][i_case].Write()
            spl_nsigma[i_lumi][i_fonll][i_case].Write()

    br_90cl.sort()
    br_90cl_cent = br_90cl[-1] - br_90cl[0]
    gr_90cl.SetPoint(gr_90cl.GetN(), lumi, br_90cl_cent)
    gr_90cl.SetPointError(
        gr_90cl.GetN()-1,
        0.,
        0.,
        br_90cl_cent-br_90cl[0],
        br_90cl[-1]-br_90cl_cent
    )

outfile.cd()
gr_90cl.Write()
if cfit:
    cfit.Write()
outfile.Close()
