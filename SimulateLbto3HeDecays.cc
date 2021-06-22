#define __HEPMC2__ //__HEPMC2__ or __HEPMC3__

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include "Pythia8/Pythia.h"

#ifdef __HEPMC2__
#include "Pythia8Plugins/HepMC2.h"
#endif

#ifdef __HEPMC3__
#include "Pythia8Plugins/HepMC3.h"
#endif

#endif

using namespace Pythia8;

enum pdgNuclei
{
    kPdgHe3 = 1000020030,
};

static const double massHe3 = 2.80839160743;

enum FONLLPred {
    kCentral, 
    kMin,
    kMax
};

//__________________________________________________________________________________________________
void SimulateLbto3HeDecays(std::string cfgFileName, int nEvents=1000000, std::string outFileNameRoot="AnalysisResults.root", std::string outFileNameHepMC="AnalysisResults.hepmc", int seed=42);
bool SimpleCoalescence(double p1[3], double p2[3], double p3[3], double pHe3[3], double pc = 0.2 /* GeV/c */);
std::array<int, 5> CountNumberOfDaughters(std::vector<int> pdg, std::vector<int> status);
TH1F* ReadFONLL(std::string txtFileName, int whichPred, double varMin, double binWidth, int nBins, std::string var = "pt");

template <typename T>
int findPtBin(std::vector<T> binMins, std::vector<T> binMaxs, T value);

//__________________________________________________________________________________________________
void SimulateLbto3HeDecays(std::string cfgFileName, int nEvents, std::string outFileNameRoot, std::string outFileNameHepMC, int seed)
{
    //__________________________________________________________
    //Load configs from yaml file
    YAML::Node config = YAML::LoadFile(cfgFileName.data());
    if (config.IsNull())
    {
        std::cerr << "\033[31mERROR: yaml config file not found! Exit\033[0m" << std::endl;
        return;
    }
    else 
    {
        std::cout << "\n\n*******************************************" << std::endl;
        std::cout << Form("\033[32mLoading configuration from file %s\033[0m\n", cfgFileName.data()) << std::endl;
    }

    std::string process = config["generator"]["process"].as<std::string>();
    std::string tune = config["generator"]["process"].as<std::string>();

    std::string FONLLPtFileName = config["FONLL"]["pt"]["filename"].as<std::string>();
    double ptMinFONLL = config["FONLL"]["pt"]["ptmin"].as<double>();
    double ptMaxFONLL = config["FONLL"]["pt"]["ptmax"].as<double>();
    int ptBinsFONLL = config["FONLL"]["pt"]["ptbins"].as<int>();
    std::vector<std::string> FONLLyFileNames = config["FONLL"]["y"]["filenames"].as<std::vector<std::string>>();
    std::vector<double> FONLLyPtMins = config["FONLL"]["y"]["ptmins"].as<std::vector<double>>();
    std::vector<double> FONLLyPtMaxs = config["FONLL"]["y"]["ptmaxs"].as<std::vector<double>>();
    double yMinFONLL = config["FONLL"]["y"]["ymin"].as<double>();
    double yMaxFONLL = config["FONLL"]["y"]["ymax"].as<double>();
    int yBinsFONLL = config["FONLL"]["y"]["ybins"].as<int>();

    double coalRadius = config["coalescence"]["radius"].as<double>();
    double coalMomRadius = config["coalescence"]["momentum_radius"].as<double>();
    bool outputRoot = config["output"]["root"].as<bool>();
    bool outputHepMC = config["output"]["hepmc"].as<bool>();
    double yMin = config["acceptance"]["ymin"].as<double>();
    double yMax = config["acceptance"]["ymax"].as<double>();

    if (FONLLyFileNames.size() != FONLLyPtMins.size() || FONLLyFileNames.size() != FONLLyPtMaxs.size())
    {
        std::cerr << "\033[31mERROR: inputs for FONLL dsigma/dy not consistent! Exit\033[0m" << std::endl;
        return;
    }

    //__________________________________________________________
    // create and configure pythia generator
    Pythia pythia;
    if(process == "SoftQCD")
        pythia.readString("SoftQCD:all = on");
    else if(process == "HardQCD")
    {
        pythia.readString("HardQCD:hardccbar = on");
        pythia.readString("HardQCD:hardbbbar = on");
    }

    // set tune
    if(tune == "Monash")
    {
        pythia.readString(Form("Tune:pp = 14"));
    }
    else if(tune == "CRMode0")
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 2.9");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.43");
        pythia.readString("ColourReconnection:timeDilationMode = 0");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation =5");
    }
    else if(tune == "CRMode2")
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.20");
        pythia.readString("ColourReconnection:timeDilationMode = 2");
        pythia.readString("ColourReconnection:timeDilationPar = 0.18");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation =5");
    }
    else if(tune == "CRMode3")
    {
        pythia.readString(Form("Tune:pp = 14"));
        pythia.readString("ColourReconnection:mode = 1");
        pythia.readString("ColourReconnection:allowDoubleJunRem = off");
        pythia.readString("ColourReconnection:m0 = 0.3");
        pythia.readString("ColourReconnection:allowJunctions = on");
        pythia.readString("ColourReconnection:junctionCorrection = 1.15");
        pythia.readString("ColourReconnection:timeDilationMode = 3");
        pythia.readString("ColourReconnection:timeDilationPar = 0.073");
        pythia.readString("StringPT:sigma = 0.335");
        pythia.readString("StringZ:aLund = 0.36");
        pythia.readString("StringZ:bLund = 0.56");
        pythia.readString("StringFlav:probQQtoQ = 0.078");
        pythia.readString("StringFlav:ProbStoUD = 0.2");
        pythia.readString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.readString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.readString("BeamRemnants:remnantMode = 1");
        pythia.readString("BeamRemnants:saturation =5");
    }

    // keep only interesting decays, to be reweighted a posteriori
    pythia.readString("5122:onMode = off");
    pythia.readString("5122:onIfMatch = 2 1 2 2101"); // bRatio="0.0120000" dominant one according to https://arxiv.org/pdf/2006.16251.pdf
    pythia.readString("5122:tau0=4.41000e-01"); // bRatio="0.0800000"
    // pythia.readString("5122:onIfMatch = 2 1 4 2101"); // bRatio="0.4411147"
    // pythia.readString("5122:onIfMatch = 2 4 1 2101"); // bRatio="0.0910000"
    // pythia.readString("5122:onIfMatch = 4 3 2 2101"); // bRatio="0.0120000"
    // pythia.readString("5122:onIfMatch = 4 3 4 2101"); // bRatio="0.0800000"

    // add 3He
    pythia.particleData.addParticle(kPdgHe3, "3He++", "3He--", 2, 6, 0, massHe3, 0., massHe3, massHe3, 1.e9);   

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.init();

    //__________________________________________________________
    // perform the simulation

#ifdef __HEPMC2__
    HepMC::Pythia8ToHepMC ToHepMC;
    HepMC::IO_GenEvent ascii_io(outFileNameHepMC.data(), std::ios::out);
#endif
#ifdef __HEPMC3__
    define HepMC output
    HepMC3::WriterAscii outFileHepMC(outFileNameHepMC);
    HepMC3::Pythia8ToHepMC3 ToHepMC;
#endif

    // define histograms
    auto hBR = new TH1F("hBR", ";;BR", 3, 0.5, 3.5);
    hBR->GetXaxis()->SetBinLabel(1, "#Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}");
    hBR->GetXaxis()->SetBinLabel(2, "#Lambda_{b}^{0} #rightarrow (>=2)p + (>=1)n / #Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}");
    hBR->GetXaxis()->SetBinLabel(3, "#Lambda_{b}^{0} #rightarrow ^{3}He + X / #Lambda_{b}^{0} #rightarrow (>=2)p + (>=1)n");
    hBR->SetBinContent(1, 0.0120000);

    auto hDecayChannel = new TH1F("hDecayChannel", "BR", 16, 0.5, 16.5);
    std::string decayLabel = "^{3}He 2#bar{p}";
    std::string decayLabelN = "^{3}He #bar{p} #bar{n}";
    hDecayChannel->GetXaxis()->SetBinLabel(1, Form("%s", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(2, Form("%s #pi^{0}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(3, Form("%s 2#pi^{0}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(4, Form("%s #pi^{+} #pi^{#minus}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(5, Form("%s 2#pi^{+} 2#pi^{#minus}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(6, Form("%s #pi^{+} #pi^{#minus} #pi^{0}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(7, Form("%s #pi^{+} #pi^{#minus} 2#pi^{0}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(8, Form("%s 2#pi^{+} 2#pi^{#minus} #pi^{0}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(9, Form("%s 2#pi^{+} 2#pi^{#minus} 2#pi^{0}", decayLabel.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(10, Form("%s #pi^{#minus}", decayLabelN.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(11, Form("%s #pi^{#minus} #pi^{0}", decayLabelN.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(12, Form("%s #pi^{#minus} 2#pi^{0}", decayLabelN.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(13, Form("%s 2#pi^{#minus} #pi^{+}", decayLabelN.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(14, Form("%s 2#pi^{#minus} #pi^{+} #pi^{0}", decayLabelN.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(15, Form("%s 2#pi^{#minus} #pi^{+} 2#pi^{0}", decayLabelN.data()));
    hDecayChannel->GetXaxis()->SetBinLabel(16, "other");

    auto hFONLLLbVsPt = ReadFONLL(FONLLPtFileName, kCentral, ptMinFONLL, (ptMaxFONLL-ptMinFONLL)/ptBinsFONLL, ptBinsFONLL);
    hFONLLLbVsPt->Scale(0.816); // f(b->B) from e+e- provides good normalisation for LHCb and CMS B-meson measurements
    hFONLLLbVsPt->SetNameTitle(Form("hFONLLLbVsPt_y_%.1f_%.1f", yMinFONLL, yMaxFONLL), Form(";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T}|_{%.1f < #it{y} < %.1f} (pb GeV^{-1} #it{c})", yMinFONLL, yMaxFONLL));

    std::vector<TH1F*> hFONLLLbVsY{};
    for(size_t iPt=0; iPt<FONLLyPtMins.size(); iPt++)
    {
        hFONLLLbVsY.push_back(ReadFONLL(FONLLyFileNames[iPt], kCentral, yMinFONLL, (yMaxFONLL-yMinFONLL)/yBinsFONLL, yBinsFONLL, "y"));
        hFONLLLbVsY[iPt]->SetNameTitle(Form("hFONLLLbVsY_pT_%.0f_%.0f", FONLLyPtMins[iPt], FONLLyPtMaxs[iPt]), ";#it{y};d#sigma/d#it{y} (pb)");
    }

    auto hLbDecayLengthVsPt = new TH2F("hLbDecayLengthVsPt", ";#it{p}_{T} (GeV/#it{c});decay length (#mum)", 2001, -0.025, 100., 1000, 0., 10000.);
    auto hHe3ProdVtxVsLbPt = new TH2F("hHe3ProdVtxVsLbPt", ";#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});^{3}He prod vtx (#mum)", 2001, -0.025, 100., 1000, 0., 10000.);
    auto hHe3ProdVtxVsHe3Pt = new TH2F("hHe3ProdVtxVsHe3Pt", ";#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});^{3}He prod vtx (#mum)", 1001, -0.025, 50., 1000, 0., 10000.);

    // f(b -> Lb) / f(b -> B) from LHCb measurement https://arxiv.org/pdf/1902.06794.pdf
    TF1* fFFLHCb = new TF1("fracLb","([4] * ([5] + exp([6] + [7] * x))) /  (([0] * ([1] + [2] * (x - [3])))  + ([4] * ([5] + exp([6] + [7] * x))) + 1)  ", 0, 50);
    double parLbA = 1;
    double parLbp1 = 0.0793;
    double parLbp2 = -1.022;
    double parLbp3 = -0.107;
    double parBsA = 1;
    double parBsp1 = 0.119;
    double parBsp2 = -0.00091;
    double parBsAvePt = 10.1;
    fFFLHCb->SetParameters(parBsA, parBsp1, parBsp2, parBsAvePt, parLbA, parLbp1, parLbp2, parLbp3);

    for(int iPt=1; iPt<hFONLLLbVsPt->GetNbinsX()+1; iPt++)
    {
        double ptCent = hFONLLLbVsPt->GetBinCenter(iPt);
        hFONLLLbVsPt->SetBinContent(iPt, hFONLLLbVsPt->GetBinContent(iPt) * fFFLHCb->Eval(ptCent>5 ? ptCent:5));
    }

    auto hHe3FromLb = (TH1F*)hFONLLLbVsPt->Clone("hHe3FromLb_y05");
    hHe3FromLb->Reset();

    auto hPtLbVsHe3FromLb = new TH2F("hPtLbVsHe3FromLb", ";#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(^{3}He) (GeV/#it{c})", 1000, 0., 100., 500., 0., 50.);
    auto hYLbVsHe3FromLb = new TH2F("hYLbVsHe3FromLb", ";#it{y}(^{3}He);#it{y}(^{3}He)", 100, -5., 5., 100., -5., 5.);
    auto hLbPtVsY = new TH2F("hLbPtVsY", ";#it{p}_{T} (GeV/#it{c});#it{y})", 1000, 0., 100., 100., -5., 5.);

    std::vector<double> pxDau, pyDau, pzDau, ptDau;
    std::vector<int> pdgDau, pdgDauAll, statusDauAll, labDau;
    double pxLb, pyLb, pzLb, ptLb;

    int pdgHb = 5122;
    double massLb = TDatabasePDG::Instance()->GetParticle(pdgHb)->Mass();
    int nEventSel = 0, nCoalesced = 0;
    auto fDecay = new TF1("fDecay", "exp(-x/0.441)", 0., 1000.);
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        ptLb = hFONLLLbVsPt->GetRandom();
        int iPtForY = findPtBin(FONLLyPtMins, FONLLyPtMaxs, ptLb);
        double yLb = hFONLLLbVsY[iPtForY]->GetRandom();
        double phiLb = gRandom->Rndm() * 2 * TMath::Pi();
        pxLb = ptLb * TMath::Cos(phiLb);
        pyLb = ptLb * TMath::Sin(phiLb);
        double mt = TMath::Sqrt(massLb * massLb + ptLb * ptLb);
        pzLb = mt * TMath::SinH(yLb);
        double pLb = TMath::Sqrt(ptLb * ptLb + pzLb * pzLb);
        double ELb = TMath::Sqrt(massLb * massLb + pLb * pLb);

        // Lb
        Particle Hb;
        Hb.id(pdgHb);
        Hb.status(81);
        Hb.m(massLb);
        Hb.xProd(0.);
        Hb.yProd(0.);
        Hb.zProd(0.);
        Hb.tProd(0.);
        Hb.e(ELb);
        Hb.px(pxLb);
        Hb.py(pyLb);
        Hb.pz(pzLb);
        Hb.tau(fDecay->GetRandom());

        pythia.event.reset();
        pythia.event.append(Hb);
        int idPart = pythia.event[1].id();
        pythia.particleData.mayDecay(idPart, true);
        pythia.moreDecays();

        if(iEvent%1000000 == 0)
            std::cout << Form("Lb decay number %10d\r", iEvent);

        int nProtons = 0, nAntiProtons = 0, nNeutrons = 0, nAntiNeutrons = 0;
        double decVtx[4] = {0., 0., 0., 0.};
        double decLen = 0.;
        int nSkipped = 0;
        for (int iPart = 2; iPart < pythia.event.size(); iPart++)
        {
            if(iPart == 2)
            {
                decVtx[0] = pythia.event.at(iPart).xProd();
                decVtx[1] = pythia.event.at(iPart).yProd();
                decVtx[2] = pythia.event.at(iPart).zProd();
                decVtx[3] = pythia.event.at(iPart).tProd();
                decLen = TMath::Sqrt(decVtx[0]*decVtx[0] + decVtx[1]*decVtx[1] + decVtx[2]*decVtx[2]);
                hLbDecayLengthVsPt->Fill(ptLb, decLen*1000);
            }
            else
            {
                pdgDauAll.push_back(pythia.event.at(iPart).id());
                statusDauAll.push_back(pythia.event.at(iPart).status());
                double prodR = TMath::Sqrt(pythia.event.at(iPart).xProd()*pythia.event.at(iPart).xProd() + pythia.event.at(iPart).yProd()*pythia.event.at(iPart).yProd() + pythia.event.at(iPart).zProd()*pythia.event.at(iPart).zProd());
                if(TMath::Abs(prodR-decLen) > coalRadius) {
                    nSkipped++;
                    continue;
                }
            }

            double px = pythia.event.at(iPart).px();
            double py = pythia.event.at(iPart).py();
            double pz = pythia.event.at(iPart).pz();
            double pt = TMath::Sqrt(px*px + py*py);

            ptDau.push_back(pt);
            pxDau.push_back(px);
            pyDau.push_back(py);
            pzDau.push_back(pz);
            pdgDau.push_back(pythia.event.at(iPart).id());
            labDau.push_back(iPart-nSkipped);

            if(pdgDau[iPart-2] == 2212)
                nProtons++;
            else if(pdgDau[iPart-2] == 2112)
                nNeutrons++;
        }
        if(nProtons >= 2 && nNeutrons >= 1)
        {
            double pCoal1[3], pCoal2[3], pCoal3[3];
            int labCoal[3] = {-1, -1, -1};
            double pCoalHe3[3] = {-999., -999., -999.};
            bool isSecond = false;
            for(size_t iPart=0; iPart<pdgDau.size(); iPart++)
            {
                if(pdgDau[iPart] == 2212) {
                    if(!isSecond) {
                        pCoal1[0] = pxDau[iPart];
                        pCoal1[1] = pyDau[iPart];
                        pCoal1[2] = pzDau[iPart];
                        isSecond = true;
                        labCoal[0] = labDau[iPart];
                    }
                    else {
                        pCoal2[0] = pxDau[iPart];
                        pCoal2[1] = pyDau[iPart];
                        pCoal2[2] = pzDau[iPart];
                        labCoal[1] = labDau[iPart];
                    }
                }
                else if(pdgDau[iPart] == 2112) {
                    pCoal3[0] = pxDau[iPart];
                    pCoal3[1] = pyDau[iPart];
                    pCoal3[2] = pzDau[iPart];
                    labCoal[2] = labDau[iPart];
                }
            }
            bool hasCoalesced = SimpleCoalescence(pCoal1, pCoal2, pCoal3, pCoalHe3, coalMomRadius);
            if(hasCoalesced) {
                double ptHe3 = TMath::Sqrt(pCoalHe3[0]*pCoalHe3[0] + pCoalHe3[1]*pCoalHe3[1]);
                double pHe3 = TMath::Sqrt(ptHe3*ptHe3 + pCoalHe3[2]*pCoalHe3[2]);
                double EHe3 = TMath::Sqrt(massHe3 * massHe3 + pHe3 * pHe3);
                double yHe3 = TMath::Log((EHe3 + pCoalHe3[2]) / (EHe3 - pCoalHe3[2])) / 2;
                hPtLbVsHe3FromLb->Fill(ptLb, ptHe3);
                hYLbVsHe3FromLb->Fill(yLb, yHe3);

                auto nParticles = CountNumberOfDaughters(pdgDauAll, statusDauAll);
                if(nParticles[3] == 2)
                {
                    if(nParticles[0] == 0 && nParticles[1] == 0 && nParticles[2] == 0)
                        hDecayChannel->Fill(1);
                    else if(nParticles[0] == 1 && nParticles[1] == 0 && nParticles[2] == 0)
                        hDecayChannel->Fill(2);
                    else if(nParticles[0] == 2 && nParticles[1] == 0 && nParticles[2] == 0)
                        hDecayChannel->Fill(3);
                    else if(nParticles[0] == 0 && nParticles[1] == 1 && nParticles[2] == 1)
                        hDecayChannel->Fill(4);
                    else if(nParticles[0] == 0 && nParticles[1] == 2 && nParticles[2] == 2)
                        hDecayChannel->Fill(5);
                    else if(nParticles[0] == 1 && nParticles[1] == 1 && nParticles[2] == 1)
                        hDecayChannel->Fill(6);
                    else if(nParticles[0] == 2 && nParticles[1] == 1 && nParticles[2] == 1)
                        hDecayChannel->Fill(7);
                    else if(nParticles[0] == 1 && nParticles[1] == 2 && nParticles[2] == 2)
                        hDecayChannel->Fill(8);
                    else if(nParticles[0] == 2 && nParticles[1] == 2 && nParticles[2] == 2)
                        hDecayChannel->Fill(9);
                }
                else if(nParticles[3] == 1 && nParticles[4] == 1)
                {
                    if(nParticles[0] == 0 && nParticles[1] == 0 && nParticles[2] == 1)
                        hDecayChannel->Fill(10);
                    else if(nParticles[0] == 1 && nParticles[1] == 1 && nParticles[2] == 0)
                        hDecayChannel->Fill(11);
                    else if(nParticles[0] == 2 && nParticles[1] == 1 && nParticles[2] == 0)
                        hDecayChannel->Fill(12);
                    else if(nParticles[0] == 0 && nParticles[1] == 2 && nParticles[2] == 1)
                        hDecayChannel->Fill(13);
                    else if(nParticles[0] == 1 && nParticles[1] == 2 && nParticles[2] == 1)
                        hDecayChannel->Fill(14);
                    else if(nParticles[0] == 2 && nParticles[1] == 2 && nParticles[2] == 1)
                        hDecayChannel->Fill(15);
                }
                else 
                    hDecayChannel->Fill(16);

                if(yHe3 >= yMin && yHe3 <= yMax) {
                    hHe3FromLb->Fill(ptHe3);
                }

                Particle He3;
                He3.id(kPdgHe3);
                He3.status(93);
                He3.m(massHe3);
                He3.xProd(decVtx[0]);
                He3.yProd(decVtx[1]);
                He3.zProd(decVtx[2]);
                He3.tProd(decVtx[3]);
                He3.e(EHe3);
                He3.px(pCoalHe3[0]);
                He3.py(pCoalHe3[1]);
                He3.pz(pCoalHe3[2]);
                He3.mother1(1); // Lb
                He3.mother2(0);
                He3.daughter1(0);
                He3.daughter2(0);
                He3.tau(1.e9); // stable

                double prodVtx = TMath::Sqrt(He3.xProd()*He3.xProd() + He3.yProd()*He3.yProd() + He3.zProd()*He3.zProd());
                hHe3ProdVtxVsLbPt->Fill(ptLb, prodVtx*1000);
                hHe3ProdVtxVsHe3Pt->Fill(ptHe3, prodVtx*1000);

                for(int iLab=0; iLab<3; iLab++)
                    pythia.event.remove(labCoal[iLab]-iLab, labCoal[iLab]-iLab); // labels shift by 1 when removing a particle
                pythia.event.append(He3);

                // write HepMC
                if(outputHepMC) {
#ifdef __HEPMC3__
                    HepMC3::GenEvent hepmcevt;
                    ToHepMC.fill_next_event(pythia.event, &hepmcevt, -1, &pythia.info);
                    outFileHepMC.write_event(hepmcevt);
#endif
#ifdef __HEPMC2__
                    HepMC::GenEvent *hepmcevt = new HepMC::GenEvent();
                    ToHepMC.fill_next_event(pythia, hepmcevt);
                    ascii_io << hepmcevt;
                    delete hepmcevt;
#endif
                }
                nCoalesced++;
            }
            nEventSel++;
        }

        ptDau.clear();
        pxDau.clear();
        pyDau.clear();
        pzDau.clear();
        pdgDau.clear();
        pdgDauAll.clear();
        statusDauAll.clear();
    }
#ifdef __HEPMC3__
    outFileHepMC.close();
#endif

    hBR->SetBinContent(2, static_cast<double>(nEventSel) / nEvents);
    hBR->SetBinContent(3, static_cast<double>(nCoalesced) / nEventSel);
    hHe3FromLb->Scale(hFONLLLbVsPt->Integral() / nEventSel * hBR->GetBinContent(1) * hBR->GetBinContent(2));
    hDecayChannel->Scale(hBR->GetBinContent(1) * hBR->GetBinContent(2) * hBR->GetBinContent(3) / hDecayChannel->Integral());

    // Save histogram on file and close file.
    if(outputRoot) {
        TFile outFile(outFileNameRoot.data(), "recreate");
        outFile.cd();
        hFONLLLbVsPt->Write();
        for(auto &histo: hFONLLLbVsY)
            histo->Write();
        hHe3FromLb->Write();
        hBR->Write();
        hDecayChannel->Write();
        hPtLbVsHe3FromLb->Write();
        hYLbVsHe3FromLb->Write();
        hLbPtVsY->Write();
        hLbDecayLengthVsPt->Write();
        hHe3ProdVtxVsLbPt->Write();
        hHe3ProdVtxVsHe3Pt->Write();
        outFile.Close();
    }
}

//__________________________________________________________________________________________________
bool SimpleCoalescence(double p1[3], double p2[3], double p3[3], double pHe3[3], double pc) {

    // returns true if coalescence is realised, false otherwise
    TLorentzVector Ptrack1, Ptrack2, Ptrack3, trackSum, Ptrack1CMS, Ptrack2CMS, Ptrack3CMS;
    Ptrack1.SetXYZM(p1[0], p1[1], p1[2], TDatabasePDG::Instance()->GetParticle(2212)->Mass());
    Ptrack2.SetXYZM(p2[0], p2[1], p2[2], TDatabasePDG::Instance()->GetParticle(2212)->Mass());
    Ptrack3.SetXYZM(p3[0], p3[1], p3[2], TDatabasePDG::Instance()->GetParticle(2112)->Mass());
    trackSum = Ptrack1 + Ptrack2 + Ptrack3;

    double beta = trackSum.Beta();
    double betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
    double betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
    double betaz = beta * cos(trackSum.Theta());

    Ptrack1CMS = Ptrack1;
    Ptrack2CMS = Ptrack2;
    Ptrack3CMS = Ptrack3;

    Ptrack1CMS.Boost(-betax, -betay, -betaz);
    Ptrack2CMS.Boost(-betax, -betay, -betaz);
    Ptrack3CMS.Boost(-betax, -betay, -betaz);

    double coalRadius = TMath::Power(2, 1./6) * pc / 2;
    if(Ptrack1CMS.P() <= coalRadius && Ptrack2CMS.P() <= coalRadius && Ptrack3CMS.P() <= coalRadius)
    {
        for(int iEl=0; iEl<3; iEl++)
            pHe3[iEl] = p1[iEl] + p2[iEl] + p3[iEl];
        return true;
    }

    for(int iEl=0; iEl<3; iEl++)
        pHe3[0] = -1.;
    return false;
}


//__________________________________________________________________________________________________
std::array<int, 5> CountNumberOfDaughters(std::vector<int> pdg, std::vector<int> status)
{
    std::array<int, 5> nPart = {0, 0, 0, 0, 0};
    for(size_t iPart=0; iPart<pdg.size(); iPart++) {
        if(status[iPart] == 91) {
            switch(pdg[iPart])
            {
                case -211:
                    nPart[1]++;
                break;
                case 211:
                    nPart[2]++;
                break;
                case -2212:
                    nPart[3]++;
                break;
                case -2112:
                    nPart[4]++;
                break;
            }
        }
        else if(status[iPart] == -91 && pdg[iPart] == 111)
            nPart[0]++;
    }

    return nPart;
}

//__________________________________________________________________________________________________
TH1F* ReadFONLL(std::string txtFileName, int whichPred, double varMin, double binWidth, int nBins, std::string var)
{
    if(txtFileName.find("txt") == string::npos && txtFileName.find("dat") == string::npos && txtFileName.find("csv") == string::npos) {
        std::cerr << "\033[31mERROR: Wrong file format! Exit.\033[0m" << std::endl;
        return nullptr;
    }
    if(var.compare("pt") != 0 && var.compare("y") != 0) {
        std::cerr << Form("\033[31mERROR: Invalid variable %s for FONLL predictions! Exit.\033[0m", var.data()) << std::endl;
        return nullptr;
    }

    std::ifstream inSet(txtFileName.data());
    if(!inSet) {
        std::cerr << "\033[31mERROR: Please check if "<< txtFileName.data() << " is the right path. Exit.\033[0m" << std::endl;
        return nullptr;
    }

    std::vector<std::string> values;
    
    auto hFONLL = new TH1F("hFONLL", "", nBins, varMin-binWidth/2, (varMin + nBins*binWidth) +binWidth/2);
    double centY = -1., minY = -1., maxY = -1.;
    int iBin = 1;
    while(!inSet.eof())
    {    
        std::vector<std::string> values;
        std::string line;
        std::getline(inSet, line);
        if(line.find("#") != std::string::npos || line.find(var.data()) != std::string::npos)
            continue;

        size_t pos = 0;
        while((pos = line.find(" ")) != std::string::npos)
        {
            values.push_back(line.substr(0, pos));
            if(values.size()==2) {
                std::stringstream convert(values[1]);
                if( !(convert >> centY) )
                    centY = -1;
            }
            else if(values.size()==3) {
                std::stringstream convert(values[2]);
                if( !(convert >> minY) )
                    minY = -1;
            }
            else if(values.size()==4) {
                std::stringstream convert(values[3]);
                if( !(convert >> maxY) )
                    maxY = -1;
            }
            line = line.substr(pos + 1);
        }
        if(whichPred == kCentral)
            hFONLL->SetBinContent(iBin, centY);
        else if(whichPred == kMin)
            hFONLL->SetBinContent(iBin, minY);
        else if(whichPred == kMax)
            hFONLL->SetBinContent(iBin, maxY);
        iBin++;
    }

    inSet.close();

    return hFONLL;
}


//__________________________________________________________________________________________________
template <typename T>
int findPtBin(std::vector<T> binMins, std::vector<T> binMaxs, T value)
{
    if (value < binMins.front())
        return 0;

    if (value >= binMaxs.back()) 
        return binMaxs.size()-1;

    return std::distance(binMins.begin(), std::upper_bound(binMins.begin(), binMins.end(), value)) - 1;
}