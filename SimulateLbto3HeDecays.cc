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
#include "Pythia8Plugins/HepMC3.h"

#endif

using namespace Pythia8;

const int kPdgHe3 = 1000020030;
const double massHe3 = 2.80839160743;

enum FONLLPred {
    kCentral, 
    kMin,
    kMax
};

//__________________________________________________________________________________________________
void SimulateLbto3HeDecays(std::string cfgFileName, int nEvents=1000000, std::string outFileNameRoot="AnalysisResults.root", std::string outFileNameHepMC="AnalysisResults.hepmc", int seed=42);
bool SimpleCoalescence(double p1[3], double p2[3], double p3[3], double pHe3[3], double pc = 0.2 /* GeV/c */);
TH1F* ReadFONLLVsPt(std::string txtFileName, int whichPred, double ptMin, double binWidth, int nBins);

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
    std::string FONLLFileName = config["FONLL"]["filename"].as<std::string>();
    double coalRadius = config["coalescence"]["radius"].as<double>();
    double coalMomRadius = config["coalescence"]["momentum_radius"].as<double>();

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
    // pythia.readString("5122:onIfMatch = 2 1 4 2101"); // bRatio="0.4411147"
    // pythia.readString("5122:onIfMatch = 2 4 1 2101"); // bRatio="0.0910000"
    // pythia.readString("5122:onIfMatch = 2 3 2 2101"); // bRatio="0.0120000"
    // pythia.readString("5122:onIfMatch = 2 3 4 2101"); // bRatio="0.0800000"

    // add 3He
    pythia.particleData.addParticle(kPdgHe3, "3He++", "3He--", 2, 6, 0, massHe3, 0., massHe3, massHe3, 1.e9);   

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.init();

    //__________________________________________________________
    // perform the simulation

    // define HepMC output
    HepMC3::WriterAscii outFileHepMC(outFileNameHepMC);
    HepMC3::Pythia8ToHepMC3 ToHepMC;

    // define histograms
    TH1F* hBR = new TH1F("hBR", ";;BR", 2, 0.5, 2.5);
    hBR->GetXaxis()->SetBinLabel(1, "#Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}");
    hBR->GetXaxis()->SetBinLabel(2, "#Lambda_{b}^{0} #rightarrow (>=2)p + (>=1)n / #Lambda_{b}^{0} #rightarrow d#bar{u}d (ud)_{0}");
    hBR->SetBinContent(1, 0.0120000);

    TH1F* hFONLLLb = ReadFONLLVsPt(FONLLFileName, kCentral, 0.025, 0.05, 2001);
    hFONLLLb->Scale(0.816); // f(b->B) from e+e- provides good normalisation for LHCb and CMS B-meson measurements
    hFONLLLb->SetNameTitle("hFONLLLb_y1", ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (pb GeV^{-1} #it{c})");

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

    for(int iPt=1; iPt<hFONLLLb->GetNbinsX()+1; iPt++)
    {
        double ptCent = hFONLLLb->GetBinCenter(iPt);
        hFONLLLb->SetBinContent(iPt, hFONLLLb->GetBinContent(iPt) * fFFLHCb->Eval(ptCent>5 ? ptCent:5));
    }

    TH1F* hHe3FromLb = (TH1F*)hFONLLLb->Clone("hHe3FromLb_y05");
    hHe3FromLb->Reset();

    TH2F* hPtLbVsHe3FromLb = new TH2F("hPtLbVsHe3FromLb", ";#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(^{3}He) (GeV/#it{c})", 1000, 0., 100., 500., 0., 50.);

    std::vector<double> pxDau, pyDau, pzDau, ptDau;
    std::vector<int> pdgDau, labDau;
    double pxLb, pyLb, pzLb, ptLb;

    int pdgHb = 5122;
    double massLb = TDatabasePDG::Instance()->GetParticle(pdgHb)->Mass();
    int nEventSel = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
        ptLb = hFONLLLb->GetRandom();
        double phiLb = gRandom->Rndm() * 2 * TMath::Pi();
        double yLb = gRandom->Rndm() * 2. - 1.; // flat in -1<y<1
        pxLb = ptLb * TMath::Cos(phiLb);
        pyLb = ptLb * TMath::Sin(phiLb);
        double mt = TMath::Sqrt(massLb * massLb + ptLb * ptLb);
        pzLb = mt * TMath::SinH(yLb);
        double pLb = TMath::Sqrt(ptLb * ptLb + pzLb * pzLb);
        double ELb = TMath::Sqrt(massLb * massLb + pLb * pLb);

        // Lb
        Particle Hb;
        Hb.id(pdgHb);
        Hb.status(11);
        Hb.m(massLb);
        Hb.xProd(0.);
        Hb.yProd(0.);
        Hb.zProd(0.);
        Hb.e(ELb);
        Hb.px(pxLb);
        Hb.py(pyLb);
        Hb.pz(pzLb);
        Hb.tau(4.41000e-01); // PDG2020

        pythia.event.reset();
        pythia.event.append(Hb);
        int idPart = pythia.event[0].id();
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
            }
            else 
            {
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
                double yHe3 = TMath::Log((EHe3 + pCoalHe3[2]) / (EHe3 - pCoalHe3[2]));
                if(TMath::Abs(yHe3) < 0.5) {
                    hHe3FromLb->Fill(ptHe3);
                    hPtLbVsHe3FromLb->Fill(ptLb, ptHe3);
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

                for(int iLab=0; iLab<3; iLab++)
                    pythia.event.remove(labCoal[iLab]-iLab, labCoal[iLab]-iLab); // labels shift by 1 when removing a particle
                pythia.event.append(He3);

                // write HepMC
                HepMC3::GenEvent hepmcevt;
                ToHepMC.fill_next_event(pythia.event, &hepmcevt, -1, &pythia.info);
                outFileHepMC.write_event(hepmcevt);
            }
            nEventSel++;
        }

        ptDau.clear();
        pxDau.clear();
        pyDau.clear();
        pzDau.clear();
        pdgDau.clear();
    }
    outFileHepMC.close();

    hBR->SetBinContent(2, static_cast<double>(nEventSel) / nEvents);
    hHe3FromLb->Scale(hFONLLLb->Integral() / nEventSel * hBR->GetBinContent(1) * hBR->GetBinContent(2));

    // Save histogram on file and close file.
    TFile outFile(outFileNameRoot.data(), "recreate");
    outFile.cd();
    hFONLLLb->Write();
    hHe3FromLb->Write();
    hBR->Write();
    outFile.Close();
}

//__________________________________________________________________________________________________
TH1F* ReadFONLLVsPt(std::string txtFileName, int whichPred, double ptMin, double binWidth, int nBins)
{
    if(txtFileName.find("txt") == string::npos && txtFileName.find("dat") == string::npos && txtFileName.find("csv") == string::npos) {
        std::cerr << "ERROR: Wrong file format! Exit." << std::endl;
        return nullptr;
    }

    std::ifstream inSet(txtFileName.data());
    if(!inSet) {
        std::cerr << "ERROR: Please check if "<< txtFileName.data() << " is the right path. Exit." << std::endl;
        return nullptr;
    }

    std::vector<std::string> values;
    
    TH1F* hFONLL = new TH1F("hFONLL", "", nBins, ptMin, ptMin + nBins*binWidth);
    double centY = -1., minY = -1., maxY = -1.;
    int iPt = 1;
    while(!inSet.eof())
    {    
        std::vector<std::string> values;
        std::string line;
        std::getline(inSet, line);
        if(line.find("#") != std::string::npos || line.find("pt") != std::string::npos)
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
            hFONLL->SetBinContent(iPt, centY);
        else if(whichPred == kMin)
            hFONLL->SetBinContent(iPt, minY);
        else if(whichPred == kMax)
            hFONLL->SetBinContent(iPt, maxY);
        iPt++;
    }

    inSet.close();

    return hFONLL;
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
