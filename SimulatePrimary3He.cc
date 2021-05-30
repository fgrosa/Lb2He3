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
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#endif

using namespace Pythia8;

const int kPdgHe3 = 1000020030;
const double massHe3 = 2.80839160743;

enum syst 
{
    kpp13TeV = 0
};

//__________________________________________________________________________________________________
void SimulatePrimary3He(int system=kpp13TeV, int nEvents=1000000, std::string outFileNameRoot="AnalysisResults.root", std::string outFileNameHepMC="AnalysisResults.hepmc", int seed=42);

//__________________________________________________________________________________________________
void SimulatePrimary3He(int system, int nEvents, std::string outFileNameRoot, std::string outFileNameHepMC, int seed)
{

    //__________________________________________________________
    // create and configure pythia generator
    Pythia pythia;
    pythia.readString("SoftQCD:all = on"); // not relevant, alla the particles will be replaced by the only He3
    pythia.readString(Form("Tune:pp = 14"));

    // add 3He
    pythia.particleData.addParticle(kPdgHe3, "3He++", "3He--", 2, 6, 0, massHe3, 0., massHe3, massHe3, 1.e9);   

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.init();

    // define HepMC output
    HepMC3::WriterAscii outFileHepMC(outFileNameHepMC);
    HepMC3::Pythia8ToHepMC3 ToHepMC;

    // define input pT shape
    TF1* fPtShape = nullptr;
    double sigmaV0 = 0.;
    if(system==kpp13TeV)
    {
        TFile* inFileSpectrum = TFile::Open("inputs/heliumSpectra_pp13TeV.root");
        fPtShape = (TF1*)inFileSpectrum->Get("fCombineHeliumSpecLevyFit_0-100");
        sigmaV0 = 58.7e9 * fPtShape->Integral(0, 100);
    }
    fPtShape->SetName("fPrimary3He");

    // define output histogram
    TH1F* hPrimary3He = new TH1F("hPrimary3He", ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} (pb GeV^{-1} #it{c})", 2001, 0., 100.05);

    double pxHe3, pyHe3, pzHe3, ptHe3;
    int nEventSel = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {

        if(iEvent%1000000 == 0)
            std::cout << Form("He3 number %10d\r", iEvent);

        ptHe3 = fPtShape->GetRandom(); // to be substitute with TSallis or BlastWave 
        double phiHe3 = gRandom->Rndm() * 2 * TMath::Pi();
        double yHe3 = gRandom->Rndm() - 0.5; // flat in -0.5<y<0.5
        pxHe3 = ptHe3 * TMath::Cos(phiHe3);
        pyHe3 = ptHe3 * TMath::Sin(phiHe3);
        double mt = TMath::Sqrt(massHe3 * massHe3 + ptHe3 * ptHe3);
        pzHe3 = mt * TMath::SinH(yHe3);
        double pHe3 = TMath::Sqrt(ptHe3 * ptHe3 + pzHe3 * pzHe3);
        double EHe3 = TMath::Sqrt(massHe3 * massHe3 + pHe3 * pHe3);

        hPrimary3He->Fill(ptHe3);
    
        // He3
        Particle He3;
        He3.id(kPdgHe3);
        He3.status(93);
        He3.m(massHe3);
        He3.xProd(0.);
        He3.yProd(0.);
        He3.zProd(0.);
        He3.tProd(0.);
        He3.e(EHe3);
        He3.px(pxHe3);
        He3.py(pyHe3);
        He3.pz(pzHe3);
        He3.mother1(1);
        He3.mother2(0);
        He3.daughter1(0);
        He3.daughter2(0);
        He3.tau(1.e9); // stable

        pythia.event.reset();
        pythia.event.append(21, 11, 0, 0, 0, 0, 0, 0, 0., 0., 0., 0., 0.); //add a dummy gluon to make hepmc happy
        pythia.event.append(He3);
        pythia.event.remove(2, pythia.event.size()-2);

        // write HepMC
        HepMC3::GenEvent hepmcevt;
        ToHepMC.fill_next_event(pythia.event, &hepmcevt, -1, &pythia.info);
        outFileHepMC.write_event(hepmcevt);
    }
    outFileHepMC.close();
    hPrimary3He->Scale(sigmaV0/nEvents);

    // Save histogram on file and close file.
    TFile outFile(outFileNameRoot.data(), "recreate");
    outFile.cd();
    fPtShape->Write();
    hPrimary3He->Write();
    outFile.Close();
}
