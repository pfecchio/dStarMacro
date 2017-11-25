#include <stdlib.h>
#include <iostream>
#include <climits>
#include <map>
#include <exception>

#include <Rtypes.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TAxis.h>
#include <TList.h>
#include <TCanvas.h>
#include <TMath.h>

#include <TLorentzVector.h>

#include <AliAnalysisTaskdStar.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"


// const float kElectronMass = 5.10998909999999971e-04;
// const float kPionMass = 1.39569997787475586e-01f;
// const float kKaonMass = 4.93676990270614624e-01;
// const float kProtonMass = 9.38271999359130859e-01f;
// const float kDeuteronMass = 1.87561297416687012f;
// const float kH3Mass = 2.80925011634826660f;
// const float kHe3Mass = 2.80923008918762207f;
// const float kHe4Mass = 3.72737908363342285f;


using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t ;


void RAnalysis(string mode);
void MCAnalysis(string mode);

float KinemFunction(float x, float y, float z);
float KinemFunctionSquare(float x, float y, float z);
double UpperMDPLimit(double x);
double LowerMDPLimit(double x);

// void PrintKLimit();


void TreeAnalysis(string mode="local"){

  RAnalysis(mode);
  // MCAnalysis(mode);

}


//==============================================================================



void RAnalysis(string mode) {

  TFile *f_input      = nullptr;
  TFile *f_output     = nullptr;

  cout << "trying to open file 'dStarTree.root' " << endl;
  // f_input  = new TFile(Form("~/workspace/dStar/dStarAnalysis/analysisresults/%s_analysis/latest/dStarTree.root",mode.data()), "READ");
    f_input  = new TFile(Form("~/workspace/dStar/dStarAnalysis/analysisresults/grid_analysis/20171123_2035/dStarTree.root",mode.data()), "READ");
  if (!f_input || f_input->IsZombie()) {
    cout << "Error opening file! dStarTree.root NOT FOUND!" << endl;
    return;
  }
  cout << "file opened successfully!" << endl;
  f_output = new TFile(Form("~/workspace/dStar/dStarAnalysis/output/out_%s/TreeAnalysis.root",mode.data()),"RECREATE");

  // varie cose per leggere il Tree
  vector<daughter_struct> *deuteron_vector  = new vector<daughter_struct>;
  vector<daughter_struct> *pion_vector      = new vector<daughter_struct>;

  TTree     *tree        = static_cast<TTree*>(f_input->Get("dStarTree"));
  TBranch   *deuteron    = tree->GetBranch("Deuteron");
  TBranch   *pion        = tree->GetBranch("Pion");

  deuteron->SetAddress(&deuteron_vector);
  pion->SetAddress(&pion_vector);

  //______________________________________________________________________________

  // SETTINGS
  const float fDalitPlotMassCutMin = 2.340;
  const float fDalitPlotMassCutMax = 2.420;

  //______________________________________________________________________________

  // creo una mappa con pdg della madre e frequenza del pdg
  map<int, int> deu_pdg_freq;
  map<int, int> pi_pdg_freq;
  vector<int> *deu_mother_pdg = new vector<int>;
  vector<int> *pi_mother_pdg  = new vector<int>;

  // istogramma con pdg e frequenza
  TH1F *deuMotherPdgFreq;
  TH1F *piMotherPdgFreq;

  // istogramma con i mother_id dei deutoni per controllo
  TH1F *deuMotherId;
  deuMotherId = new TH1F("deu_mother_id", "deuteron mother_id ;mother_id", 500, 0, 500);
  deuMotherId->SetDirectory(0);

  // istogramma con carica dei deutoni per controllo
  TH1F *deuCharge;
  deuCharge = new TH1F("deu_charge", "deu charge", 4, -2, 2);
  deuCharge->SetDirectory(0);

  //______________________________________________________________________________

  // istogrammi con spettri in pT di deutoni e pioni
  TH1F  *deuPt = new TH1F("deu_pT", "Deuteron #it{p}_{T} ;#it{p}_{T} #it{Gev/c}", 100, 0, 10.0);
  deuPt->SetDirectory(0);
  TH1F  *piPt  = new TH1F("pion_pT", "Pion #it{p}_{T} ;#it{p}_{T} #it{Gev/c}", 100, 0, 10.0);
  piPt->SetDirectory(0);

  //______________________________________________________________________________

  // istogramma massa invariate dStar
  TH2F *hMInv;

  // istogrammi massa invariante e pT
  TH2F *hMInvPipPimRhoSameMother;
  TH2F *hMInvPipPimRho;
  TH2F *hMInvPipPimdStar;
  TH2F *hMInvPipPimdStarSameMother;
  TH2F *hMInvPipPimAll;

  // plot massa invariante (projection no pT)
  TH1F *hMInvPipPimRhoSameMotherProj;
  TH1F *hMInvPipPimRhoProj;
  TH1F *hMInvPipPimdStarProj;
  TH1F *hMInvPipPimdStarSameMotherProj;
  TH1F *hMInvPipPimAllProj;

  // Dalitz Plot
  TH2F *hDalitzPlot;

  //______________________________________________________________________________

  // inizializzo gli istogrammi

  // istogramma massa invariante dStar
  hMInv = new TH2F("Minv_dStar", "titolo", 2950, 2.100, 8.000, 50, 0.0, 10);
  hMInv->SetDirectory(0);

  // istogramma massa invariante e pT
  hMInvPipPimRhoSameMother = new TH2F("Minv_pprhosamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.0, 50, 0.0, 10);
  hMInvPipPimRhoSameMother->SetDirectory(0);

  hMInvPipPimRho = new TH2F("Minv_pprho", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimRho->SetDirectory(0);

  hMInvPipPimdStar = new TH2F("Minv_ppdStar", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimdStar->SetDirectory(0);

  hMInvPipPimdStarSameMother = new TH2F("Minv_ppdStarsamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimdStarSameMother->SetDirectory(0);

  hMInvPipPimAll = new TH2F("Minv_pp", "#it{M_{inv}}  #pi^{+} #pi^{-}  ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimAll->SetDirectory(0);

  // plot massa invariante (no pT)
  hMInvPipPimRhoSameMotherProj = new TH1F("Minv_pprhosamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.0);
  hMInvPipPimRhoSameMotherProj->SetDirectory(0);

  hMInvPipPimRhoProj = new TH1F("Minv_pprho_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimRhoProj->SetDirectory(0);

  hMInvPipPimdStarProj = new TH1F("Minv_ppdStar_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimdStarProj->SetDirectory(0);

  hMInvPipPimdStarSameMotherProj = new TH1F("Minv_ppdStarsamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimdStarSameMotherProj->SetDirectory(0);

  hMInvPipPimAllProj = new TH1F("Minv_pp_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimAllProj->SetDirectory(0);

  // Dalitz Plot
  hDalitzPlot = new TH2F("dalitzplot_dStar", "Dalitz Plot d*(2380) #rightarrow d #pi^{+} #pi^{-} ;#it{M_{inv}}  #pi^{+} #pi^{-}; #it{M_{inv}} d #pi^{-}", 190, 0.05, 1.0, 800, 4.0, 8.0);
  hDalitzPlot->SetDirectory(0);

  //==============================================================================
  // Loop sulle entries

  for (int i=0; i < tree->GetEntries(); i++) {

    deuteron_vector->clear();
    pion_vector->clear();
    tree->GetEntry(i);

    //______________________________________________________________________________

    // istogrammi di controllo (carica deu, mother_pdg, spettri pT)
    for (const auto &deu : *deuteron_vector) {
      FourVector_t  deuvec = deu.vec;
      const unsigned char prop = deu.properties;
      // if (prop & p) deuPt->Fill(deuvec.Pt());
      deuPt->Fill(deuvec.Pt());
      deuMotherId->Fill(deu.mother_id);
      if (deu.mother_pdg == 900010020) {
        unsigned char prop = deu.properties;
        (prop & c) ? deuCharge->Fill(0.5) : deuCharge->Fill(-0.5);
      }
      int pdg = deu.mother_pdg;
      auto it1 = std::find(deu_mother_pdg->begin(), deu_mother_pdg->end(), pdg);
      if (it1 == deu_mother_pdg->end()) {
        deu_pdg_freq[pdg] = 1;
        deu_mother_pdg->push_back(pdg);
      } else {
        deu_pdg_freq[pdg] += 1;
      }
    }

    // loop sul vettore contenente i pioni (sia + che -)
    for (const auto &pi : *pion_vector) {
      FourVector_t pivec = pi.vec;
      piPt->Fill(pivec.Pt());

      int pdg = pi.mother_pdg;
      auto it2 = std::find(pi_mother_pdg->begin(), pi_mother_pdg->end(), pdg);
      if (it2 == pi_mother_pdg->end()) {
        pi_pdg_freq[pdg] = 1;
        pi_mother_pdg->push_back(pdg);
      } else {
        pi_pdg_freq[pdg] += 1;
      }
    }

    //______________________________________________________________________________

    // massa invariante π+ π- da ρ imponendo stessa madre in un caso e senza imporlo nell altro
    for (const auto &pip : *pion_vector) {
      // FourVector_t piplus, piminus;
      unsigned char pip_prop = pip.properties;
      if ((pip_prop & c) && (pip.mother_pdg == 113)) {
        for (const auto &pim : *pion_vector) {
          unsigned char pim_prop = pim.properties;
          if (!(pim_prop & c) && (pim.mother_pdg == 113)) {
            FourVector_t piplus = pip.vec;
            FourVector_t piminus = pim.vec;
            piminus += piplus;
            hMInvPipPimRho->Fill(piminus.M(), piminus.Pt());
            if (pim.mother_id == pip.mother_id) hMInvPipPimRhoSameMother->Fill(piminus.M(), piminus.Pt());
          }
        }
      }
    }

    // // massa invariante π+ π- da dStar imponendo stessa madre in un caso e senza imporlo nell altro
    for (const auto &pip : *pion_vector) {
      // FourVector_t piplus, piminus;
      unsigned char pip_prop = pip.properties;
      if ((pip_prop & c) && (pip.mother_pdg == 900010020)) {
        for (const auto &pim : *pion_vector) {
          unsigned char pim_prop = pim.properties;
          if (!(pim_prop & c) && (pim.mother_pdg == 900010020)) {
            FourVector_t piplus = pip.vec;
            FourVector_t piminus = pim.vec;
            piminus += piplus;
            hMInvPipPimdStar->Fill(piminus.M(), piminus.Pt());
            if (pim.mother_id == pip.mother_id) hMInvPipPimdStarSameMother->Fill(piminus.M(), piminus.Pt());
          }
        }
      }
    }

    // // massa invariante tutti π+ π-
    for (const auto &pip : *pion_vector) {
      // FourVector_t piplus, piminus;
      unsigned char pip_prop = pip.properties;
      if ((pip_prop & c) && (pip.mother_pdg != 900010020)) {
        for (const auto &pim : *pion_vector) {
          unsigned char pim_prop = pim.properties;
          if (!(pim_prop & c) && (pim.mother_pdg != 900010020)) {
            FourVector_t piplus = pip.vec;
            FourVector_t piminus = pim.vec;
            piminus += piplus;
            hMInvPipPimAll->Fill(piminus.M(), piminus.Pt());
          }
        }
      }
    }

    //______________________________________________________________________________

    // Dalitz Plot
    for (const auto &deu : *deuteron_vector) {
      const unsigned char pdeu = deu.properties;
      if (deu.mother_pdg != 900010020 || !(pdeu & c)) continue;
      for (const auto& pim : *pion_vector) {
        unsigned char ppim = pim.properties;
        if (pim.mother_id != deu.mother_id || (ppim & c)) continue;
        for (const auto& pip : *pion_vector) {
          unsigned char ppip = pip.properties;
          if ((pip.mother_id != deu.mother_id) || !(ppip & c)) continue;
          FourVector_t deuvec = deu.vec;
          FourVector_t pimvec = pim.vec;
          FourVector_t pipvec = pip.vec;
          pipvec += pimvec;
          pimvec += deuvec;
          deuvec += pipvec;
          if (deuvec.M() < fDalitPlotMassCutMin || deuvec.M() > fDalitPlotMassCutMax) continue;
          hDalitzPlot->Fill(pipvec.M2(), pimvec.M2());
        }
      }
    }

    //______________________________________________________________________________

    // creo l istogramma della massa invariante
    for (const auto &deu : *deuteron_vector) {
      FourVector_t deuvec = deu.vec;
      for (const auto &pip : *pion_vector) {
        unsigned char ppip = pip.properties;
        if (!(ppip & c)) continue;
        FourVector_t pipvec = pip.vec;
        for (const auto &pim : *pion_vector) {
          unsigned char ppim = pim.properties;
          if (ppim & c) continue;
          FourVector_t pimvec = pim.vec;
          pimvec += pipvec;
          pimvec += deuvec;
          hMInv->Fill(pimvec.M(), pimvec.Pt());
        }
      }
    }

  }

  // fine loop sulle entries
  //==============================================================================


  // creo istogrammi con motherPDG-freq
  deuMotherPdgFreq = new TH1F("deuteron's mother freq (pdg)", "titolo", deu_mother_pdg->size(), 0, deu_mother_pdg->size());
  deuMotherPdgFreq->SetDirectory(0);
  piMotherPdgFreq = new TH1F("pion's mother freq (pdg)", "titolo", pi_mother_pdg->size(), 0, pi_mother_pdg->size());
  piMotherPdgFreq->SetDirectory(0);
  sort(deu_mother_pdg->begin(), deu_mother_pdg->end());
  sort(pi_mother_pdg->begin(), pi_mother_pdg->end());

  for(unsigned i=0; i<deu_mother_pdg->size(); i++){
    int d_pdg = deu_mother_pdg->at(i);
    deuMotherPdgFreq->SetBinContent(i+1, deu_pdg_freq[d_pdg]);
    deuMotherPdgFreq->GetXaxis()->SetBinLabel(i+1, Form("%i", d_pdg));
  }

  for(unsigned j=0; j<deu_mother_pdg->size(); j++){
    int p_pdg = pi_mother_pdg->at(j);
    piMotherPdgFreq->SetBinContent(j+1, pi_pdg_freq[p_pdg]);
    piMotherPdgFreq->GetXaxis()->SetBinLabel(j+1, Form("%i", p_pdg));
  }

  deuMotherPdgFreq->Write();
  piMotherPdgFreq->Write();

  //______________________________________________________________________________

  // proietto i TH2 per avere gli isto solo Minv e scrivo tutto
  hMInvPipPimRhoSameMotherProj = (TH1F*)hMInvPipPimRhoSameMother->ProjectionX();
  hMInvPipPimRhoSameMotherProj->SetDirectory(0);
  hMInvPipPimRhoProj = (TH1F*)hMInvPipPimRho->ProjectionX();
  hMInvPipPimRhoProj->SetDirectory(0);
  hMInvPipPimdStarProj = (TH1F*)hMInvPipPimdStar->ProjectionX();
  hMInvPipPimdStarProj->SetDirectory(0);
  hMInvPipPimdStarSameMotherProj = (TH1F*)hMInvPipPimdStarSameMother->ProjectionX();
  hMInvPipPimdStarSameMotherProj->SetDirectory(0);
  hMInvPipPimAllProj = (TH1F*)hMInvPipPimAll->ProjectionX();
  hMInvPipPimAllProj->SetDirectory(0);

  //______________________________________________________________________________

  // definisco le funzioni per i limiti cinematici del dalitz plot
  // TF2 *fKLimit1 = new TF2("kl1", "");
  // TF2 *fKLimit2 = new TF2("kl2", "",);

  //______________________________________________________________________________

  deuPt->Write();
  piPt->Write();
  deuMotherId->Write();
  deuCharge->Write();
  hMInvPipPimRhoSameMother->Write();
  hMInvPipPimRho->Write();
  hMInvPipPimdStar->Write();
  hMInvPipPimdStarSameMother->Write();
  hMInvPipPimAll->Write();
  hMInvPipPimRhoSameMotherProj->Write();
  hMInvPipPimRhoProj->Write();
  hMInvPipPimdStarProj->Write();
  hMInvPipPimdStarSameMotherProj->Write();
  hMInvPipPimAllProj->Write();
  hDalitzPlot->Write();
  hMInv->Write();

  //______________________________________________________________________________

  f_output->Write();
  f_output->Close();
  f_input->Close();

}


//______________________________________________________________________________



void MCAnalysis(string mode) {

  TFile *f_MCinput      = nullptr;
  TFile *f_MCoutput     = nullptr;

  cout << "trying to open file 'dStarMCTree.root' " << endl;
  f_MCinput  = new TFile(Form("~/workspace/dStar/dStarAnalysis/analysisresults/%s_analysis/latest/dStarMCTree.root",mode.data()), "READ");
  if (!f_MCinput || f_MCinput->IsZombie()) {
    cout << "Error opening file! dStarMCTree.root NOT FOUND!" << endl;
    return;
  }
  cout << "file opened successfully!" << endl;
  f_MCoutput = new TFile(Form("~/workspace/dStar/dStarAnalysis/output/out_%s/MCTreeAnalysis.root",mode.data()),"RECREATE");

  // varie cose per leggere il Tree
  vector<daughter_struct> *deuteron_vector  = new vector<daughter_struct>;
  vector<daughter_struct> *pion_vector      = new vector<daughter_struct>;

  TTree     *tree        = static_cast<TTree*>(f_MCinput->Get("dStarMCTree"));
  TBranch   *deuteron    = tree->GetBranch("MCDeuteron");
  TBranch   *pion        = tree->GetBranch("MCPion");

  deuteron->SetAddress(&deuteron_vector);
  pion->SetAddress(&pion_vector);

  //______________________________________________________________________________

  // istogramma con carica dei deutoni per controllo
  TH1F *deuCharge = new TH1F("deu_charge", "deu charge", 2, -1, 1);
  deuCharge->SetDirectory(0);

  //______________________________________________________________________________

  // istogrammi con spettri in pT di deutoni e pioni
  TH1F  *deuPt = new TH1F("deu_pT", "Deuteron #it{p}_{T} ;#it{p}_{T} #it{Gev/c}", 20, 0, 10.0);
  deuPt->SetDirectory(0);
  TH1F  *piPt  = new TH1F("pion_pT", "Pion #it{p}_{T} ;#it{p}_{T} #it{Gev/c}", 20, 0, 10.0);
  piPt->SetDirectory(0);

  //______________________________________________________________________________

  // istogramma massa invariate dStar
  TH2F *hMInv;

  // istogrammi massa invariante e pT
  TH2F *hMInvPipPimRhoSameMother;
  TH2F *hMInvPipPimRho;
  TH2F *hMInvPipPimdStar;
  TH2F *hMInvPipPimdStarSameMother;
  TH2F *hMInvPipPimAll;

  // plot massa invariante (projection no pT)
  TH1F *hMInvPipPimRhoSameMotherProj;
  TH1F *hMInvPipPimRhoProj;
  TH1F *hMInvPipPimdStarProj;
  TH1F *hMInvPipPimdStarSameMotherProj;
  TH1F *hMInvPipPimAllProj;

  // Dalitz Plot
  TH2F *hDalitzPlot;

  //______________________________________________________________________________

  // inizializzo gli istogrammi

  // istogramma massa invariante dStar
  hMInv = new TH2F("Minv_dStar", "titolo", 2950, 2.100, 8.000, 50, 0.0, 10);
  hMInv->SetDirectory(0);

  // istogramma massa invariante e pT
  hMInvPipPimRhoSameMother = new TH2F("Minv_pprhosamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.0, 50, 0.0, 10);
  hMInvPipPimRhoSameMother->SetDirectory(0);

  hMInvPipPimRho = new TH2F("Minv_pprho", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimRho->SetDirectory(0);

  hMInvPipPimdStar = new TH2F("Minv_ppdStar", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimdStar->SetDirectory(0);

  hMInvPipPimdStarSameMother = new TH2F("Minv_ppdStarsamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimdStarSameMother->SetDirectory(0);

  hMInvPipPimAll = new TH2F("Minv_pp", "#it{M_{inv}}  #pi^{+} #pi^{-}  ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  hMInvPipPimAll->SetDirectory(0);

  // plot massa invariante (no pT)
  hMInvPipPimRhoSameMotherProj = new TH1F("Minv_pprhosamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.0);
  hMInvPipPimRhoSameMotherProj->SetDirectory(0);

  hMInvPipPimRhoProj = new TH1F("Minv_pprho_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimRhoProj->SetDirectory(0);

  hMInvPipPimdStarProj = new TH1F("Minv_ppdStar_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimdStarProj->SetDirectory(0);

  hMInvPipPimdStarSameMotherProj = new TH1F("Minv_ppdStarsamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimdStarSameMotherProj->SetDirectory(0);

  hMInvPipPimAllProj = new TH1F("Minv_pp_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  hMInvPipPimAllProj->SetDirectory(0);

  // Dalitz Plot
  hDalitzPlot = new TH2F("dalitzplot_dStar", "Dalitz Plot d*(2380) #rightarrow d #pi^{+} #pi^{-} ;#it{M_{inv}}  #pi^{+} #pi^{-}; #it{M_{inv}} d #pi^{-}", 190, 0.05, 1.0, 800, 4.0, 8.0);
  hDalitzPlot->SetDirectory(0);

  //==============================================================================
  // Loop sulle entries

  for (int i=0; i < tree->GetEntries(); i++) {

    deuteron_vector->clear();
    pion_vector->clear();
    tree->GetEntry(i);

    //______________________________________________________________________________

    // istogrammi di controllo (carica deu, mother_pdg, spettri pT)
    for (const auto &deu : *deuteron_vector) {
      FourVector_t  deuvec = deu.vec;
      deuPt->Fill(deuvec.Pt());
      if (deu.mother_pdg == 900010020) {
        unsigned char prop = deu.properties;
        (prop & c) ? deuCharge->Fill(0.5) : deuCharge->Fill(-0.5);
      }
    }

    // loop sul vettore contenente i pioni (sia + che -)
    for (const auto &pi : *pion_vector) {
      FourVector_t pivec = pi.vec;
      piPt->Fill(pivec.Pt());
    }

    //______________________________________________________________________________

    // // massa invariante π+ π- da ρ imponendo stessa madre in un caso e senza imporlo nell altro
    // for (const auto &pip : *pion_vector) {
    //   FourVector_t piplus, piminus;
    //   unsigned char pip_prop = pip.properties;
    //   if ((pip_prop & c) && (pip.mother_pdg == 113)) {
    //     piplus = pip.vec;
    //     for (const auto &pim : *pion_vector) {
    //       unsigned char pim_prop = pim.properties;
    //       if (!(pim_prop & c) && (pim.mother_pdg == 113)) {
    //         piminus = pim.vec;
    //         piminus += piplus;
    //         hMInvPipPimRho->Fill(piminus.M(), piminus.Pt());
    //         if (pim.mother_id == pip.mother_id) hMInvPipPimRhoSameMother->Fill(piminus.M(), piminus.Pt());
    //       }
    //     }
    //   }
    // }
    //
    // // massa invariante π+ π- da dStar imponendo stessa madre in un caso e senza imporlo nell altro
    // for (const auto &pip : *pion_vector) {
    //   FourVector_t piplus, piminus;
    //   unsigned char pip_prop = pip.properties;
    //   if ((pip_prop & c) && (pip.mother_pdg == 900010020)) {
    //     piplus = pip.vec;
    //     for (const auto &pim : *pion_vector) {
    //       unsigned char pim_prop = pim.properties;
    //       if (!(pim_prop & c) && (pim.mother_pdg == 900010020)) {
    //         piminus = pim.vec;
    //         piminus += piplus;
    //         hMInvPipPimdStar->Fill(piminus.M(), piminus.Pt());
    //         if (pim.mother_id == pip.mother_id) hMInvPipPimdStarSameMother->Fill(piminus.M(), piminus.Pt());
    //       }
    //     }
    //   }
    // }
    //
    // // massa invariante tutti π+ π-
    // for (const auto &pip : *pion_vector) {
    //   FourVector_t piplus, piminus;
    //   unsigned char pip_prop = pip.properties;
    //   if ((pip_prop & c) && (pip.mother_pdg != 900010020)) {
    //     piplus = pip.vec;
    //     for (const auto &pim : *pion_vector) {
    //       unsigned char pim_prop = pim.properties;
    //       if (!(pim_prop & c) && (pim.mother_pdg != 900010020)) {
    //         piminus = pim.vec;
    //         piminus += piplus;
    //         hMInvPipPimAll->Fill(piminus.M(), piminus.Pt());
    //       }
    //     }
    //   }
    // }

    //______________________________________________________________________________

    // Dalitz Plot non necessario, lo faccio già nel Task
    // for (const auto &deu : *deuteron_vector) {
    //   const unsigned char deu_prop = deu.properties;
    //   if (deu.mother_pdg != 900010020 || !(deu_prop & c)) continue;
    //   for (const auto &pim : *pion_vector) {
    //     if (pim.mother_id != deu.mother_id) continue;
    //     const unsigned char pim_prop = pim.properties;
    //     if (pim_prop & c) continue;
    //     for (const auto &pip : *pion_vector) {
    //       if (pip.mother_id != deu.mother_id) continue;
    //       const unsigned char pip_prop = pip.properties;
    //       if (!(pip_prop & c)) continue;
    //       FourVector_t deu_vec = deu.vec;
    //       FourVector_t pim_vec = pim.vec;
    //       FourVector_t pip_vec = pip.vec;
    //       pip_vec += pim_vec;
    //       pim_vec += deu_vec;
    //       hDalitzPlot->Fill(pip_vec.M2(), pim_vec.M2());
    //     }
    //   }
    // }

    //______________________________________________________________________________

    // // creo l istogramma della massa invariante
    // FourVector_t v1, v2, v3;
    //
    // for (auto &deu : *deuteron_vector) {
    //   v1 = deu.vec;
    //   for (auto &pip : *pion_vector) {
    //     unsigned char pip_prop = pip.properties;
    //     if (pip_prop & c) {
    //       v2 = pip.vec;
    //       for (auto &pim : *pion_vector) {
    //         unsigned char pim_prop = pim.properties;
    //         if (!(pim_prop & c)) {
    //           v3 = pim.vec;
    //           v3 += v2;
    //           v3 += v1;
    //           hMInv->Fill(v3.M(), v3.Pt());
    //         }
    //       }
    //     }
    //   }
    // }

  }

  // fine loop sulle entries
  //==============================================================================


  // proietto i TH2 per avere gli isto solo Minv e scrivo tutto
  hMInvPipPimRhoSameMotherProj = (TH1F*)hMInvPipPimRhoSameMother->ProjectionX();
  hMInvPipPimRhoSameMotherProj->SetDirectory(0);
  hMInvPipPimRhoProj = (TH1F*)hMInvPipPimRho->ProjectionX();
  hMInvPipPimRhoProj->SetDirectory(0);
  hMInvPipPimdStarProj = (TH1F*)hMInvPipPimdStar->ProjectionX();
  hMInvPipPimdStarProj->SetDirectory(0);
  hMInvPipPimdStarSameMotherProj = (TH1F*)hMInvPipPimdStarSameMother->ProjectionX();
  hMInvPipPimdStarSameMotherProj->SetDirectory(0);
  hMInvPipPimAllProj = (TH1F*)hMInvPipPimAll->ProjectionX();
  hMInvPipPimAllProj->SetDirectory(0);

  deuPt->Write();
  piPt->Write();
  deuCharge->Write();
  hMInvPipPimRhoSameMother->Write();
  hMInvPipPimRho->Write();
  hMInvPipPimdStar->Write();
  hMInvPipPimdStarSameMother->Write();
  hMInvPipPimAll->Write();
  hMInvPipPimRhoSameMotherProj->Write();
  hMInvPipPimRhoProj->Write();
  hMInvPipPimdStarProj->Write();
  hMInvPipPimdStarSameMotherProj->Write();
  hMInvPipPimAllProj->Write();
  hDalitzPlot->Write();
  // hMInv->Write();

  //______________________________________________________________________________

  f_MCoutput->Write();
  f_MCoutput->Close();
  f_MCinput->Close();

}


//______________________________________________________________________________


float KinemFunction(float x, float y, float z) {return x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;}


//______________________________________________________________________________


float KinemFunctionSquare(float x, float y, float z) { return TMath::Sqrt(x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z);}


//______________________________________________________________________________


double UpperMDPLimit(double s1) {

  // FourVector_t vdpim = vdeu + vpim;
  // FourVector_t vdpip = vdeu + vpip;
  // FourVector_t vpp = vpim + vpip;

  // float s1 = vpp.M2();
  // float s2 = vdpip.M2();
  // float s3 = vdpim.M2();

  const float kPionMass = 1.39569997787475586e-01f;
  const float kDeuteronMass = 1.87561297416687012f;
  const float kDeuteronMass2 = kDeuteronMass * kDeuteronMass;
  const float kPionMass2 = kPionMass * kPionMass;

  float s = 5.856;
  float a = KinemFunctionSquare(s1,s,kDeuteronMass2)*KinemFunctionSquare(s1,kPionMass2,kPionMass2);

  return kDeuteronMass2 + kPionMass2 + (1./(2.*s1))*((s - s1 - kDeuteronMass2)*(s1 - kPionMass2 - kPionMass2) + a);
}


//______________________________________________________________________________


double LowerMDPLimit(double s1) {

  const float kPionMass = 1.39569997787475586e-01f;
  const float kDeuteronMass = 1.87561297416687012f;
  const float kDeuteronMass2 = kDeuteronMass * kDeuteronMass;
  const float kPionMass2 = kPionMass * kPionMass;

  float s = 5.856;

  return kDeuteronMass2 + kPionMass2 + (1./(2.*s1))*((s - s1 - kDeuteronMass2)*(s1 - kPionMass2 - kPionMass2) - KinemFunctionSquare(s1,s,kDeuteronMass2)*KinemFunctionSquare(s1,kPionMass2,kPionMass2));
}


//______________________________________________________________________________


void PrintKLimit() {

  TFile *f = new TFile("func.root","RECREATE");

  TCanvas *fcanvas = new TCanvas("fcanvas","c",500,500);
  TF1 *f1 = new TF1("f1", "UpperMDPLimit(x)", 0.2, 0.22);
  f1->Draw();
  fcanvas->Write();

  f->Write();
  f->Close();

}
