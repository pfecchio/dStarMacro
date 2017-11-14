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
#include <TAxis.h>
#include <TList.h>
#include <TCanvas.h>
#include <TMath.h>

#include <TLorentzVector.h>

#include <AliAnalysisTaskdStar.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"


using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t ;


void RAnalysis();
void MCAnalysis();


void TreeAnalysis(char selection='B'){

  // Ranalysis();
  MCAnalysis();

}


//==============================================================================



void RAnalysis() {

  TFile *f_input      = nullptr;
  TFile *f_output     = nullptr;

  cout << "trying to open -> dStarTree.root" << endl;
  f_input  = new TFile("~/workspace/dStar/dStarAnalysis/analysisresults/local_analysis/20171109_1750/dStarTree.root", "READ");
  if (!f_input || f_input->IsZombie()) {
    cout << "Error opening file! dStarTree.root NOT FOUND!" << endl;
    return;
  }
  f_output = new TFile("~/workspace/dStar/dStarAnalysis/output/out_local/TreeAnalysis.root","RECREATE");

  // varie cose per leggere il Tree
  vector<daughter_struct> *deuteron_vector  = new vector<daughter_struct>;
  vector<daughter_struct> *pion_vector      = new vector<daughter_struct>;

  TTree     *tree        = static_cast<TTree*>(f_input->Get("dStarTree"));
  TBranch   *deuteron    = tree->GetBranch("Deuteron");
  TBranch   *pion        = tree->GetBranch("Pi");

  deuteron->SetAddress(&deuteron_vector);
  pion->SetAddress(&pion_vector);

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
  TH1F  *deuPt = new TH1F("deu_pT", "Deuteron #it{p}_{T} ;#it{p}_{T} #it{Gev/c}", 20, 0, 10.0);
  deuPt->SetDirectory(0);
  TH1F  *piPt  = new TH1F("pion_pT", "Pion #it{p}_{T} ;#it{p}_{T} #it{Gev/c}", 20, 0, 10.0);
  piPt->SetDirectory(0);

  //______________________________________________________________________________

  // istogramma massa invariate dStar
  TH2F *mInv;

  // istogrammi massa invariante e pT
  TH2F *mInvPipPimRhoSameMother;
  TH2F *mInvPipPimRho;
  TH2F *mInvPipPimdStar;
  TH2F *mInvPipPimdStarSameMother;
  TH2F *mInvPipPimAll;

  // plot massa invariante (projection no pT)
  TH1F *mInvPipPimRhoSameMotherProj;
  TH1F *mInvPipPimRhoProj;
  TH1F *mInvPipPimdStarProj;
  TH1F *mInvPipPimdStarSameMotherProj;
  TH1F *mInvPipPimAllProj;

  // Dalitz Plot
  TH2F *DalitzPlot;

  //______________________________________________________________________________

  // inizializzo gli istogrammi

  // istogramma massa invariante dStar
  mInv = new TH2F("Minv_dStar", "titolo", 2950, 2.100, 8.000, 50, 0.0, 10);
  mInv->SetDirectory(0);

  // istogramma massa invariante e pT
  mInvPipPimRhoSameMother = new TH2F("Minv_pprhosamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.0, 50, 0.0, 10);
  mInvPipPimRhoSameMother->SetDirectory(0);

  mInvPipPimRho = new TH2F("Minv_pprho", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimRho->SetDirectory(0);

  mInvPipPimdStar = new TH2F("Minv_ppdStar", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimdStar->SetDirectory(0);

  mInvPipPimdStarSameMother = new TH2F("Minv_ppdStarsamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimdStarSameMother->SetDirectory(0);

  mInvPipPimAll = new TH2F("Minv_pp", "#it{M_{inv}}  #pi^{+} #pi^{-}  ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimAll->SetDirectory(0);

  // plot massa invariante (no pT)
  mInvPipPimRhoSameMotherProj = new TH1F("Minv_pprhosamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.0);
  mInvPipPimRhoSameMotherProj->SetDirectory(0);

  mInvPipPimRhoProj = new TH1F("Minv_pprho_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimRhoProj->SetDirectory(0);

  mInvPipPimdStarProj = new TH1F("Minv_ppdStar_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimdStarProj->SetDirectory(0);

  mInvPipPimdStarSameMotherProj = new TH1F("Minv_ppdStarsamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimdStarSameMotherProj->SetDirectory(0);

  mInvPipPimAllProj = new TH1F("Minv_pp_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimAllProj->SetDirectory(0);

  // Dalitz Plot
  DalitzPlot = new TH2F("dalitzplot_dStar", "Dalitz Plot d*(2380) #rightarrow d #pi^{+} #pi^{-} ;#it{M_{inv}}  #pi^{+} #pi^{-}; #it{M_{inv}} d #pi^{-}", 2000, 0.0, 20.0, 2000, 0., 20.0);
  DalitzPlot->SetDirectory(0);

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
      deuMotherId->Fill(deu.mother_id);
      if (deu.mother_pdg == 900010020){
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
      FourVector_t piplus, piminus;
      unsigned char pip_prop = pip.properties;
      if ((pip_prop & c) && (pip.mother_pdg == 113)) {
        // cout << "pip ok" << endl;
        piplus = pip.vec;
        for (const auto &pim : *pion_vector) {
          unsigned char pim_prop = pim.properties;
          if (!(pim_prop & c) && (pim.mother_pdg == 113)) {
            // cout << "pim ok" << endl;
            piminus = pim.vec;
            piminus += piplus;
            mInvPipPimRho->Fill(piminus.M(), piminus.Pt());
            if (pim.mother_id == pip.mother_id) mInvPipPimRhoSameMother->Fill(piminus.M(), piminus.Pt());
          }
        }
      }
    }

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
    //         mInvPipPimdStar->Fill(piminus.M(), piminus.Pt());
    //         if (pim.mother_id == pip.mother_id) mInvPipPimdStarSameMother->Fill(piminus.M(), piminus.Pt());
    //       }
    //     }
    //   }
    // }

    // massa invariante tutti π+ π-
    for (const auto &pip : *pion_vector) {
      FourVector_t piplus, piminus;
      unsigned char pip_prop = pip.properties;
      if ((pip_prop & c) && (pip.mother_pdg != 900010020)) {
        piplus = pip.vec;
        for (const auto &pim : *pion_vector) {
          unsigned char pim_prop = pim.properties;
          if (!(pim_prop & c) && (pim.mother_pdg != 900010020)) {
            piminus = pim.vec;
            piminus += piplus;
            mInvPipPimAll->Fill(piminus.M(), piminus.Pt());
          }
        }
      }
    }

    //______________________________________________________________________________

    // Dalitz Plot
    // for (const auto &deu : *deuteron_vector) {
    //   FourVector_t deuteron = {0.f,0.f,0.f,0.f};
    //   unsigned char deu_prop = deu.properties;
    //   if (deu.mother_pdg == 900010020 && (deu_prop & c)) {
    //     deuteron = deu.vec;
    //     for (const auto &pip : *pion_vector) {
    //       FourVector_t piplus = {0.f,0.f,0.f,0.f};
    //       unsigned char pip_prop = pip.properties;
    //       if (pip.mother_id == deu.mother_id && (pip_prop & c)) {
    //         piplus = pip.vec;
    //         for (const auto &pim : *pion_vector) {
    //           FourVector_t piminus = {0.f,0.f,0.f,0.f};
    //           unsigned char pim_prop = pim.properties;
    //           if (pim.mother_id == deu.mother_id && !(pim_prop & c)) {
    //             piminus = pim.vec;
    //             piplus += piminus;
    //             piminus += deuteron;
    //             DalitzPlot->Fill(piplus.M(), piminus.M());
    //           }
    //         }
    //       }
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
    //           mInv->Fill(v3.M(), v3.Pt());
    //         }
    //       }
    //     }
    //   }
    // }

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
  cout << "Fino a Qua Tutto Bene!" << endl;

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
  mInvPipPimRhoSameMotherProj = (TH1F*)mInvPipPimRhoSameMother->ProjectionX();
  mInvPipPimRhoSameMotherProj->SetDirectory(0);
  mInvPipPimRhoProj = (TH1F*)mInvPipPimRho->ProjectionX();
  mInvPipPimRhoProj->SetDirectory(0);
  mInvPipPimdStarProj = (TH1F*)mInvPipPimdStar->ProjectionX();
  mInvPipPimdStarProj->SetDirectory(0);
  mInvPipPimdStarSameMotherProj = (TH1F*)mInvPipPimdStarSameMother->ProjectionX();
  mInvPipPimdStarSameMotherProj->SetDirectory(0);
  mInvPipPimAllProj = (TH1F*)mInvPipPimAll->ProjectionX();
  mInvPipPimAllProj->SetDirectory(0);

  deuPt->Write();
  piPt->Write();
  deuMotherId->Write();
  deuCharge->Write();
  mInvPipPimRhoSameMother->Write();
  mInvPipPimRho->Write();
  mInvPipPimdStar->Write();
  mInvPipPimdStarSameMother->Write();
  mInvPipPimAll->Write();
  mInvPipPimRhoSameMotherProj->Write();
  mInvPipPimRhoProj->Write();
  mInvPipPimdStarProj->Write();
  mInvPipPimdStarSameMotherProj->Write();
  mInvPipPimAllProj->Write();
  DalitzPlot->Write();
  // mInv->Write();

  //______________________________________________________________________________

  f_output->Write();
  f_output->Close();
  f_input->Close();

}


//______________________________________________________________________________



void MCAnalysis() {

  TFile *f_MCinput      = nullptr;
  TFile *f_MCoutput     = nullptr;

  cout << "trying to open -> dStarMCTree.root" << endl;
  f_MCinput  = new TFile("~/workspace/dStar/dStarAnalysis/analysisresults/local_analysis/20171109_1750/dStarMCTree.root", "READ");
  if (!f_MCinput || f_MCinput->IsZombie()) {
    cout << "Error opening file! dStarMCTree.root NOT FOUND!" << endl;
    return;
  }
  f_MCoutput = new TFile("~/workspace/dStar/dStarAnalysis/output/out_local/MCTreeAnalysis.root","RECREATE");

  // varie cose per leggere il Tree
  vector<daughter_struct> *deuteron_vector  = new vector<daughter_struct>;
  vector<daughter_struct> *pion_vector      = new vector<daughter_struct>;

  TTree     *tree        = static_cast<TTree*>(f_MCinput->Get("dStarMCTree"));
  TBranch   *deuteron    = tree->GetBranch("Deuteron");
  TBranch   *pion        = tree->GetBranch("Pi");

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
  TH2F *mInv;

  // istogrammi massa invariante e pT
  TH2F *mInvPipPimRhoSameMother;
  TH2F *mInvPipPimRho;
  TH2F *mInvPipPimdStar;
  TH2F *mInvPipPimdStarSameMother;
  TH2F *mInvPipPimAll;

  // plot massa invariante (projection no pT)
  TH1F *mInvPipPimRhoSameMotherProj;
  TH1F *mInvPipPimRhoProj;
  TH1F *mInvPipPimdStarProj;
  TH1F *mInvPipPimdStarSameMotherProj;
  TH1F *mInvPipPimAllProj;

  // Dalitz Plot
  TH2F *DalitzPlot;

  //______________________________________________________________________________

  // inizializzo gli istogrammi

  // istogramma massa invariante dStar
  mInv = new TH2F("Minv_dStar", "titolo", 2950, 2.100, 8.000, 50, 0.0, 10);
  mInv->SetDirectory(0);

  // istogramma massa invariante e pT
  mInvPipPimRhoSameMother = new TH2F("Minv_pprhosamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.0, 50, 0.0, 10);
  mInvPipPimRhoSameMother->SetDirectory(0);

  mInvPipPimRho = new TH2F("Minv_pprho", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimRho->SetDirectory(0);

  mInvPipPimdStar = new TH2F("Minv_ppdStar", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimdStar->SetDirectory(0);

  mInvPipPimdStarSameMother = new TH2F("Minv_ppdStarsamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimdStarSameMother->SetDirectory(0);

  mInvPipPimAll = new TH2F("Minv_pp", "#it{M_{inv}}  #pi^{+} #pi^{-}  ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
  mInvPipPimAll->SetDirectory(0);

  // plot massa invariante (no pT)
  mInvPipPimRhoSameMotherProj = new TH1F("Minv_pprhosamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.0);
  mInvPipPimRhoSameMotherProj->SetDirectory(0);

  mInvPipPimRhoProj = new TH1F("Minv_pprho_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimRhoProj->SetDirectory(0);

  mInvPipPimdStarProj = new TH1F("Minv_ppdStar_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimdStarProj->SetDirectory(0);

  mInvPipPimdStarSameMotherProj = new TH1F("Minv_ppdStarsamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimdStarSameMotherProj->SetDirectory(0);

  mInvPipPimAllProj = new TH1F("Minv_pp_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
  mInvPipPimAllProj->SetDirectory(0);

  // Dalitz Plot
  DalitzPlot = new TH2F("dalitzplot_dStar", "Dalitz Plot d*(2380) #rightarrow d #pi^{+} #pi^{-} ;#it{M_{inv}}  #pi^{+} #pi^{-}; #it{M_{inv}} d #pi^{-}", 2000, 0.0, 20.0, 2000, 0., 20.0);
  DalitzPlot->SetDirectory(0);

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
    //     // cout << "pip ok" << endl;
    //     piplus = pip.vec;
    //     for (const auto &pim : *pion_vector) {
    //       unsigned char pim_prop = pim.properties;
    //       if (!(pim_prop & c) && (pim.mother_pdg == 113)) {
    //         // cout << "pim ok" << endl;
    //         piminus = pim.vec;
    //         piminus += piplus;
    //         mInvPipPimRho->Fill(piminus.M(), piminus.Pt());
    //         if (pim.mother_id == pip.mother_id) mInvPipPimRhoSameMother->Fill(piminus.M(), piminus.Pt());
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
    //         mInvPipPimdStar->Fill(piminus.M(), piminus.Pt());
    //         if (pim.mother_id == pip.mother_id) mInvPipPimdStarSameMother->Fill(piminus.M(), piminus.Pt());
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
    //         mInvPipPimAll->Fill(piminus.M(), piminus.Pt());
    //       }
    //     }
    //   }
    // }

    //______________________________________________________________________________

    // Dalitz Plot
    // for (const auto &deu : *deuteron_vector) {
    //   FourVector_t deuteron = {0.f,0.f,0.f,0.f};
    //   unsigned char deu_prop = deu.properties;
    //   if (deu.mother_pdg == 900010020 && (deu_prop & c)) {
    //     deuteron = deu.vec;
    //     for (const auto &pip : *pion_vector) {
    //       FourVector_t piplus = {0.f,0.f,0.f,0.f};
    //       unsigned char pip_prop = pip.properties;
    //       if (pip.mother_id == deu.mother_id && (pip_prop & c)) {
    //         piplus = pip.vec;
    //         for (const auto &pim : *pion_vector) {
    //           FourVector_t piminus = {0.f,0.f,0.f,0.f};
    //           unsigned char pim_prop = pim.properties;
    //           if (pim.mother_id == deu.mother_id && !(pim_prop & c)) {
    //             piminus = pim.vec;
    //             piplus += piminus;
    //             piminus += deuteron;
    //             DalitzPlot->Fill(piplus.M(), piminus.M());
    //           }
    //         }
    //       }
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
    //           mInv->Fill(v3.M(), v3.Pt());
    //         }
    //       }
    //     }
    //   }
    // }

  }

  // fine loop sulle entries
  //==============================================================================


  // proietto i TH2 per avere gli isto solo Minv e scrivo tutto
  mInvPipPimRhoSameMotherProj = (TH1F*)mInvPipPimRhoSameMother->ProjectionX();
  mInvPipPimRhoSameMotherProj->SetDirectory(0);
  mInvPipPimRhoProj = (TH1F*)mInvPipPimRho->ProjectionX();
  mInvPipPimRhoProj->SetDirectory(0);
  mInvPipPimdStarProj = (TH1F*)mInvPipPimdStar->ProjectionX();
  mInvPipPimdStarProj->SetDirectory(0);
  mInvPipPimdStarSameMotherProj = (TH1F*)mInvPipPimdStarSameMother->ProjectionX();
  mInvPipPimdStarSameMotherProj->SetDirectory(0);
  mInvPipPimAllProj = (TH1F*)mInvPipPimAll->ProjectionX();
  mInvPipPimAllProj->SetDirectory(0);

  deuPt->Write();
  piPt->Write();
  deuCharge->Write();
  mInvPipPimRhoSameMother->Write();
  mInvPipPimRho->Write();
  mInvPipPimdStar->Write();
  mInvPipPimdStarSameMother->Write();
  mInvPipPimAll->Write();
  mInvPipPimRhoSameMotherProj->Write();
  mInvPipPimRhoProj->Write();
  mInvPipPimdStarProj->Write();
  mInvPipPimdStarSameMotherProj->Write();
  mInvPipPimAllProj->Write();
  DalitzPlot->Write();
  // mInv->Write();

  //______________________________________________________________________________

  f_MCoutput->Write();
  f_MCoutput->Close();
  f_MCinput->Close();

}















//
// // creo una mappa con pdg della madre e frequenza del pdg
// map<int, int> deu_pdg_freq;
// map<int, int> pi_pdg_freq;
//
// vector<int> *deu_mother_pdg = new vector<int>;
// vector<int> *pi_mother_pdg  = new vector<int>;

// // creo un istogramma con pdg e frequenza
// TH1F *deuMotherPdgFreq;
// TH1F *piMotherPdgFreq;

//
// // istogramma massa invariate dStar
// TH2F *mInv;
//
// // istogrammi massa invariante e pT
// TH2F *mInvPipPimRhoSameMother;
// TH2F *mInvPipPimRho;
// TH2F *mInvPipPimdStar;
// TH2F *mInvPipPimdStarSameMother;
// TH2F *mInvPipPimAll;

//
// // plot massa invariante (projection no pT)
// TH1F *mInvPipPimRhoSameMotherProj;
// TH1F *mInvPipPimRhoProj;
// TH1F *mInvPipPimdStarProj;
// TH1F *mInvPipPimdStarSameMotherProj;
// TH1F *mInvPipPimAllProj;

//
// // Dalitz Plot
// TH2F *DalitzPlot;

//
// // istogramma con i mother_id dei deutoni per controllo
// TH1F *deuMotherId;



// // istogramma massa invariante dStar
// mInv = new TH2F("Minv_dStar", "titolo", 2950, 2.100, 8.000, 50, 0.0, 10);
// mInv->SetDirectory(0);
//
// // istogramma massa invariante e pT
// mInvPipPimRhoSameMother = new TH2F("Minv_pprhosamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.0, 50, 0.0, 10);
// mInvPipPimRhoSameMother->SetDirectory(0);
//
// mInvPipPimRho = new TH2F("Minv_pprho", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
// mInvPipPimRho->SetDirectory(0);
//
// mInvPipPimdStar = new TH2F("Minv_ppdStar", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
// mInvPipPimdStar->SetDirectory(0);
//
// mInvPipPimdStarSameMother = new TH2F("Minv_ppdStarsamemother", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
// mInvPipPimdStarSameMother->SetDirectory(0);
//
// mInvPipPimAll = new TH2F("Minv_pp", "#it{M_{inv}}  #pi^{+} #pi^{-}  ;#it{M_{inv}}  #it{GeV/c^{2}}; #it{p}_{T} #it{Gev/c}", 800, 0.0, 8.000, 50, 0.0, 10);
// mInvPipPimAll->SetDirectory(0);
//
// // plot massa invariante (no pT)
// mInvPipPimRhoSameMotherProj = new TH1F("Minv_pprhosamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.0);
// mInvPipPimRhoSameMotherProj->SetDirectory(0);
//
// mInvPipPimRhoProj = new TH1F("Minv_pprho_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from #rho ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
// mInvPipPimRhoProj->SetDirectory(0);
//
// mInvPipPimdStarProj = new TH1F("Minv_ppdStar_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
// mInvPipPimdStarProj->SetDirectory(0);
//
// mInvPipPimdStarSameMotherProj = new TH1F("Minv_ppdStarsamemother_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} from d* (same mother imposed) ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
// mInvPipPimdStarSameMotherProj->SetDirectory(0);
//
// mInvPipPimAllProj = new TH1F("Minv_pp_proj", "#it{M_{inv}}  #pi^{+} #pi^{-} ;#it{M_{inv}}  #it{GeV/c^{2}}", 800, 0.0, 8.000);
// mInvPipPimAllProj->SetDirectory(0);
//
// // Dalitz Plot
// DalitzPlot = new TH2F("dalitzplot_dStar", "Dalitz Plot d*(2380) #rightarrow d #pi^{+} #pi^{-} ;#it{M_{inv}}  #pi^{+} #pi^{-}; #it{M_{inv}} d #pi^{-}", 150, 0.25, 1.0, 600, 2.0, 5.0);
// DalitzPlot->SetDirectory(0);
//
// // istogramma con i mother_id dei deutoni per controllo
// deuMotherId = new TH1F("deu_mother_id", "deuteron mother_id ;mother_id", 500, 0, 500);
// deuMotherId->SetDirectory(0);
//




//______________________________________________________________________________

// for (int i=0; i < tree->GetEntries(); i++) {
//
//   deuteron_vector->clear();
//   pion_vector->clear();
//   tree->GetEntry(i);
//
//   // loop sul vettore contenente i deutoni
//   for (const auto &deu : *deuteron_vector) {
//     FourVector_t  deuvec = deu.vec;
//     deuPt->Fill(deuvec.Pt());
//     deuMotherId->Fill(deu.mother_id);
//     if (deu.mother_pdg == 900010020){
//       unsigned char prop = deu.properties;
//       (prop & c) ? deuCharge->Fill(0.5) : deuCharge->Fill(-0.5);
//     }
//     int pdg = deu.mother_pdg;
//     auto it1 = std::find(deu_mother_pdg->begin(), deu_mother_pdg->end(), pdg);
//     if (it1 == deu_mother_pdg->end()) {
//       deu_pdg_freq[pdg] = 1;
//       deu_mother_pdg->push_back(pdg);
//     } else {
//       deu_pdg_freq[pdg] += 1;
//     }
//   }
//
//   // loop sul vettore contenente i pioni (sia + che -)
//   for (const auto &pi : *pion_vector) {
//     FourVector_t pivec = pi.vec;
//     piPt->Fill(pivec.Pt());
//
//     int pdg = pi.mother_pdg;
//     auto it2 = std::find(pi_mother_pdg->begin(), pi_mother_pdg->end(), pdg);
//     if (it2 == pi_mother_pdg->end()) {
//       pi_pdg_freq[pdg] = 1;
//       pi_mother_pdg->push_back(pdg);
//     } else {
//       pi_pdg_freq[pdg] += 1;
//     }
//   }

//______________________________________________________________________________

// // massa invariante π+ π- da ρ imponendo stessa madre in un caso e senza imporlo nell altro
// for (const auto &pip : *pion_vector) {
//   FourVector_t piplus, piminus;
//   unsigned char pip_prop = pip.properties;
//   if ((pip_prop & c) && (pip.mother_pdg == 113)) {
//     // cout << "pip ok" << endl;
//     piplus = pip.vec;
//     for (const auto &pim : *pion_vector) {
//       unsigned char pim_prop = pim.properties;
//       if (!(pim_prop & c) && (pim.mother_pdg == 113)) {
//         // cout << "pim ok" << endl;
//         piminus = pim.vec;
//         piminus += piplus;
//         mInvPipPimRho->Fill(piminus.M(), piminus.Pt());
//         if (pim.mother_id == pip.mother_id) mInvPipPimRhoSameMother->Fill(piminus.M(), piminus.Pt());
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
//         mInvPipPimdStar->Fill(piminus.M(), piminus.Pt());
//         if (pim.mother_id == pip.mother_id) mInvPipPimdStarSameMother->Fill(piminus.M(), piminus.Pt());
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
//         mInvPipPimAll->Fill(piminus.M(), piminus.Pt());
//       }
//     }
//   }
// }
//
// //______________________________________________________________________________
//
// // Dalitz Plot
// for (const auto &deu : *deuteron_vector) {
//   FourVector_t deuteron = {0.f,0.f,0.f,0.f};
//   unsigned char deu_prop = deu.properties;
//   if (deu.mother_pdg == 900010020 && (deu_prop & c)) {
//     deuteron = deu.vec;
//     for (const auto &pip : *pion_vector) {
//       FourVector_t piplus = {0.f,0.f,0.f,0.f};
//       unsigned char pip_prop = pip.properties;
//       if (pip.mother_id == deu.mother_id && (pip_prop & c)) {
//         piplus = pip.vec;
//         for (const auto &pim : *pion_vector) {
//           FourVector_t piminus = {0.f,0.f,0.f,0.f};
//           unsigned char pim_prop = pim.properties;
//           if (pim.mother_id == deu.mother_id && !(pim_prop & c)) {
//             piminus = pim.vec;
//             piplus += piminus;
//             piminus += deuteron;
//             DalitzPlot->Fill(piplus.M(), piminus.M());
//           }
//         }
//       }
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
//           mInv->Fill(v3.M(), v3.Pt());
//         }
//       }
//     }
//   }
// }



//______________________________________________________________________________

// // creo gli istogrammi con pdg e frequenza delle particelle madri
// deuMotherPdgFreq = new TH1F("deuteron's mother freq (pdg)", "titolo", deu_mother_pdg->size(), 0, deu_mother_pdg->size());
// piMotherPdgFreq = new TH1F("pion's mother freq (pdg)", "titolo", pi_mother_pdg->size(), 0, pi_mother_pdg->size());
// deuMotherPdgFreq->SetDirectory(0);
// piMotherPdgFreq->SetDirectory(0);
//
// sort(deu_mother_pdg->begin(), deu_mother_pdg->end());
// sort(pi_mother_pdg->begin(), pi_mother_pdg->end());
//
// for(unsigned i=0; i<deu_mother_pdg->size(); i++){
//   int d_pdg = deu_mother_pdg->at(i);
//   deuMotherPdgFreq->SetBinContent(i+1, deu_pdg_freq[d_pdg]);
//   deuMotherPdgFreq->GetXaxis()->SetBinLabel(i+1, Form("%i", d_pdg));
// }
//
// for(unsigned j=0; j<deu_mother_pdg->size(); j++){
//   int p_pdg = pi_mother_pdg->at(j);
//   piMotherPdgFreq->SetBinContent(j+1, pi_pdg_freq[p_pdg]);
//   piMotherPdgFreq->GetXaxis()->SetBinLabel(j+1, Form("%i", p_pdg));
// }
//
// deuMotherPdgFreq->Write();
// piMotherPdgFreq->Write();

//______________________________________________________________________________





// switch (selection) {
//   case 'A':
//   f_output->Write();
//   f_output->Close();
//   f_input->Close();
//   break;
//
//   case 'B':
//   f_outputMC->Write();
//   f_outputMC->Close();
//   f_inputMC->Close();
//   break;
//
//   default:
//   f_output->Write();
//   f_outputMC->Write();
//   f_input->Close();
//   f_inputMC->Close();
//   f_output->Close();
//   f_outputMC->Close();
// }

// }
