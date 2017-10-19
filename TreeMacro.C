#include <iostream>
#include <climits>
#include <map>

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

// #include <PdgToName.h>

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t ;


void TreeAnalysis(){

  // file input-output
  TFile *f_input         = new TFile("../20171016_2327_grid/dStarTree.root", "READ");
  TFile *f_output        = new TFile("macroOut/TreeAnalysis.root", "RECREATE");

  // varie cose per leggere il Tree
  TH1F  *deuPt = new TH1F("Deuteron #it{p}_{T}", "titolo", 40, 0, 20);
  TH1F  *piPt  = new TH1F("Pion #it{p}_{T}", "titolo", 20, 0, 10);

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

  // creo un istogramma con pdg e frequenza
  TH1F *deuMotherPdgFreq;
  TH1F *piMotherPdgFreq;

  // creo l' istogramma della massa invariate della dStar
  TH2F *mInv = new TH2F("m_invariante", "titolo", 2950, 2.100, 8.000, 50, 0.0, 10);
  mInv->SetDirectory(0);

  //______________________________________________________________________________

  for(int i=0; i < tree->GetEntries(); i++){

    deuteron_vector->clear();
    pion_vector->clear();
    tree->GetEntry(i);


    // loop sul vettore contenente i deutoni
    for(auto &deu : *deuteron_vector) {
      FourVector_t  deuvec;
      deuPt->Fill(deuvec.Pt());

      int pdg = deu.mother_pdg;
      auto it1 = std::find(deu_mother_pdg->begin(), deu_mother_pdg->end(), pdg);
      if(it1 == deu_mother_pdg->end()){
        deu_pdg_freq[pdg] = 1;
        deu_mother_pdg->push_back(pdg);
      } else {
        deu_pdg_freq[pdg] += 1;
      }
    }

    // loop sul vettore contenente i pioni (sia + che -)
    for(auto &pi : *pion_vector) {
      FourVector_t pivec = pi.vec;
      piPt->Fill(pivec.Pt());

      int pdg = pi.mother_pdg;
      auto it2 = std::find(pi_mother_pdg->begin(), pi_mother_pdg->end(), pdg);
      if(it2 == pi_mother_pdg->end()){
        pi_pdg_freq[pdg] = 1;
        pi_mother_pdg->push_back(pdg);
      } else {
        pi_pdg_freq[pdg] += 1;
      }
    }

    //______________________________________________________________________________

    // creo l istogramma della massa invariante
    FourVector_t v1, v2, v3;

    for(auto &deu : *deuteron_vector) {
      v1 = deu.vec;
      for(auto &pip : *pion_vector) {
        if(pip.charge) {
          v2 = pip.vec;
          for(auto &pim : *pion_vector) {
            if(!pim.charge) {
              v3 = pim.vec;
              v3 += v2;
              v3 += v1;
              mInv->Fill(v3.M(), v3.Pt());
            }
          }
        }
      }
    }
  }

  //______________________________________________________________________________

  // creo gli istogrammi con pdg e frequenza delle particelle madri
  deuMotherPdgFreq = new TH1F("deuteron's mother freq (pdg)", "titolo", deu_mother_pdg->size(), 0, deu_mother_pdg->size());
  piMotherPdgFreq = new TH1F("pion's mother freq (pdg)", "titolo", pi_mother_pdg->size(), 0, pi_mother_pdg->size());
  deuMotherPdgFreq->SetDirectory(0);
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

  mInv->Write();

  f_output->Write();
  f_output->Close();
  f_input->Close();


  // cout << endl;
  // cout << "deuteron mothers pdg list: " << endl;
  // for(auto &d : *deu_mother_pdg){
  //   cout << d << "   ";
  // }
  // cout << endl;
  // cout << "______________________________________";
  // cout << endl;
  // cout << endl;
  //
  // cout << "pion mothers pdg list: " << endl;
  // for(auto &p : *pi_mother_pdg){
  //   cout << p << "   ";
  // }
  // cout << endl;
  // cout << "______________________________________";
  // cout << endl;
  // cout << endl;

  // for(auto &ii : deu_pdg_freq){
  //   cout << ii.first << " ==> " << ii.second << endl;
  // }
  //
  // cout << endl;
  // cout << "______________________________________";
  // cout << endl;
  // cout << endl;
  //
  // for(auto &ij : pi_pdg_freq){
  //   cout << ij.first << " ==> " << ij.second << endl;
  // }
  //
  // cout << endl;
  // cout << "______________________________________";
  // cout << endl;
  // cout << endl;

}
