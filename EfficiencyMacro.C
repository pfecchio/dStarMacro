#include <iostream>
#include <climits>

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TCanvas.h>
#include <TMath.h>


using namespace std;


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     macro for calculate reconstruction efficiency for dStar decay         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void Efficiency(){

  TFile *f_input         = new TFile("../20171010_2315_grid/AnalysisResults.root", "READ");
  TFile *f_output        = new TFile("macroOut/Efficiency.root", "RECREATE");

  if(f_input->IsOpen()){
    printf("File opened successfully\n");
  }else{
    printf("File NOT opened successfully\n");
    exit(0);
  }

  // creo un istogramma per l' efficienza vs pT per materia e un per antimateria
  TH1D *hEfficiencyA = new TH1D("Efficiency_a", "titolo", 10, 0, 10);
  TH1D *hEfficiencyM = new TH1D("Efficiency_m", "titolo", 10, 0, 10);

  // gli istogrammi sono contenuti in una lista, vado a leggerli da li
  TList *hList = static_cast<TList*>(f_input->Get("pfecchio_dStar"));

  TH2F *hTotal_a         = static_cast<TH2F*>(hList->FindObject("fTotal_dStar(2380)_a"));
  TH2F *hTotal_m         = static_cast<TH2F*>(hList->FindObject("fTotal_dStar(2380)_m"));
  hTotal_a->SetDirectory(0);
  hTotal_m->SetDirectory(0);

  TH2F *hReconstructed_a = static_cast<TH2F*>(hList->FindObject("fRec_dStar(2380)_a_ITS_TPC"));
  TH2F *hReconstructed_m = static_cast<TH2F*>(hList->FindObject("fRec_dStar(2380)_m_ITS_TPC"));
  hReconstructed_a->SetDirectory(0);
  hReconstructed_m->SetDirectory(0);

  // creo 2 vettori di TH1F che contengono le proiezioni sul pT degli istogrammi total e recopnstructed
  TH1D *hTotalProjection_a[10];
  TH1D *hTotalProjection_m[10];
  TH1D *hReconstructedProjection_a[10];
  TH1D *hReconstructedProjection_m[10];

  double eff_a;
  double eff_m;
  double total_a;
  double total_m;
  double reconstructed_a;
  double reconstructed_m;

  for(int i=0; i<10; i++){

    int bin_min = i*2;
    int bin_max = (i*2)+2;

    hTotalProjection_a[i] = hTotal_a->ProjectionX(Form("tot_a #it{p}_{T} [%i, %i] GeV/#it{c}", i, i+1), bin_min, bin_max);
    hTotalProjection_m[i] = hTotal_m->ProjectionX(Form("tot_m #it{p}_{T} [%i, %i] GeV/#it{c}", i, i+1), bin_min, bin_max);
    hTotalProjection_a[i]->SetDirectory(0);
    hTotalProjection_m[i]->SetDirectory(0);

    hReconstructedProjection_a[i] = hReconstructed_a->ProjectionX(Form("rec_a #it{p}_{T} [%i, %i] GeV/#it{c}", i, i+1), bin_min, bin_max);
    hReconstructedProjection_m[i] = hReconstructed_m->ProjectionX(Form("rec_m #it{p}_{T} [%i, %i] GeV/#it{c}", i, i+1), bin_min, bin_max);
    hReconstructedProjection_a[i]->SetDirectory(0);
    hReconstructedProjection_m[i]->SetDirectory(0);

    total_a = hTotalProjection_a[i]->GetEntries();
    total_m = hTotalProjection_m[i]->GetEntries();
    reconstructed_a = hReconstructedProjection_a[i]->GetEntries();
    reconstructed_m = hReconstructedProjection_m[i]->GetEntries();

    eff_a = reconstructed_a / total_a;
    eff_m = reconstructed_m / total_m;

    hEfficiencyA->SetBinContent(i+1, eff_a);
    hEfficiencyM->SetBinContent(i+1, eff_m);

    double error_eff_a = TMath::Sqrt(eff_a*(1-eff_a)/total_a);
    double error_eff_m = TMath::Sqrt(eff_m*(1-eff_m)/total_m);

    hEfficiencyA->SetBinError(i+1, error_eff_a);
    hEfficiencyM->SetBinError(i+1, error_eff_m);

  }

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);


  hEfficiencyA->SetOption("HISTE1");
  hEfficiencyA->SetLineColor(kBlue);
  // hEfficiencyA->SetFillColorAlpha(kBlue, 0.10);
  hEfficiencyA->SetLineWidth(2.75);
  hEfficiencyA->SetTitle("dStar_A Reconstruction Efficiency");
  hEfficiencyA->GetXaxis()->SetTitle("#it{p}_{T} GeV/#it{c}");
  hEfficiencyA->GetYaxis()->SetTitle("Efficiency #epsilon");

  hEfficiencyM->SetOption("HISTE1");
  hEfficiencyM->SetLineColor(kBlue);
  // hEfficiencyM->SetFillColorAlpha(kBlue, 0.10);
  hEfficiencyM->SetLineWidth(2.75);
  hEfficiencyM->SetTitle("dStar_M Reconstruction Efficiency");
  hEfficiencyM->GetXaxis()->SetTitle("#it{p}_{T} GeV/#it{c}");
  hEfficiencyM->GetYaxis()->SetTitle("Efficiency #epsilon");

  // TCanvas *c = new TCanvas("M_A efficiency","matter-antimatter efficiency", 1200, 500);
  // c->Divide(2,1);
  // c->cd(1);
  // hEfficiencyA->Draw();
  // c->cd(2);
  // hEfficiencyM->Draw();
  // c->Write();
  // c->Close();

  f_output->Write();
  f_output->Close();
  f_input->Close();

}
