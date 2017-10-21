#include <iostream>
#include <climits>

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TMath.h>


using namespace std;


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     macro for calculate reconstruction efficiency for dStar               //
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

  // creo gli istogrammi per l' efficienza vs pT per materia e un per antimateria
  // per le 3 config (ITS TPC TOF)
  string tpctofMC[3] = {"TPC","TPC_TOF","TPC_(TOF)"};

  TH1F *hEfficiencyA[3];
  TH1F *hEfficiencyM[3];

  for(int iT=0; iT<3; iT++) {
    hEfficiencyA[iT] = new TH1F(Form("Efficiency_a (%s)", tpctofMC[iT].data()), "", 10, 0, 10);
    hEfficiencyM[iT] = new TH1F(Form("Efficiency_m (%s)", tpctofMC[iT].data()), "", 10, 0, 10);
  }

  // gli istogrammi sono contenuti in una lista, vado a leggerli da li
  TList *hList = static_cast<TList*>(f_input->Get("pfecchio_dStar"));

  TH2F *hTotal_a         = static_cast<TH2F*>(hList->FindObject("fTotal_dStar(2380)_a"));
  TH2F *hTotal_m         = static_cast<TH2F*>(hList->FindObject("fTotal_dStar(2380)_m"));
  hTotal_a->SetDirectory(0);
  hTotal_m->SetDirectory(0);

  TH1F *hTotalProjectionIntMass_a;
  TH1F *hTotalProjectionIntMass_m;
  TH1F *hTotalProjectionIntSpec_a;
  TH1F *hTotalProjectionIntSpec_m;

  TH2F *hReconstructed_a[3];
  TH2F *hReconstructed_m[3];

  for(int iT=0; iT<3; iT++){
    hReconstructed_a[iT] = static_cast<TH2F*>(hList->FindObject(Form("fRec_dStar(2380)_a_ITS_%s", tpctofMC[iT].data())));
    hReconstructed_m[iT] = static_cast<TH2F*>(hList->FindObject(Form("fRec_dStar(2380)_m_ITS_%s", tpctofMC[iT].data())));
    hReconstructed_a[iT]->SetDirectory(0);
    hReconstructed_m[iT]->SetDirectory(0);
  }

  // creo 2 vettori di TH1F che contengono le proiezioni sul pT degli istogrammi total e reconstructed
  TH1F *hTotalProjection_a[10];
  TH1F *hTotalProjection_m[10];
  TH1F *hReconstructedProjection_a[10][3];
  TH1F *hReconstructedProjection_m[10][3];

  double eff_a;
  double eff_m;
  double total_a;
  double total_m;
  double reconstructed_a;
  double reconstructed_m;

  for(int i=0; i<10; i++){

    int bin_min = i*2;
    int bin_max = (i*2)+2;

    hTotalProjection_a[i] = (TH1F*)hTotal_a->ProjectionX(Form("tot_a #it{p}_{T} [%i, %i] GeV/#it{c}", i, i+1), bin_min, bin_max);
    hTotalProjection_m[i] = (TH1F*)hTotal_m->ProjectionX(Form("tot_m #it{p}_{T} [%i, %i] GeV/#it{c}", i, i+1), bin_min, bin_max);
    hTotalProjection_a[i]->SetDirectory(0);
    hTotalProjection_m[i]->SetDirectory(0);
    total_a = hTotalProjection_a[i]->GetEntries();
    total_m = hTotalProjection_m[i]->GetEntries();

    for(int iT=0; iT<3; iT++){

      hReconstructedProjection_a[i][iT] = (TH1F*)hReconstructed_a[iT]->ProjectionX(Form("rec_a #it{p}_{T} [%i, %i] (%s) GeV/#it{c}", i, i+1, tpctofMC[iT].data()), bin_min, bin_max);
      hReconstructedProjection_m[i][iT] = (TH1F*)hReconstructed_m[iT]->ProjectionX(Form("rec_m #it{p}_{T} [%i, %i] (%s) GeV/#it{c}", i, i+1, tpctofMC[iT].data()), bin_min, bin_max);
      hReconstructedProjection_a[i][iT]->SetDirectory(0);
      hReconstructedProjection_m[i][iT]->SetDirectory(0);

      reconstructed_a = hReconstructedProjection_a[i][iT]->GetEntries();
      reconstructed_m = hReconstructedProjection_m[i][iT]->GetEntries();

      eff_a = reconstructed_a / total_a;
      eff_m = reconstructed_m / total_m;

      hEfficiencyA[iT]->SetBinContent(i+1, eff_a);
      hEfficiencyM[iT]->SetBinContent(i+1, eff_m);

      double error_eff_a = TMath::Sqrt(eff_a*(1-eff_a)/total_a);
      double error_eff_m = TMath::Sqrt(eff_m*(1-eff_m)/total_m);

      hEfficiencyA[iT]->SetBinError(i+1, error_eff_a);
      hEfficiencyM[iT]->SetBinError(i+1, error_eff_m);

    }
  }

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(1111);


  for(int iT=0; iT<3; iT++){
    hEfficiencyA[iT]->SetOption("PE");
    hEfficiencyA[iT]->SetLineColor(kBlue);
    hEfficiencyM[iT]->SetLineWidth(1);
    hEfficiencyA[iT]->SetStats(kFALSE);
    hEfficiencyA[iT]->SetMinimum(0.);
    hEfficiencyA[iT]->SetMaximum(1.1);
    hEfficiencyA[iT]->SetTitle(Form("dStar_A Reconstruction Efficiency %s", tpctofMC[iT].data()));
    hEfficiencyA[iT]->GetXaxis()->SetTitle("#it{p}_{T} GeV/#it{c}");
    hEfficiencyA[iT]->GetYaxis()->SetTitle("Acc #times #epsilon (%)");

    hEfficiencyM[iT]->SetOption("PE");
    hEfficiencyM[iT]->SetLineColor(kBlue);
    hEfficiencyM[iT]->SetLineWidth(1);
    hEfficiencyM[iT]->SetStats(kFALSE);
    hEfficiencyM[iT]->SetMinimum(0.);
    hEfficiencyM[iT]->SetMaximum(1.1);
    hEfficiencyM[iT]->SetTitle(Form("dStar_M Reconstruction Efficiency %s", tpctofMC[iT].data()));
    hEfficiencyM[iT]->GetXaxis()->SetTitle("#it{p}_{T} GeV/#it{c}");
    hEfficiencyM[iT]->GetYaxis()->SetTitle("Acc #times #epsilon (%)");
    if(iT != 0) hEfficiencyM[iT]->SetMaximum(0.1);
  }

  // // total production projection in mass for evaluating signal shape
  // hTotalProjectionIntMass_a = (TH1F*)hTotal_a->ProjectionX("fMassTotSpec_A");
  // hTotalProjectionIntMass_m = (TH1F*)hTotal_m->ProjectionX("fMassTotSpec_M");
  // hTotalProjectionIntMass_a->SetDirectory(0);
  // hTotalProjectionIntMass_m->SetDirectory(0);
  //
  // // total production projection pT spectrum for evaluating signal shape
  // hTotalProjectionIntSpec_a = (TH1F*)hTotal_a->ProjectionY("fSpecTotSpec_A");
  // hTotalProjectionIntSpec_m = (TH1F*)hTotal_m->ProjectionY("fSpecTotSpec_M");
  // hTotalProjectionIntSpec_a->SetDirectory(0);
  // hTotalProjectionIntSpec_m->SetDirectory(0);
  // hTotalProjectionIntSpec_a->Rebin(4);
  // hTotalProjectionIntSpec_m->Rebin(4);
  //
  // // histo make up
  // hTotalProjectionIntMass_a->SetMarkerStyle(20);
  // hTotalProjectionIntMass_a->SetMarkerColor(kBlue);
  // hTotalProjectionIntMass_a->SetLineColor(kBlue);
  // hTotalProjectionIntMass_a->SetMarkerSize(1);
  // hTotalProjectionIntMass_a->SetStats(false);
  // hTotalProjectionIntMass_a->SetMinimum(0.);
  // hTotalProjectionIntMass_m->SetMarkerStyle(20);
  // hTotalProjectionIntMass_m->SetMarkerColor(kBlue);
  // hTotalProjectionIntMass_m->SetLineColor(kBlue);
  // hTotalProjectionIntMass_m->SetMarkerSize(1);
  // hTotalProjectionIntMass_m->SetStats(false);
  // hTotalProjectionIntMass_m->SetMinimum(0.);
  //
  //
  // hTotalProjectionIntMass_a->Write();
  // hTotalProjectionIntMass_m->Write();

  // creo il canvas con i vari plot dell efficienza fatti bene
  TCanvas* canEff[3];
  // TLegend* mLeg = new TLegend(0.17,0.52,0.45,0.73);

  for(int iT=0; iT<3; iT++){

    int altezza    = 900;
    int larghezza  = 700;
    canEff[iT] = new TCanvas(Form("cEfficiency_ITC_%s",tpctofMC[iT].data()),Form("cEfficiency_ITS_%s",tpctofMC[iT].data()), altezza, larghezza);
    canEff[iT]->Divide(2,1);

    hEfficiencyA[iT]->SetTitle("d^{*} (2380) A");
    canEff[iT]->cd(1);
    hEfficiencyA[iT]->SetMarkerStyle(20);
    hEfficiencyA[iT]->SetMarkerSize(1);
    hEfficiencyA[iT]->Draw("PE");

    hEfficiencyM[iT]->SetTitle("d^{*} (2380) M");
    canEff[iT]->cd(2);
    hEfficiencyM[iT]->SetMarkerStyle(20);
    hEfficiencyM[iT]->SetMarkerSize(1);
    hEfficiencyM[iT]->Draw("PE");


    TPaveText pave(0.2,0.6,0.4,0.65,"blNDC");
    pave.AddText(tpctofMC[iT].data());
    pave.SetBorderSize(1);
    pave.SetFillColor(0);
    pave.Draw();
    // canEff[iT]->Print(Form("Efficiency_ITS_%s.pdf",tpctofMC[iT].data()));
    canEff[iT]->Write();
  }

  f_output->Write();
  f_output->Close();
  f_input->Close();

}
