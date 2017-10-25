#include <iostream>
#include <climits>

#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TRatioPlot.h>
#include <TPaveStats.h>
#include <TGaxis.h>
#include <TMath.h>


using namespace std;


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     macro for calculate reconstruction efficiency for dStar               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void Efficiency(){

  TFile *f_input         = new TFile("../analysisresults/grid_analysis/1020_1700_grid/AnalysisResults.root", "READ");
  TFile *f_output        = new TFile("../output/out_grid/Efficiency.root", "RECREATE");


  if(f_input->IsOpen()){
    printf("File opened successfully\n");
  }else{
    printf("File NOT opened successfully\n");
    exit(0);
  }

  gStyle->SetOptStat(1101);
  gStyle->SetTextFont(42);
  gStyle->SetPadTickY(1);


  // creo gli istogrammi per l' efficienza vs pT per materia e un per antimateria
  // per le 3 config (ITS TPC TOF)
  string tpctofMC[3] = {"TPC","TPC_TOF","TPC_(TOF)"};

  TH1F *hEfficiencyA[3];
  TH1F *hEfficiencyM[3];


  for(int iT=0; iT<3; iT++) {
    hEfficiencyA[iT] = new TH1F(Form("Efficiency_a (%s)", tpctofMC[iT].data()), "", 10, 0, 10);
    hEfficiencyM[iT] = new TH1F(Form("Efficiency_m (%s)", tpctofMC[iT].data()), "", 10, 0, 10);
    hEfficiencyA[iT]->SetDirectory(0);
    hEfficiencyM[iT]->SetDirectory(0);
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

  // efficiency calculation
  for(int i=0; i<10; i++){

    int bin_min = i*2;
    int bin_max = (i*2)+2;

    hTotalProjection_a[i] = (TH1F*)hTotal_a->ProjectionX(Form("tot_a #it{p}_{T} [%i, %i] #it{GeV/c}", i, i+1), bin_min, bin_max);
    hTotalProjection_m[i] = (TH1F*)hTotal_m->ProjectionX(Form("tot_m #it{p}_{T} [%i, %i] #it{GeV/c}", i, i+1), bin_min, bin_max);
    hTotalProjection_a[i]->SetDirectory(0);
    hTotalProjection_m[i]->SetDirectory(0);
    total_a = hTotalProjection_a[i]->GetEntries();
    total_m = hTotalProjection_m[i]->GetEntries();

    // efficiency calculation for 3 PID configuration
    for(int iT=0; iT<3; iT++){

      hReconstructedProjection_a[i][iT] = (TH1F*)hReconstructed_a[iT]->ProjectionX(Form("rec_a #it{p}_{T} [%i, %i] (%s) #it{Gev/c}", i, i+1, tpctofMC[iT].data()), bin_min, bin_max);
      hReconstructedProjection_m[i][iT] = (TH1F*)hReconstructed_m[iT]->ProjectionX(Form("rec_m #it{p}_{T} [%i, %i] (%s) #it{Gev/c}", i, i+1, tpctofMC[iT].data()), bin_min, bin_max);
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

  double xoffset = 0.9;
  double yoffset = 1;
  double titlesize = 0.045;

  for(int iT=0; iT<3; iT++){

    hEfficiencyA[iT]->SetOption("PE");
    hEfficiencyA[iT]->SetStats(kFALSE);
    hEfficiencyA[iT]->SetMarkerStyle(21);
    hEfficiencyA[iT]->SetMarkerSize(0.8);
    hEfficiencyA[iT]->SetMarkerColor(kBlue);
    hEfficiencyA[iT]->SetLineColor(kBlue);
    hEfficiencyA[iT]->SetLineWidth(2);
    hEfficiencyA[iT]->SetMinimum(-0.02);
    hEfficiencyA[iT]->SetMaximum(1.1);
    hEfficiencyA[iT]->SetTitle("");
    // hEfficiencyA[iT]->SetTitle(Form("dStar Reconstruction Efficiency ITS_%s  a", tpctofMC[iT].data()));
    hEfficiencyA[iT]->GetXaxis()->SetTitle("#it{p}_{T} #it{GeV/c}");
    hEfficiencyA[iT]->GetYaxis()->SetTitle("Acc #times #epsilon");
    hEfficiencyA[iT]->GetXaxis()->SetTitleSize(titlesize);
    hEfficiencyA[iT]->GetYaxis()->SetTitleSize(titlesize);
    hEfficiencyA[iT]->GetXaxis()->SetTitleOffset(xoffset);
    hEfficiencyA[iT]->GetYaxis()->SetTitleOffset(yoffset);


    hEfficiencyM[iT]->SetOption("PE");
    hEfficiencyM[iT]->SetStats(kFALSE);
    hEfficiencyM[iT]->SetMarkerStyle(21);
    hEfficiencyM[iT]->SetMarkerSize(0.8);
    hEfficiencyM[iT]->SetMarkerColor(kBlue);
    hEfficiencyM[iT]->SetLineColor(kBlue);
    hEfficiencyM[iT]->SetLineWidth(2);
    hEfficiencyM[iT]->SetMinimum(-0.02);
    hEfficiencyM[iT]->SetMaximum(1.1);
    // hEfficiencyM[iT]->SetTitle(Form("dStar Reconstruction Efficiency ITS_%s  m", tpctofMC[iT].data()));
    hEfficiencyM[iT]->SetTitle("");
    hEfficiencyM[iT]->GetXaxis()->SetTitle("#it{p}_{T} #it{GeV/c}");
    hEfficiencyM[iT]->GetYaxis()->SetTitle("Acc #times #epsilon");
    hEfficiencyM[iT]->GetXaxis()->SetTitleSize(titlesize);
    hEfficiencyM[iT]->GetYaxis()->SetTitleSize(titlesize);
    hEfficiencyM[iT]->GetXaxis()->SetTitleOffset(xoffset);
    hEfficiencyM[iT]->GetYaxis()->SetTitleOffset(yoffset);

    if(iT != 0) {
      hEfficiencyM[iT]->SetMaximum(0.31);
      hEfficiencyA[iT]->SetMaximum(0.31);
    }
    hEfficiencyM[iT]->Write();
    hEfficiencyA[iT]->Write();
  }


  // total production projection in mass for evaluating signal shape
  hTotalProjectionIntMass_a = (TH1F*)hTotal_a->ProjectionX("fMassTotSpec_A");
  hTotalProjectionIntMass_m = (TH1F*)hTotal_m->ProjectionX("fMassTotSpec_M");
  hTotalProjectionIntMass_a->SetDirectory(0);
  hTotalProjectionIntMass_m->SetDirectory(0);

  // total production projection pT spectrum for evaluating signal shape
  hTotalProjectionIntSpec_a = (TH1F*)hTotal_a->ProjectionY("fSpecTotSpec_A");
  hTotalProjectionIntSpec_m = (TH1F*)hTotal_m->ProjectionY("fSpecTotSpec_M");
  hTotalProjectionIntSpec_a->SetDirectory(0);
  hTotalProjectionIntSpec_m->SetDirectory(0);
  hTotalProjectionIntSpec_a->Rebin(4);
  hTotalProjectionIntSpec_m->Rebin(4);

  // histo make up
  hTotalProjectionIntMass_a->SetLineColor(kBlue);
  hTotalProjectionIntMass_a->SetLineWidth(2);
  hTotalProjectionIntMass_a->SetFillColor(kBlue);
  hTotalProjectionIntMass_a->SetFillStyle(3005);
  hTotalProjectionIntMass_a->SetMinimum(0.);
  hTotalProjectionIntMass_a->GetXaxis()->SetTitle("M #it{GeV/c^{2}}");
  hTotalProjectionIntMass_a->GetYaxis()->SetTitle("Counts");
  hTotalProjectionIntMass_a->GetYaxis()->SetTitleOffset(1.3);
  hTotalProjectionIntMass_a->SetName("d^{*} (A) invariant mass");

  //imparare a lavorare con la statbox
  hTotalProjectionIntMass_m->SetLineColor(kBlue);
  hTotalProjectionIntMass_m->SetLineWidth(2);
  hTotalProjectionIntMass_m->SetFillColor(kBlue);
  hTotalProjectionIntMass_m->SetFillStyle(3005);
  hTotalProjectionIntMass_m->SetMinimum(0.);
  hTotalProjectionIntMass_m->GetXaxis()->SetTitle("M #it{GeV/c^{2}}");
  hTotalProjectionIntMass_m->GetYaxis()->SetTitle("Counts");

  TCanvas *canMass = new TCanvas("cMass_a", "", 400, 400);
  hTotalProjectionIntMass_a->Draw();

  canMass->SaveAs("mass_signal.pdf");
  canMass->Write();
  canMass->Close();
  hTotalProjectionIntMass_a->Write();
  hTotalProjectionIntMass_m->Write();


  // creo i canvas con i vari plot dell efficienza fatti bene
  int altezza    = 400;
  int larghezza  = 800;

  TCanvas* canEff[3];

  for(int iT=0; iT<3; iT++){

    canEff[iT] = new TCanvas(Form("cEfficiency_ITC_%s",tpctofMC[iT].data()), Form("cEfficiency_ITS_%s",tpctofMC[iT].data()), larghezza, altezza);
    canEff[iT]->Divide(2,1);

    canEff[iT]->cd(1);
    hEfficiencyA[iT]->Draw("PE");

    canEff[iT]->cd(2);
    hEfficiencyM[iT]->Draw("PE");

    canEff[iT]->SaveAs(Form("eff_ITS_%s.pdf",tpctofMC[iT].data()));
    canEff[iT]->Write();
    canEff[iT]->Close();
  }

  // compared efficiency ITS-TPC Vs ITS-TPC-TOF
  TCanvas *canEffCon = new TCanvas("cEfficiency_con", "cEfficiency_con", larghezza, altezza);
  TLegend* legEffA = new TLegend(0.182, 0.678, 0.462, 0.890);
  legEffA->SetTextSize(0.040);
  TLegend* legEffM = new TLegend(0.182, 0.678, 0.462, 0.890);
  legEffM->SetTextSize(0.040);
  canEffCon->Divide(2,1);


  canEffCon->cd(1);
  hEfficiencyA[0]->Draw("PE");
  hEfficiencyA[1]->SetMarkerStyle(22);
  hEfficiencyA[1]->SetMarkerColor(kRed);
  hEfficiencyA[1]->SetLineColor(kRed);
  hEfficiencyA[1]->SetMarkerSize(1);
  hEfficiencyA[1]->Draw("PESAME");
  legEffA->SetBorderSize(0);
  legEffA->SetHeader("Reconstruction efficiency A");
  legEffA->AddEntry(hEfficiencyA[0], "ITS-TPC");
  legEffA->AddEntry(hEfficiencyA[1], "ITS-TPC-TOF");
  legEffA->Draw();


  canEffCon->cd(2);
  hEfficiencyM[0]->Draw("PE");
  hEfficiencyM[1]->SetMarkerStyle(22);
  hEfficiencyM[1]->SetMarkerColor(kRed);
  hEfficiencyM[1]->SetLineColor(kRed);
  hEfficiencyM[1]->SetMarkerSize(1);
  hEfficiencyM[1]->Draw("PESAME");
  legEffM->SetBorderSize(0);
  legEffM->SetHeader("Reconstruction efficiency M");
  legEffM->AddEntry(hEfficiencyM[0], "ITS-TPC");
  legEffM->AddEntry(hEfficiencyM[1], "ITS-TPC-TOF");
  legEffM->Draw();


  canEffCon->SaveAs("eff_compared.pdf");
  canEffCon->Write();
  canEffCon->Close();

  // compared efficiency ITS-TPC Vs ITS-TPC-TOF Vs ITS-TPC-(TOF)
  TCanvas *canEffConDeu3 = new TCanvas("cEfficiency_con3", "cEfficiency_con3", larghezza, altezza);
  TLegend* legEffA3 = new TLegend(0.180, 0.623, 0.488, 0.833);
  legEffA3->SetTextSize(0.040);
  TLegend* legEffM3 = new TLegend(0.180, 0.623, 0.488, 0.833);
  legEffM3->SetTextSize(0.040);
  canEffConDeu3->Divide(2,1);

  canEffConDeu3->cd(1);
  hEfficiencyA[0]->Draw("PE");
  hEfficiencyA[1]->SetMarkerStyle(22);
  hEfficiencyA[1]->SetMarkerColor(kRed);
  hEfficiencyA[1]->SetLineColor(kRed);
  hEfficiencyA[1]->SetMarkerSize(1);
  hEfficiencyA[1]->Draw("PESAME");
  hEfficiencyA[2]->SetMarkerStyle(20);
  hEfficiencyA[2]->SetMarkerColor(kMagenta);
  hEfficiencyA[2]->SetLineColor(kMagenta);
  hEfficiencyA[2]->SetMarkerSize(1);
  hEfficiencyA[2]->Draw("PESAME");
  legEffA3->SetBorderSize(0);
  legEffA3->SetHeader("Reconstruction efficiency A");
  legEffA3->AddEntry(hEfficiencyM[0], "ITS-TPC");
  legEffA3->AddEntry(hEfficiencyM[1], "ITS-TPC-TOF");
  legEffA3->AddEntry(hEfficiencyM[2], "ITS-TPC-TOF_{deuteron}");
  legEffA3->Draw();

  canEffConDeu3->cd(2);
  hEfficiencyM[0]->Draw("PE");
  hEfficiencyM[1]->SetMarkerStyle(22);
  hEfficiencyM[1]->SetMarkerColor(kRed);
  hEfficiencyM[1]->SetLineColor(kRed);
  hEfficiencyM[1]->SetMarkerSize(1);
  hEfficiencyM[1]->Draw("PESAME");
  hEfficiencyM[2]->SetMarkerStyle(20);
  hEfficiencyM[2]->SetMarkerColor(kMagenta);
  hEfficiencyM[2]->SetLineColor(kMagenta);
  hEfficiencyM[2]->SetMarkerSize(1);
  hEfficiencyM[2]->Draw("PESAME");
  legEffM3->SetBorderSize(0);
  legEffM3->SetHeader("Reconstruction efficiency M");
  legEffM3->AddEntry(hEfficiencyM[0], "ITS-TPC");
  legEffM3->AddEntry(hEfficiencyM[1], "ITS-TPC-TOF");
  legEffM3->AddEntry(hEfficiencyM[2], "ITS-TPC-TOF_{deuteron}");
  legEffM3->Draw();

  canEffConDeu3->SaveAs("eff_compared3.pdf");
  canEffConDeu3->Write();
  canEffConDeu3->Close();

  // efficiency ratio TOF A/M
  TCanvas *canEffRatio = new TCanvas("cEfficiency_ratio", "cEfficiency_con", 400, 400);

  hEfficiencyM[1]->SetMarkerStyle(21);
  hEfficiencyM[1]->SetMarkerColor(kBlue);
  hEfficiencyM[1]->SetLineColor(kBlue);
  hEfficiencyM[1]->SetMarkerSize(1);
  hEfficiencyM[1]->SetName("Efficiency M #scale[0.8]{(ITS-TPC-TOF)}");
  hEfficiencyA[1]->SetMarkerStyle(22);
  hEfficiencyA[1]->SetMarkerColor(kRed);
  hEfficiencyA[1]->SetLineColor(kRed);
  hEfficiencyA[1]->SetMarkerSize(1);
  hEfficiencyA[1]->SetName("Efficiency A #scale[0.8]{(ITS-TPC-TOF)}");

  // canEffRatio->SetTicks(0, 1);
  auto rp = new TRatioPlot(hEfficiencyM[1], hEfficiencyA[1], "diff");

  rp->SetH1DrawOpt("PE");
  rp->SetH2DrawOpt("PE");
  rp->Draw();

  rp->GetLowerRefYaxis()->SetTitle("difference");
  rp->GetLowerRefYaxis()->SetTitleOffset(1.3);
  rp->GetUpperRefYaxis()->SetTitleOffset(1.2);
  rp->GetUpperRefYaxis()->SetRangeUser(-0.01, 0.31);
  rp->GetLowerRefYaxis()->SetLabelSize(0.025);
  rp->GetLowerRefYaxis()->SetRangeUser(-0.012, 0.074);
  rp->GetLowerRefGraph()->SetLineColor(kBlack);
  rp->GetLowerRefGraph()->SetMarkerColor(kBlack);
  TLegend *legEffDiff = rp->GetUpperPad()->BuildLegend(0.439, 0.682, 0.747, 0.892);
  legEffDiff->SetBorderSize(0);
  legEffDiff->SetTextSize(0.050);
  legEffDiff->SetHeader("");
  canEffRatio->Update();

  canEffRatio->SaveAs("eff_diff.pdf");
  canEffRatio->Write();
  // canEffRatio->Close();


  f_output->Write();
  f_output->Close();
  f_input->Close();

}
