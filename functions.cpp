#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>

#include "functions.hpp"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

void SetStyle() {
  gStyle->SetHistLineColor(kBlue + 2);
  gStyle->SetHistFillColor(kBlue);
  gStyle->SetHistFillStyle(1001); // solid
  gStyle->SetFrameFillColor(kYellow);
  gStyle->SetFrameFillStyle(3005); // bent lines in background
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetTitleYOffset(0.9);
  gStyle->SetOptStat(2210); // show histogram #entries, mean and rms with errors
  gStyle->SetOptFit(1111);  // show fit probability, fit chi/ndf, fit parameters
}

int DataAnalysis(char const *filename) {

  TFile *KDecay = new TFile(filename, "UPDATE");
  // histos must be fitted and overwritten

  TH1F *hIndex = (TH1F *)KDecay->Get("ParticleIndex");
  TH1F *hTheta = (TH1F *)KDecay->Get("ThetaAngle");
  TH1F *hPhi = (TH1F *)KDecay->Get("PhiAngle");
  TH1F *hImpulse = (TH1F *)KDecay->Get("ParticleImpulse");
  TH1F *hTransvImpulse = (TH1F *)KDecay->Get("ParticleTransvImpulse");
  TH1F *hEnergy = (TH1F *)KDecay->Get("ParticleEnergy");
  TH1F *hTotalInvMass = (TH1F *)KDecay->Get("TotalInvMass");
  TH1F *hUnlikeInvMass = (TH1F *)KDecay->Get("UnlikeInvMass");
  TH1F *hLikeInvMass = (TH1F *)KDecay->Get("LikeInvMass");
  TH1F *hUnlikeKPInvMass = (TH1F *)KDecay->Get("UnlikeKPInvMass");
  TH1F *hLikeKPInvMass = (TH1F *)KDecay->Get("LikeKPInvMass");
  TH1F *hDaughtersInvMass = (TH1F *)KDecay->Get("DaughtersInvMass");

  hUnlikeInvMass->Sumw2();
  hLikeInvMass->Sumw2();
  hUnlikeKPInvMass->Sumw2();
  hLikeKPInvMass->Sumw2();

  TH1F *hResult1 = new TH1F("Result1", "Comparison histogram #1", 500, 0., 4.);
  TH1F *hResult2 = new TH1F("Result2", "Comparison histogram #2", 500, 0., 4.);
  // hResult2 is more efficient: only K and P particle are considered,
  // so it should be more precise

  hResult1->Add(hUnlikeInvMass, hLikeInvMass, 1, -1);
  hResult2->Add(hUnlikeKPInvMass, hLikeKPInvMass, 1, -1);

  hResult1->Write();
  hResult2->Write();

  TF1 *fitfunc0 = new TF1("UnifDistrib", "[0]", 0., 4.);
  TF1 *fitfunc1 = new TF1("ExponDistrib", "exp([0]-x/[1])", 0., 4.);
  // parameter [1] is exponential mean
  fitfunc1->SetParLimits(0, 0., 100.);
  fitfunc1->SetParLimits(1, 0., 10.);
  TF1 *fitfunc2 = new TF1("GaussDistrib", "gaus", 0., 2.);
  fitfunc2->SetParLimits(0, 100., 10000.);
  fitfunc2->SetParLimits(1, 0., 1.5);
  fitfunc2->SetParLimits(2, 0., 0.5);

  std::cout << "-------------- Histogram #Entries ---------------\n";
  Int_t nEventGen = 1E5;
  Int_t nPartGen = 100;
  std::cout << "Expected value for the entries: " << nEventGen * nPartGen
            << '\n'
            << "Generated particles histo #entries: " << hIndex->GetEntries()
            << '\n'
            << "Theta histo #entries: " << hTheta->GetEntries() << '\n'
            << "Phi histo #entries: " << hPhi->GetEntries() << '\n'
            << "Impulse histo #entries: " << hImpulse->GetEntries() << '\n'
            << "Transverse Impulse histo #entries: "
            << hTransvImpulse->GetEntries() << '\n'
            << "Energy histo #entries: " << hEnergy->GetEntries() << '\n'
            << "---- The following should have a lot more entries ---\n"
            << "Total Invariant Mass histo #entries: "
            << hTotalInvMass->GetEntries() << '\n'
            << "Opposite charges Invariant Mass histo #entries: "
            << hUnlikeInvMass->GetEntries() << '\n'
            << "Same charges Invariant Mass histo #entries: "
            << hLikeInvMass->GetEntries() << '\n'
            << "Unlike charge Kaon/Pion Pair Invariant Mass histo #entries: "
            << hUnlikeKPInvMass->GetEntries() << '\n'
            << "Like charge Kaon/Pion Pair Invariant Mass histo #entries: "
            << hLikeKPInvMass->GetEntries() << '\n'
            << "Decay daughters Invariant Mass histo #entries: "
            << hDaughtersInvMass->GetEntries() << '\n'
            << "Comparison histo#1 #entries: " << hResult1->GetEntries() << "\n"
            << "Comparison histo#2 #entries: " << hResult2->GetEntries()
            << "\n";
  std::cout << "----------- Particle generation proportions -----------\n";
  for (Int_t i = 1; i != 8; ++i) {
    std::cout << "Bin content is " << hIndex->GetBinContent(i) << " +/- "
              << hIndex->GetBinError(i) << "\n"
              << "Proportion obtained is "
              << hIndex->GetBinContent(i) / hIndex->GetEntries() << "\n";
  }
  std::cout << "----------- Fitting Uniform Distributions -----------\n";
  hTheta->Fit("UnifDistrib", "Q");
  TF1 *UniformTheta = hTheta->GetFunction("UnifDistrib");
  hTheta->Write();
  std::cout << "Theta Fit Chi^2/NDF: " << UniformTheta->GetChisquare() << " / "
            << UniformTheta->GetNDF() << '\n'
            << "Theta Fit Probability: " << UniformTheta->GetProb() << '\n'
            << "Theta Fit Parameter: " << UniformTheta->GetParameter(0)
            << " +/- " << UniformTheta->GetParError(0) << '\n';

  hPhi->Fit("UnifDistrib", "Q");
  TF1 *UniformPhi = hPhi->GetFunction("UnifDistrib");
  hPhi->Write();
  std::cout << "Phi Fit Chi^2/NDF: " << UniformPhi->GetChisquare() << " / "
            << UniformPhi->GetNDF() << '\n'
            << "Phi Fit Probability: " << UniformPhi->GetProb() << '\n'
            << "Phi Fit Parameter: " << UniformPhi->GetParameter(0) << " +/- "
            << UniformPhi->GetParError(0) << '\n';

  std::cout << "----------- Fitting Exponential Distributions -----------\n";
  hImpulse->Fit("ExponDistrib", "Q");
  TF1 *ExpImpulse = hImpulse->GetFunction("ExponDistrib");
  hImpulse->Write();
  std::cout << "Impulse Fit Chi^2/NDF: " << ExpImpulse->GetChisquare() << " / "
            << ExpImpulse->GetNDF() << '\n'
            << "Impulse Fit Probability: " << ExpImpulse->GetProb() << '\n'
            << "Impulse Fit Par0: " << ExpImpulse->GetParameter(0) << " +/- "
            << ExpImpulse->GetParError(0) << '\n'
            << "Impulse Fit Par1: " << ExpImpulse->GetParameter(1) << " +/- "
            << ExpImpulse->GetParError(1) << '\n';

  std::cout << "----------- Fitting Gaussian Distributions -----------\n";
  hResult1->Fit("GaussDistrib", "QR");
  TF1 *GaussSignal1 = hResult1->GetFunction("GaussDistrib");
  hResult1->Write();
  std::cout << "Fitting first signal (UnlikeCharges):\n"
            << "First signal Fit Chi^2/NDF: " << GaussSignal1->GetChisquare()
            << " / " << GaussSignal1->GetNDF() << '\n'
            << "First signal Fit Probability: " << GaussSignal1->GetProb()
            << '\n'
            << "First signal Fit Par0: " << GaussSignal1->GetParameter(0)
            << " +/- " << GaussSignal1->GetParError(0) << '\n'
            << "First signal Fit Par1: " << GaussSignal1->GetParameter(1)
            << " +/- " << GaussSignal1->GetParError(1) << '\n'
            << "First signal Fit Par2: " << GaussSignal1->GetParameter(2)
            << " +/- " << GaussSignal1->GetParError(2) << '\n';
  hResult2->Fit("GaussDistrib", "QR");
  TF1 *GaussSignal2 = hResult2->GetFunction("GaussDistrib");
  hResult2->Write();
  std::cout << "Fitting second signal (UnlikeKPCharges):\n"
            << "Second signal Fit Chi^2/NDF: " << GaussSignal2->GetChisquare()
            << " / " << GaussSignal2->GetNDF() << '\n'
            << "Second signal Fit Probability: " << GaussSignal2->GetProb()
            << '\n'
            << "Second signal Fit Par0: " << GaussSignal2->GetParameter(0)
            << " +/- " << GaussSignal2->GetParError(0) << '\n'
            << "Second signal Fit Par1: " << GaussSignal2->GetParameter(1)
            << " +/- " << GaussSignal2->GetParError(1) << '\n'
            << "Second signal Fit Par2: " << GaussSignal2->GetParameter(2)
            << " +/- " << GaussSignal2->GetParError(2) << '\n';
  hDaughtersInvMass->Fit("GaussDistrib", "QR");
  TF1 *GaussGeneration = hDaughtersInvMass->GetFunction("GaussDistrib");
  hDaughtersInvMass->Write();
  std::cout << "K* data obtained by decay:\n"
            << "K* mass: " << GaussGeneration->GetParameter(1) << " +/- "
            << GaussGeneration->GetParError(1) << '\n'
            << "K* width: " << GaussGeneration->GetParameter(2) << " +/- "
            << GaussGeneration->GetParError(2) << '\n'
            << "Signal amplitude: " << GaussGeneration->GetParameter(0)
            << " +/- " << GaussGeneration->GetParError(0) << '\n'
            << "Signal Fit Chi^2/NDF: " << GaussGeneration->GetChisquare()
            << " / " << GaussGeneration->GetNDF() << '\n'
            << "Signal Fit Probability: " << GaussGeneration->GetProb() << '\n';

  hTheta->AddDirectory(kFALSE);
  hPhi->AddDirectory(kFALSE);
  hImpulse->AddDirectory(kFALSE);
  hTransvImpulse->AddDirectory(kFALSE);
  hEnergy->AddDirectory(kFALSE);
  hIndex->AddDirectory(kFALSE);
  hTotalInvMass->AddDirectory(kFALSE);
  hUnlikeInvMass->AddDirectory(kFALSE);
  hLikeInvMass->AddDirectory(kFALSE);
  hUnlikeKPInvMass->AddDirectory(kFALSE);
  hLikeKPInvMass->AddDirectory(kFALSE);
  hDaughtersInvMass->AddDirectory(kFALSE);
  hResult1->AddDirectory(kFALSE);
  hResult2->AddDirectory(kFALSE);
  // If the histogram or graph is made persistent, the list of associated
  // functions (e.g. fit functions) is also persistent

  KDecay->Close();
  return 0;
}

void DrawHistos(char const *filename) {
  SetStyle();

  TFile *KDecay = new TFile(filename);
  // open KDecay.root in read mode by default

  TH1F *hIndex = (TH1F *)KDecay->Get("ParticleIndex");
  TH1F *hTheta = (TH1F *)KDecay->Get("ThetaAngle");
  TH1F *hPhi = (TH1F *)KDecay->Get("PhiAngle");
  TH1F *hImpulse = (TH1F *)KDecay->Get("ParticleImpulse");
  TH1F *hTransvImpulse = (TH1F *)KDecay->Get("ParticleTransvImpulse");
  TH1F *hEnergy = (TH1F *)KDecay->Get("ParticleEnergy");
  TH1F *hTotalInvMass = (TH1F *)KDecay->Get("TotalInvMass");
  TH1F *hUnlikeInvMass = (TH1F *)KDecay->Get("UnlikeInvMass");
  TH1F *hLikeInvMass = (TH1F *)KDecay->Get("LikeInvMass");
  TH1F *hUnlikeKPInvMass = (TH1F *)KDecay->Get("UnlikeKPInvMass");
  TH1F *hLikeKPInvMass = (TH1F *)KDecay->Get("LikeKPInvMass");
  TH1F *hDaughtersInvMass = (TH1F *)KDecay->Get("DaughtersInvMass");
  TH1F *hResult1 = (TH1F *)KDecay->Get("Result1");
  TH1F *hResult2 = (TH1F *)KDecay->Get("Result2");

  hTheta->AddDirectory(kFALSE);
  hPhi->AddDirectory(kFALSE);
  hImpulse->AddDirectory(kFALSE);
  hEnergy->AddDirectory(kFALSE);
  hIndex->AddDirectory(kFALSE);
  hTotalInvMass->AddDirectory(kFALSE);
  hUnlikeInvMass->AddDirectory(kFALSE);
  hLikeInvMass->AddDirectory(kFALSE);
  hUnlikeKPInvMass->AddDirectory(kFALSE);
  hLikeKPInvMass->AddDirectory(kFALSE);
  hDaughtersInvMass->AddDirectory(kFALSE);
  hResult1->AddDirectory(kFALSE);
  hResult2->AddDirectory(kFALSE);

  TCanvas *c0 = new TCanvas("canvas0", "Particle and impulse generation data",
                            200, 10, 1000, 600);
  c0->Divide(2, 2);
  c0->cd(1);
  hIndex->GetXaxis()->SetTitle("Type of Particle Generated");
  hIndex->GetYaxis()->SetTitle("Counts");
  hIndex->Draw();
  c0->cd(2);
  hImpulse->GetXaxis()->SetTitle("Impulse (GeV)");
  hImpulse->GetYaxis()->SetTitle("Counts");
  hImpulse->Draw();
  TLegend *leg_exp = new TLegend(0.1, 0.7, 0.48, 0.9, "Impulse distribution");
  leg_exp->AddEntry(hImpulse, "Data");
  leg_exp->AddEntry(hImpulse->GetFunction("ExponDistrib"), "Exponential fit");
  leg_exp->Draw("same");
  c0->cd(3);
  hTheta->GetXaxis()->SetTitle("Theta angle");
  hTheta->GetYaxis()->SetTitle("Counts");
  hTheta->Draw();
  TLegend *leg_theta =
      new TLegend(0.1, 0.7, 0.48, 0.9, "Theta angle distribution");
  leg_theta->AddEntry(hTheta, "Data");
  leg_theta->AddEntry(hTheta->GetFunction("UnifDistrib"), "Uniform Fit");
  leg_theta->Draw("same");
  c0->cd(4);
  hPhi->GetXaxis()->SetTitle("Phi angle");
  hPhi->GetYaxis()->SetTitle("Counts");
  hPhi->Draw();
  TLegend *leg_phi = new TLegend(0.1, 0.7, 0.48, 0.9, "Phi angle distribution");
  leg_phi->AddEntry(hPhi, "Data");
  leg_phi->AddEntry(hPhi->GetFunction("UnifDistrib"), "Uniform Fit");
  leg_phi->Draw("same");

  TCanvas *c1 = new TCanvas("canvas1", "K* signal", 200, 10, 1000, 600);
  c1->Divide(3, 1);
  c1->cd(1);
  hDaughtersInvMass->GetXaxis()->SetRangeUser(0.5, 1.5);
  hDaughtersInvMass->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hDaughtersInvMass->GetYaxis()->SetTitle("Counts");
  hDaughtersInvMass->Draw();
  TLegend *leg_gauss0 =
      new TLegend(0.1, 0.8, 0.3, 0.9, "True K* Invariant Mass Distribution");
  leg_gauss0->AddEntry(hDaughtersInvMass, "Data");
  leg_gauss0->AddEntry(hDaughtersInvMass->GetFunction("GaussDistrib"),
                       "Gaussian Fit");
  leg_gauss0->Draw("same");
  c1->cd(2);
  hResult1->GetXaxis()->SetRangeUser(0.5, 1.5);
  hResult1->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hResult1->SetMinimum(0);
  hResult1->GetYaxis()->SetTitle("Counts");
  hResult1->Draw("BAR");
  TLegend *leg_gauss1 = new TLegend(0.1, 0.8, 0.3, 0.9,
                                    "K* Invariant Mass with opposite charges");
  leg_gauss1->AddEntry(hResult1, "Data");
  leg_gauss1->AddEntry(hResult1->GetFunction("GaussDistrib"), "Gaussian Fit");
  leg_gauss1->Draw("same");
  c1->cd(3);
  hResult2->GetXaxis()->SetRangeUser(0.5, 1.5);
  hResult2->GetXaxis()->SetTitle("Invariant Mass (GeV/c^2)");
  hResult2->SetMinimum(0);
  hResult2->GetYaxis()->SetTitle("Counts");
  hResult2->Draw("BAR");
  TLegend *leg_gauss2 =
      new TLegend(0.1, 0.8, 0.3, 0.9, "K* Invariant Mass with K/P pair");
  leg_gauss2->AddEntry(hResult2, "Data");
  leg_gauss2->AddEntry(hResult2->GetFunction("GaussDistrib"), "Gaussian Fit");
  leg_gauss2->Draw("same");

  c0->UseCurrentStyle();
  c1->UseCurrentStyle();

  c0->Print("GenerationData.pdf");
  c0->Print("GenerationData.root");
  c0->Print("GenerationData.C");
  c1->Print("K*Decay.pdf");
  c1->Print("K*Decay.root");
  c1->Print("K*Decay.C");

  KDecay->Close();

  return;
}