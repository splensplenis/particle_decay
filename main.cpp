// to compile and execute, on root command line:
// root [0] gROOT->LoadMacro("particle_type.cpp+")
// root [1] gROOT->LoadMacro("resonance_type.cpp+")
// root [2] gROOT->LoadMacro("particle.cpp+")
// root [3] gROOT->LoadMacro("functions.cpp+")
// root [4] .L main.cpp+
// root [5] KDecay()
// note: launch commands 0-3 just once to create shared libraries

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"

#include <iostream>
#include <stdexcept>
#include <vector>

#include "functions.hpp"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"

int KDecay(long int seed = 0) {

  R__LOAD_LIBRARY(particle_type_cpp.so);
  R__LOAD_LIBRARY(resonance_type_cpp.so);
  R__LOAD_LIBRARY(particle_cpp.so);
  R__LOAD_LIBRARY(functions_cpp.so);

  if (seed == 0.) {
    gRandom->SetSeed();
    std::cout << "Random generation seed number is " << gRandom->GetSeed()
              << '\n';
  } else {
    gRandom->SetSeed(seed);
    std::cout << "Generating according to seed number " << seed << '\n';
  }

  TH1F *hIndex =
      new TH1F("ParticleIndex", "Generated particles histo", 7, 0., 7.);
  TH1F *hTheta = new TH1F("ThetaAngle", "Theta histo", 1000, 0., TMath::Pi());
  TH1F *hPhi = new TH1F("PhiAngle", "Phi histo", 1000, 0., 2 * TMath::Pi());
  TH1F *hImpulse = new TH1F("ParticleImpulse", "Impulse histo", 1E3, 0., 4.);
  TH1F *hTransvImpulse = new TH1F("ParticleTransvImpulse",
                                  "Transverse Impulse histo", 1E3, 0., 4.);
  TH1F *hEnergy = new TH1F("ParticleEnergy", "Energy histo", 1E3, 0., 4.);

  TH1F *hTotalInvMass =
      new TH1F("TotalInvMass", "Total Invariant Mass histo", 500, 0., 4.);
  TH1F *hUnlikeInvMass = new TH1F(
      "UnlikeInvMass", "Opposite charges Invariant Mass histo", 500, 0., 4.);
  TH1F *hLikeInvMass =
      new TH1F("LikeInvMass", "Same charges Invariant Mass histo", 500, 0., 4.);
  TH1F *hUnlikeKPInvMass =
      new TH1F("UnlikeKPInvMass",
               "Unlike charge Kaon/Pion Pair Invariant Mass histo", 500, 0, 4.);
  TH1F *hLikeKPInvMass =
      new TH1F("LikeKPInvMass",
               "Like charge Kaon/Pion Pair Invariant Mass histo", 500, 0., 4.);
  TH1F *hDaughtersInvMass = new TH1F(
      "DaughtersInvMass", "Decay daughters Invariant Mass histo", 500, 0., 4.);

  try {
    struct DecayPair {
      Particle pair1{};
      Particle pair2{};
    };
    struct GenerationData{
      Double_t theta{};
      Double_t phi{};
      Double_t impulse{};
    };

    Particle::AddToTable("K*", 0.89166, 0, 0.050); // fIndex == 0
    Particle::AddToTable("Pi+", 0.13957, +1);      // fIndex == 1
    Particle::AddToTable("Pi-", 0.13957, -1);      // fIndex == 2
    Particle::AddToTable("K+", 0.49367, +1);       // fIndex == 3
    Particle::AddToTable("K-", 0.49367, -1);       // fIndex == 4
    Particle::AddToTable("P+", 0.93827, +1);       // fIndex == 5
    Particle::AddToTable("P-", 0.93827, -1);       // fIndex == 6
    // mass and width unit is GeV/c^2, charge is measured with electron charge

    std::vector<Particle> EventParticles{};
    std::vector<DecayPair> DecayParticles{};
    std::vector<GenerationData> EventData{};

    Double_t theta, phi = 0;
    Double_t px, py, pz = 0;
    Double_t mean_impulse = 1;

    Int_t nEventGen = 1E5;
    Int_t nPartGen = 100;

    for (Int_t i = 0; i < nEventGen; ++i) {
      EventParticles.clear();
      DecayParticles.clear();
      EventData.clear();

      // random particle generation and K*decay
      for (Int_t j = 0; j < nPartGen; ++j) {

        theta = gRandom->Rndm() * TMath::Pi();
        phi = gRandom->Rndm() * 2 * TMath::Pi();
        Double_t impulse = gRandom->Exp(mean_impulse);
        EventData.push_back(GenerationData{theta, phi, impulse});
        px = impulse * TMath::Sin(theta) * TMath::Cos(phi);
        py = impulse * TMath::Sin(theta) * TMath::Sin(phi);
        pz = impulse * TMath::Cos(theta);
        Particle particle{};
        particle.SetParticleImpulse(px, py, pz);

        Double_t probability = gRandom->Rndm();
        if (probability < 0.4) {
          particle.SetIndex("Pi+");
        } else if (probability < 0.8) {
          particle.SetIndex("Pi-");
        } else if (probability < 0.85) {
          particle.SetIndex("K+");
        } else if (probability < 0.9) {
          particle.SetIndex("K-");
        } else if (probability < 0.945) {
          particle.SetIndex("P+");
        } else if (probability < 0.99) {
          particle.SetIndex("P-");
        } else {
          particle.SetIndex("K*");
          Double_t occurrency = gRandom->Rndm();
          Particle decaydaug1{}, decaydaug2{};
          if (occurrency < 0.5) {
            decaydaug1.SetIndex("Pi+");
            decaydaug2.SetIndex("K-");
          } else {
            decaydaug1.SetIndex("Pi-");
            decaydaug2.SetIndex("K+");
          }
          particle.Decay2body(decaydaug1, decaydaug2);
          DecayPair pair{decaydaug1, decaydaug2};
          DecayParticles.push_back(pair);
          EventParticles.push_back(decaydaug1);
          EventParticles.push_back(decaydaug2);
        }

        EventParticles.push_back(particle);
      }

      // filling histograms
      for (auto data : EventData) {
        hTheta->Fill(data.theta);
        hPhi->Fill(data.phi);
        hImpulse->Fill(data.impulse);
      }
      for (auto particle : EventParticles) {
        hEnergy->Fill(particle.GetParticleEnergy());
        hIndex->Fill(particle.GetIndex());
        hTransvImpulse->Fill(
            sqrt(particle.GetParticlePx() * particle.GetParticlePx() +
                 particle.GetParticlePy() * particle.GetParticlePy()));
      }
      for (auto pair : DecayParticles) {
        hDaughtersInvMass->Fill((pair.pair1).InvariantMass(pair.pair2));
      }

      for (size_t j = 0; j < EventParticles.size(); ++j) {
        for (size_t k = j + 1; k < EventParticles.size(); ++k) {
          Particle &first = EventParticles[j];
          Particle &second = EventParticles[k];
          hTotalInvMass->Fill(first.InvariantMass(second));
          if (first.GetCharge() * second.GetCharge() < 0) {
            hUnlikeInvMass->Fill(first.InvariantMass(second));
          }
          if (first.GetCharge() * second.GetCharge() > 0) {
            hLikeInvMass->Fill(first.InvariantMass(second));
          }
          if ((first.GetName() == "K+" && second.GetName() == "Pi-") ||
              (first.GetName() == "K-" && second.GetName() == "Pi+") ||
              (first.GetName() == "Pi-" && second.GetName() == "K+") ||
              (first.GetName() == "Pi+" && second.GetName() == "K-")) {
            hUnlikeKPInvMass->Fill(first.InvariantMass(second));
          }
          if ((first.GetName() == "K+" && second.GetName() == "Pi+") ||
              (first.GetName() == "K-" && second.GetName() == "Pi-") ||
              (first.GetName() == "Pi+" && second.GetName() == "K+") ||
              (first.GetName() == "Pi-" && second.GetName() == "K-")) {
            hLikeKPInvMass->Fill(first.InvariantMass(second));
          }
        }
      }
    }
  } catch (std::invalid_argument &ia) {
    std::cerr << "Invalid argument error: " << ia.what() << '\n';
  } catch (std::runtime_error &re) {
    std::cerr << "Runtime error: " << re.what() << '\n';
  }

  TFile *KDecay = new TFile("KDecay.root", "RECREATE");
  hIndex->Write();
  hTheta->Write();
  hPhi->Write();
  hImpulse->Write();
  hTransvImpulse->Write();
  hEnergy->Write();
  hTotalInvMass->Write();
  hUnlikeInvMass->Write();
  hLikeInvMass->Write();
  hUnlikeKPInvMass->Write();
  hLikeKPInvMass->Write();
  hDaughtersInvMass->Write();

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

  KDecay->Close();

  DataAnalysis("KDecay.root");
  DrawHistos("KDecay.root");

  return 0;
}