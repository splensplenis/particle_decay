#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"
#include <algorithm>
#include <cmath>   //M_PI
#include <cstdlib> //RAND_MAX
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<ParticleType *> Particle::fTable{};

Particle::Particle(std::string particle_name, double px, double py, double pz)
    : fPx{px}, fPy{py}, fPz{pz} {
  fIndex = FindParticle(particle_name);
  if (fIndex == -1) {
    throw std::invalid_argument{
        "Particle name does not correspond to any knwon Particle Type"};
  }
}
Particle::Particle() : fPx{}, fPy{}, fPz{} { fIndex = -1; }

int Particle::FindParticle(std::string particle_name) {
  for (size_t i = 0.; i != fTable.size(); ++i) {
    if (particle_name == fTable[i]->GetName()) {
      return i;
    }
  }
  return -1;
}

int Particle::GetIndex() const { return fIndex; }
void Particle::SetIndex(std::string name) {
  if (FindParticle(name) != -1) {
    fIndex = FindParticle(name);
  } else {
    throw std::invalid_argument{
        "Particle name does not correspond to any known Particle Type"};
  }
}
void Particle::SetIndex(int index) {
  if (index > -1 && index < fTable.size()) {
    fIndex = index;
  } else {
    throw std::invalid_argument{
        "Particle name does not correspond to any known Particle Type"};
  }
}

void Particle::AddToTable(std::string particle_name, double mass, int charge,
                          double width) {
  if (FindParticle(particle_name) == -1) {
    ParticleType *addentry =
        new ResonanceType(particle_name, mass, charge, width);
    fTable.push_back(addentry);
  } else {
    // particle already present, do nothing
  }
}

void Particle::PrintTable() {
  for (size_t i = 0; i != fTable.size(); ++i) {
    std::cout << "------\n";
    fTable[i]->Print();
  }
}
void Particle::PrintParticle() const {
  std::cout << "Particle has index " << GetIndex() << '\n'
            << "and is of type " << GetName() << '\n';
  std::cout << "Px is " << GetParticlePx() << '\n'
            << "Py is " << GetParticlePy() << '\n'
            << "Pz is " << GetParticlePz() << '\n';
}
std::string Particle::GetName() const { return fTable[GetIndex()]->GetName(); }
double Particle::GetMass() const { return fTable[GetIndex()]->GetMass(); }
int Particle::GetCharge() const { return fTable[GetIndex()]->GetCharge(); }
double Particle::GetParticlePx() const { return fPx; }
double Particle::GetParticlePy() const { return fPy; }
double Particle::GetParticlePz() const { return fPz; }
double Particle::GetParticleImpulse() const {
  return sqrt(GetParticlePx() * GetParticlePx() +
              GetParticlePy() * GetParticlePy() +
              GetParticlePz() * GetParticlePz());
}
void Particle::SetParticleImpulse(double px, double py, double pz) {
  fPx = px;
  fPy = py;
  fPz = pz;
}
double Particle::GetParticleEnergy() const {
  return sqrt(GetMass() * GetMass() + GetParticleImpulse());
}
double Particle::InvariantMass(Particle &other) const {
  double delta_E = GetParticleEnergy() + other.GetParticleEnergy();
  double delta_px = (GetParticlePx() + other.GetParticlePx());
  double delta_py = (GetParticlePy() + other.GetParticlePy());
  double delta_pz = (GetParticlePz() + other.GetParticlePz());
  return sqrt(delta_E * delta_E - delta_px * delta_px - delta_py * delta_py -
              delta_pz * delta_pz);
}

void Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (GetMass() == 0.0) {
    throw std::invalid_argument{
        "Decayment cannot be preformed if mass is zero"};
  }

  double massMot = GetMass();       // mass Mother particle
  double massDau1 = dau1.GetMass(); // mass daughter particle
  double massDau2 = dau2.GetMass();

  if (GetIndex() > -1) {
    // gaussian random numbers

    float x1, x2, w, y1, y2;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fTable[GetIndex()]->GetWidth() * y1; // add width effect
  }

  if (massMot < massDau1 + massDau2) {
    throw std::runtime_error{
        "Decayment cannot be preformed because mass is too low in "
        "this channel"};
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetParticleImpulse(pout * sin(theta) * cos(phi),
                          pout * sin(theta) * sin(phi), pout * cos(theta));
  dau2.SetParticleImpulse(-pout * sin(theta) * cos(phi),
                          -pout * sin(theta) * sin(phi), -pout * cos(theta));

  double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);
}

void Particle::Boost(double bx, double by, double bz) {

  double energy = GetParticleEnergy(); // boost this Lorentz vector

  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}