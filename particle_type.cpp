#include "particle_type.hpp"
#include <iostream>
#include <stdexcept>

ParticleType::ParticleType(std::string name, double mass, int charge)
    : fName{name}, fMass{mass}, fCharge{charge} {
  if (mass < 0) {
    throw std::invalid_argument{"Mass of a particle cannot be negative"};
  }
}

const std::string ParticleType::GetName() const { return fName; }
double ParticleType::GetMass() const { return fMass; }
int ParticleType::GetCharge() const { return fCharge; }

void ParticleType::Print() const {
  std::cout << "Particle name is " << fName << '\n'
            << "Particle mass is " << fMass << " GeV/c^2\n"
            << "Particle charge is " << fCharge << " e\n";
}