#include "resonance_type.hpp"
#include <iostream>

double ResonanceType::GetWidth() const { return fWidth; }
void ResonanceType::Print() const {
  ParticleType::Print();
  if (fWidth > 0) {
    std::cout << "Particle (resonance) width is " << fWidth << " GeV/c^2\n";
  }
}