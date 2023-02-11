#ifndef RESONANCE_TYPE
#define RESONANCE_TYPE

#include "particle_type.hpp"
#include <string>

class ResonanceType : public ParticleType {
  public:
    ResonanceType(std::string name, double mass, double charge, double width) 
    : ParticleType(name, mass, charge), fWidth{width} {}
    //~ResonanceType() = default;
    
    double GetWidth() const;

    void Print() const;
  private:
    double fWidth; 
};

#endif