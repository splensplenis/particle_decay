#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "particle_type.hpp"
#include "resonance_type.hpp"
#include <vector>
#include <string>

class Particle {
public:
  Particle(std::string particle_name, double px = 0., double py = 0., double pz = 0.);
  Particle();

  int GetIndex() const;
  void SetIndex(std::string particletype_name);
  void SetIndex(int code_number);

  static void AddToTable(std::string particle_name, double mass, int charge,
                         double width = 0.);
  static void PrintTable();

  void PrintParticle() const;
  std::string GetName() const;
  double GetMass() const;
  int GetCharge() const;

  double GetParticlePx() const;
  double GetParticlePy() const;
  double GetParticlePz() const;
  double GetParticleImpulse() const;
  void SetParticleImpulse(double, double, double);
  double GetParticleEnergy() const;

  double InvariantMass(Particle &other) const;

  void Decay2body(Particle &dau1, Particle &dau2) const;

private:
  static std::vector<ParticleType *> fTable;
  int fIndex;
  double fPx;
  double fPy;
  double fPz;

  static int FindParticle(std::string particle_name);
  void Boost(double bx, double by, double bz);
};

#endif