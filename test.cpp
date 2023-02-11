//compile with: 
//g++ -Wall -Wextra -o test particle_type.cpp resonance_type.cpp particle.cpp test.cpp

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "particle.hpp"
#include "particle_type.hpp"
#include "resonance_type.hpp"
#include <iostream>
#include <vector>

TEST_CASE("Testing Particles") {
  // constructors test: cannot instantiate a ParticleType (abstract class)
  ResonanceType *p1 = new ResonanceType("M", 20, 1, 0);
  ResonanceType *p2 = new ResonanceType("N", 35, 2, 0.4);
  std::cout << "Testing Particle ('N', 35, 2, 0.4):\n";
  p2->Print();

  // run-time polimorphism test
  std::cout << "Testing both particles with Print method:\n";
  std::vector<ParticleType *> vecpar{p1, p2};
  for (size_t i = 0; i != vecpar.size(); ++i) {
    vecpar[i]->Print();
  }

  delete p1;
  delete p2;
}

TEST_CASE("Testing Table") {
  Particle::AddToTable("A", 1.5, 1);
  Particle::AddToTable("B", 2.4, -2);
  Particle::AddToTable("C", 3.6, 3);
  Particle::AddToTable("R", 0.01, -1, 0.1);

  // Particle::PrintTable();

  Particle par0 = Particle("C", 0.5, -0.5, 1);
  CHECK(par0.GetIndex() == 2);
  CHECK(par0.GetParticlePx() == 0.5);
  CHECK(par0.GetParticlePy() == -0.5);
  CHECK(par0.GetParticlePz() == 1);
  CHECK(par0.GetParticleEnergy() == doctest::Approx(3.8).epsilon(0.01));

  Particle par1 = Particle("B");
  CHECK(par1.GetParticlePx() == 0);
  CHECK(par1.GetParticlePy() == 0);
  CHECK(par1.GetParticlePz() == 0);
  par1.SetParticleImpulse(1, -2, 0.01);
  CHECK(par1.GetParticlePx() == 1);
  CHECK(par1.GetParticlePy() == -2);
  CHECK(par1.GetParticlePz() == 0.01);
  par1.PrintParticle();
  CHECK(par1.InvariantMass(par0) == doctest::Approx(6.37).epsilon(0.01));
  CHECK(par0.InvariantMass(par1) == doctest::Approx(6.37).epsilon(0.01));

  Particle res0 = Particle("R");
  res0.PrintParticle();

  // testing for a particle unknown to the Table
  CHECK_THROWS(Particle("X", 0, 0, 1));

  //testing for index exceptions
  res0.SetIndex("A");
  CHECK(res0.GetIndex() == 0);
  CHECK_THROWS(res0.SetIndex(4));
  CHECK_THROWS(res0.SetIndex(11));
  CHECK_THROWS(res0.SetIndex("X"));

  SUBCASE("Testing Decay errors") {
    Particle::AddToTable("D0", 0, 1, 0);
    Particle::AddToTable("D1", 1, 2, 0);

    Particle mot0 = Particle("D0", 0.5, 0.6, 1);
    Particle mot1 = Particle("D1", 0.5, 0.6, 1);
    Particle dau1 = Particle("A");
    Particle dau2 = Particle("B");
    CHECK_THROWS(mot0.Decay2body(dau1, dau2)); // mot0 has no mass
    CHECK_THROWS(mot1.Decay2body(dau1, dau2));
    // mot1 mass is lower than dau1+dau2 mass
    //(width is 0 so no width effect can raise mot1 mass value)
    //so decay cannot happen
  }
}
