# particle_decay
Simulazione di urti tra particelle e decadimento di una risonanza, analisi dati.

// to compile and execute, on root command line:
// root [0] gROOT->LoadMacro("particle_type.cpp+")
// root [1] gROOT->LoadMacro("resonance_type.cpp+")
// root [2] gROOT->LoadMacro("particle.cpp+")
// root [3] gROOT->LoadMacro("functions.cpp+")
// root [4] .L main.cpp+
// root [5] KDecay()
// note: launch commands 0-3 just once to create shared libraries

Risultati sulla generazione casuale (metodo Montecarlo):
[GenerationData.pdf](https://github.com/splensplenis/particle_decay/files/10714216/GenerationData.pdf)

Risultati sull'analisi dati dei decadimenti:
[KDecay.pdf](https://github.com/splensplenis/particle_decay/files/10714218/K.Decay.pdf)

(root software is necessary)
