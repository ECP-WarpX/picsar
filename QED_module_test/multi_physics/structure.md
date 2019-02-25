#PICSAR multi-physics library

##Physical modules

###"Simple" replacements of the Boris pusher:
- Classical Landau-Lifshitz Radiation Reaction
- Modified Landau-Lifshitz Radiation Reaction

###Processes invovling leptons and photons
- Nonlinear Breit-Wheeler
- Quantum synchrotron emission

###Processes involving ions, leptons and photons
- Bremsstrahlung
- Bethe-Heitler
- Trident

###Processes requiring to modify Maxwell solver
- Vacuum polarization

##Technical requirements
Must comply with strict C++11. Could optionally use features from C++17 STL
(e.g. special functions). Will have a Fortran interface. Will use cmake as
building system. Will be incapsulated in the namespace `picsar::multiphysics` .

##Dependencies
- Boost library (will be optional)

##Normalized units convention
The co

---

##Modules in detail

###Classical Landau-Lifshitz Radiation Reaction

This module should provide a replacement for the Boris pusher algorithm. It
implements classical radiation reaction according to the algorithm proposed
by Matteo Tamburini. Basically, Radiation Reaction is modeled as an additional force
term acting on each particle. The only additional requirement with respect to the
"regular" Boris pusher is the lengthscale λ.

####Implementation
This module can be implemented as a collection of functions, defined in a header file and
implemented in a cpp file.
~~~~
void boris_plus_landau_lifshitz_push(
  std::vector<double>& px, std::vector<double>& py, std::vector<double>& pz,
  const std::vector<double>& ex, const std::vector<double>& ey, const std::vector<double>& ez,
  const std::vector<double>& bx, const std::vector<double>& by, const std::vector<double>& bz,
  double mass, double charge, double dt, double lambda);
~~~~


####References
- M.Tamburini et al. "Radiation reaction effects on radiation pressure acceleration".
New Jouranl of Physics 12, 123005 (2010)
[https://iopscience.iop.org/article/10.1088/1367-2630/12/12/123005/meta]
- M.Tamburini. "Radiation reaction effects in superintense laser-plasma interaction".
PhD thesis, Università di Pisa, Italy 2011.
[https://etd.adm.unipi.it/t/etd-11232011-111059/]
- M.Vranic et al. "Classical radiation reaction in particle-in-cell simulations".
Computer Physics Communications 204, 141-151 (2016)
[https://www.sciencedirect.com/science/article/pii/S001046551630090X]
