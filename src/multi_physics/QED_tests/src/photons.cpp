#include "photons.h"

using namespace testbed;

photons::photons(std::string name):species(name){
    header = std::string("photons");
    mass = 0.0;
    charge = 0.0;
}

photons::~photons(){}

void photons::push_positions(ttime dt){
    int num_particles = pos[0].size();

    for (int i = 0; i < num_particles; i++){
        double inv_ptotal = 1.0/sqrt(mom[0][i]*mom[0][i] + mom[1][i]*mom[1][i] + mom[2][i]*mom[2][i]);
        double coeff = dt*inv_ptotal;

        pos[0][i] += coeff*mom[0][i];
        pos[1][i] += coeff*mom[1][i];
        pos[2][i] += coeff*mom[2][i];
    }
}
void photons::push_momenta(ttime){
    return;
}
