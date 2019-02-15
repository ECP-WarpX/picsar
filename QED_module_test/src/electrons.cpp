#include "electrons.h"

using namespace testbed;

electrons::electrons(std::string name):leptons(name){
    header = std::string("electrons");
    init_mass_and_charge();
}

electrons::~electrons(){}

void electrons::init_mass_and_charge(){
    mass = 1.0;
    charge = - 1.0;
}
